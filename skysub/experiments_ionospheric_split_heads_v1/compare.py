"""Ionospheric coefficient scores for the split-heads notebook: deployed, control, split heads, dedicated MLPs."""

from pathlib import Path

import numpy as np
import pandas as pd

from mlp_predictor import serialization, trainer

LINES = ('ATOM_Or', 'ATOM_Orc_OI0777', 'ATOM_Orc_OI0845', 'ATOM_N', 'ATOM_Og')
DEDICATED = {'ATOM_Or': 'oi_mlp_or', 'ATOM_Orc_OI0777': 'oi_mlp_orc', 'ATOM_Orc_OI0845': 'oi_mlp_orc',
             'ATOM_Og': 'oi_mlp_og'}
FLOOR = 1e-4  # internal 1e14 units, as in experiments_oi_dedicated_mlp_v1


def _columns(triplet):
    names = [str(n) for n in triplet['coef_names']]
    return [names.index(n) for n in LINES]


def predict_lines(ensembles, triplet, rows):
    """Ensemble-mean science coefficients of LINES on ``rows`` for each saved ensemble that exists."""
    cols, out = _columns(triplet), {}
    have = [str(n) for n in triplet['ctx_names']]
    for label, path in ensembles.items():
        if not Path(path).exists():
            print(f'{label}: {path} not found; skipped')
            continue
        ens = serialization.load_ensemble(path, device='cpu')
        ctx = [have.index(n) for n in ens['ctx_names']]
        arm = {a: np.asarray(triplet[f'coef_{a}'], np.float32)[rows] for a in ('near', 'far')}
        cx = {a: np.asarray(triplet[f'ctx_{a}'], np.float32)[rows][:, ctx] for a in ('near', 'far', 'sci')}
        pred = trainer.predict_sci_coefficients_default(
            ens, coef_near_phys=arm['near'], coef_far_phys=arm['far'],
            ctx_near_phys=cx['near'], ctx_far_phys=cx['far'], ctx_sci_phys=cx['sci'])
        out[label] = np.asarray(pred, np.float64)[:, cols]
    return out


def _oi_index(oi_dir, triplet, rows):
    table = np.load(Path(oi_dir) / 'oi_table.npz')
    pos = {int(r): i for i, r in enumerate(table['row_index'])}
    return np.array([pos.get(int(r), -1) for r in np.asarray(triplet['row_index'])[rows]])


def dedicated_predictions(oi_dir, triplet, rows):
    """The dedicated per-mechanism networks of experiments_oi_dedicated_mlp_v1, aligned by corpus row."""
    idx = _oi_index(oi_dir, triplet, rows)
    pred = np.load(Path(oi_dir) / 'run_v1' / 'predictions.npz')
    oi_names = [str(n) for n in pred['oi_names']]
    out = np.full((len(rows), len(LINES)), np.nan)
    for j, line in enumerate(LINES):
        if line in DEDICATED:
            col = 10 ** pred[DEDICATED[line]][:, oi_names.index(line)]
            out[idx >= 0, j] = col[idx[idx >= 0]]
    print(f'dedicated: {np.sum(idx >= 0)}/{len(rows)} test rows found in the [O I] experiment table')
    return {'dedicated': out}


def clean_subset(oi_dir, triplet, rows):
    """Rows the [O I] experiment kept: science [O I] within 0.3 dex of the sky arms, tile 11111 excluded."""
    idx = _oi_index(oi_dir, triplet, rows)
    clean = np.load(Path(oi_dir) / 'run_v1' / 'predictions.npz')['clean']
    return (idx >= 0) & clean[np.maximum(idx, 0)]


def _errors(preds, triplet, rows):
    truth = np.asarray(triplet['coef_sci'], np.float64)[rows][:, _columns(triplet)]
    log_true = np.log10(np.maximum(truth, FLOOR))
    return {k: np.log10(np.maximum(v, FLOOR)) - log_true for k, v in preds.items()}, truth > 0


def score_tables(preds, triplet, rows, clean):
    """Median |log10 error| per method and line, on all test rows and on the clean subset."""
    err, positive = _errors(preds, triplet, rows)
    tables = {}
    for label, keep in (('all test rows', np.ones(len(rows), bool)), ('clean test rows', clean)):
        tables[f'{label} (n={keep.sum()})'] = pd.DataFrame(
            {line: {k: np.nanmedian(np.abs(e[keep & positive[:, j], j])) for k, e in err.items()}
             for j, line in enumerate(LINES)})
    return tables


def bootstrap_table(preds, triplet, rows, clean, pairs, n_draw=2000, seed=0):
    """Relative MAD change of method a against b per line, with a 95% interval from resampling test nights."""
    err, positive = _errors(preds, triplet, rows)
    night = np.floor(np.asarray(triplet['obstime_mjd'], np.float64)[rows] - 0.5).astype(int)
    rng = np.random.default_rng(seed)
    records = []
    for subset, keep in (('all', np.ones(len(rows), bool)), ('clean', clean)):
        nights = np.unique(night[keep])
        draws = [rng.choice(nights, nights.size) for _ in range(n_draw)]
        for a, b in pairs:
            if a not in err or b not in err:
                continue
            for j, line in enumerate(LINES):
                ok = keep & positive[:, j] & np.isfinite(err[a][:, j]) & np.isfinite(err[b][:, j])
                if not ok.any():
                    continue
                ea, eb = np.abs(err[a][:, j]), np.abs(err[b][:, j])
                by_night = {n: np.flatnonzero(ok & (night == n)) for n in nights}
                ratio = []
                for d in draws:
                    s = np.concatenate([by_night[n] for n in d])
                    ratio.append(np.median(ea[s]) / np.median(eb[s]) - 1)
                records.append(dict(subset=subset, line=line, pair=f'{a} vs {b}',
                                    change=np.median(ea[ok]) / np.median(eb[ok]) - 1,
                                    lo=np.percentile(ratio, 2.5), hi=np.percentile(ratio, 97.5),
                                    better_rows=np.mean(ea[ok] < eb[ok])))
    return pd.DataFrame(records).set_index(['subset', 'line', 'pair'])
