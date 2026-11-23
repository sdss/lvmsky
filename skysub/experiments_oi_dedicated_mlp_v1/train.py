"""Dedicated [O I] MLPs against the deployed ensemble, on the table written by ``prepare.py``.

    python -m experiments_oi_dedicated_mlp_v1.train --table FILE --out-dir DIR

Every variant sees the same inputs: log10 of the four O I coefficients in both sky
arms plus the deployed ensemble's 42 context features for the science, near and far
pointings.  Variants differ only in which science coefficients they predict:

    or    [O I] 6300/6364            (ATOM_Or)
    orc   O I 7774 + 8446            (ATOM_Orc_OI0777, ATOM_Orc_OI0845)
    og    [O I] 5577                 (ATOM_Og)
    all   all four in one network

The target is log10(sci) - log10(near), weighted by the science coefficient error in
log space with a 0.02 dex floor.  Rows whose science ATOM_Or differs from the
geometric mean of the two sky arms by more than ``--outlier-dex`` (and tile 11111)
are removed from training, validation and the primary test score.  The seed
networks and the feature scaling are saved to ``oi_models.pt`` for ``exposure.py``.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from torch import nn

VARIANTS = {
    'or': ('ATOM_Or',),
    'orc': ('ATOM_Orc_OI0777', 'ATOM_Orc_OI0845'),
    'og': ('ATOM_Og',),
    'all': ('ATOM_Or', 'ATOM_Orc_OI0777', 'ATOM_Orc_OI0845', 'ATOM_Og'),
}
FLOOR = 1e-4          # coefficient floor before log10 (internal 1e14 units)
SIGMA_FLOOR_DEX = 0.02
EXCLUDED_TILES = (11111,)
EPOCHS = 400


def log10(x):
    return np.log10(np.maximum(x, FLOOR))


def clean_mask(t, outlier_dex):
    i = list(t['oi_names']).index('ATOM_Or')
    sci, near, far = (t[f'coef_{a}'][:, i] for a in ('sci', 'near', 'far'))
    ratio = log10(sci) - 0.5 * (log10(near) + log10(far))
    positive = (sci > 0) & (near > 0) & (far > 0)
    return positive & (np.abs(ratio) <= outlier_dex) & ~np.isin(t['tile_id'], EXCLUDED_TILES)


def features(t):
    x = np.hstack([log10(t['coef_near']), log10(t['coef_far']), t['ctx_sci'], t['ctx_near'], t['ctx_far']])
    return np.nan_to_num(x.astype(np.float64))


class Net(nn.Module):
    def __init__(self, n_in, n_out, hidden=(128, 64), dropout=0.1):
        super().__init__()
        layers, d = [], n_in
        for h in hidden:
            layers += [nn.Linear(d, h), nn.LayerNorm(h), nn.GELU(), nn.Dropout(dropout)]
            d = h
        self.body = nn.Sequential(*layers, nn.Linear(d, n_out))

    def forward(self, x):
        return self.body(x)


def fit_one(x_tr, y_tr, w_tr, x_va, y_va, w_va, seed, epochs=EPOCHS, patience=40, lr=1e-3, wd=1e-4, batch=256):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    net = Net(x_tr.shape[1], y_tr.shape[1])
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=wd)
    tensors = [torch.as_tensor(a, dtype=torch.float32) for a in (x_tr, y_tr, w_tr, x_va, y_va, w_va)]
    xt, yt, wt, xv, yv, wv = tensors
    best, best_state, wait = np.inf, None, 0
    for _ in range(epochs):
        net.train()
        for b in np.array_split(rng.permutation(len(xt)), max(1, len(xt) // batch)):
            opt.zero_grad()
            loss = (wt[b] * (net(xt[b]) - yt[b]) ** 2).mean()
            loss.backward()
            opt.step()
        net.eval()
        with torch.no_grad():
            val = float((wv * (net(xv) - yv) ** 2).mean())
        if val < best - 1e-6:
            best, best_state, wait = val, {k: v.clone() for k, v in net.state_dict().items()}, 0
        else:
            wait += 1
            if wait >= patience:
                break
    net.load_state_dict(best_state)
    return net, best


def scores(err, truth, pred):
    """Log-space error statistics plus the amplitude-weighted relative RMS."""
    a = np.abs(err)
    rel = np.sqrt(np.sum((pred - truth) ** 2) / np.sum(truth ** 2))
    return dict(n=int(err.size), bias_dex=float(np.median(err)), mad_dex=float(np.median(a)),
                rms_dex=float(np.sqrt(np.mean(err ** 2))), p90_dex=float(np.percentile(a, 90)),
                weighted_rel_rms=float(rel))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--table', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--seeds', type=int, default=5)
    p.add_argument('--outlier-dex', type=float, default=0.3)
    a = p.parse_args()

    t = dict(np.load(a.table, allow_pickle=False))
    names = list(t['oi_names'])
    clean = clean_mask(t, a.outlier_dex)
    split = {k: t[f'{k}_idx'] for k in ('train', 'val', 'test')}
    tr, va = (np.intersect1d(split[k], np.flatnonzero(clean)) for k in ('train', 'val'))
    te_all = split['test'][(t['coef_sci'][split['test']] > 0).all(axis=1)]
    te = np.intersect1d(te_all, np.flatnonzero(clean))

    x = features(t)
    med, iqr = np.median(x[tr], 0), np.subtract(*np.percentile(x[tr], [75, 25], 0))
    keep = iqr > 1e-8
    x = np.clip((x[:, keep] - med[keep]) / iqr[keep], -10, 10)

    log_sci, log_near, log_far = (log10(t[f'coef_{k}']) for k in ('sci', 'near', 'far'))
    sigma = t['err_sci'] / (np.maximum(t['coef_sci'], FLOOR) * np.log(10))
    weight = 1.0 / (np.nan_to_num(sigma, nan=1.0) ** 2 + SIGMA_FLOOR_DEX ** 2)
    weight /= np.median(weight[tr], axis=0)
    y = log_sci - log_near

    preds = {'near_copy': log_near, 'sky_geomean': 0.5 * (log_near + log_far),
             'deployed_ensemble': log10(t['deployed_pred'])}
    history, saved = {}, {}
    for variant, lines in VARIANTS.items():
        cols = [names.index(n) for n in lines]
        out = np.zeros((len(x), len(cols)))
        history[variant], saved[variant] = [], dict(lines=list(lines), states=[])
        for seed in range(a.seeds):
            net, best = fit_one(x[tr], y[tr][:, cols], weight[tr][:, cols],
                                x[va], y[va][:, cols], weight[va][:, cols], seed)
            with torch.no_grad():
                out += net(torch.as_tensor(x, dtype=torch.float32)).numpy() / a.seeds
            history[variant].append(best)
            saved[variant]['states'].append(net.state_dict())
        full = np.full(log_sci.shape, np.nan)
        full[:, cols] = out + log_near[:, cols]
        preds[f'oi_mlp_{variant}'] = full

    results = {}
    for label, rows in (('test_clean', te), ('test_all', te_all)):
        results[label] = {}
        for k, name in enumerate(names):
            truth = t['coef_sci'][rows, k]
            results[label][name] = {
                method: scores(pr[rows, k] - log_sci[rows, k], truth, 10 ** pr[rows, k])
                for method, pr in preds.items() if np.isfinite(pr[rows, k]).all()
            }
    summary = dict(rows=dict(total=int(len(x)), clean=int(clean.sum()), train=int(tr.size), val=int(va.size),
                             test_clean=int(te.size), test_all=int(te_all.size)),
                   outlier_dex=a.outlier_dex, seeds=a.seeds, best_val_loss=history, results=results)
    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / 'results.json').write_text(json.dumps(summary, indent=1))
    torch.save(dict(oi_names=names, ctx_names=[str(n) for n in t['ctx_names']], floor=FLOOR,
                    median=med, iqr=iqr, keep=keep, variants=saved), out_dir / 'oi_models.pt')
    np.savez_compressed(out_dir / 'predictions.npz', oi_names=np.array(names), clean=clean,
                        test_idx=te_all, **{k: v for k, v in preds.items()})
    for name in names:
        print(f'\n{name}  (test_clean n={te.size}; dex)')
        print(f"{'method':22s} {'bias':>8s} {'MAD':>8s} {'RMS':>8s} {'p90':>8s} {'wRelRMS':>8s}")
        for method, s in results['test_clean'][name].items():
            print(f"{method:22s} {s['bias_dex']:8.4f} {s['mad_dex']:8.4f} {s['rms_dex']:8.4f} "
                  f"{s['p90_dex']:8.4f} {s['weighted_rel_rms']:8.4f}")


if __name__ == '__main__':
    main()
