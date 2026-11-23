"""Build the [O I] experiment table: filtered corpus triplet, split, and deployed-ensemble predictions.

Run from ``skysub/`` (``mlp_predictor`` and ``sky_decomp`` are top-level packages there):

    python -m experiments_oi_dedicated_mlp_v1.prepare --corpus-dir DIR --ensemble PATH --out FILE

The triplet, filters and context augments follow cell 3 of
``notebook_sky_interpolation_triplet_dual_encoder_group_mlp_split_zodi_module.ipynb``,
including the moon/zodi reversal gate on the full ``*_decomp_*`` files.  The split is the notebook's
night-level, moon-phase-stratified split with seed 42 on THIS filtered set, so it is
not guaranteed to equal the deployed ensemble's split.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mlp_predictor import data, serialization, trainer  # noqa: E402
from mlp_predictor.config import DataConfig  # noqa: E402
from mlp_predictor.ml_utils import moon_phase_deg_from_ctx, split_indices_by_moon_phase  # noqa: E402

OI_NAMES = ('ATOM_Or', 'ATOM_Orc_OI0777', 'ATOM_Orc_OI0845', 'ATOM_Og')
STEM = 'lvmsframe_median_stack_1.3.2_knox_gaia1over100'


def native_wave(path):
    with fits.open(path) as h:
        wave = np.asarray(h['WAVE'].data, np.float64)
    return wave if wave.ndim == 1 else wave[0]


def build_triplet(cfg):
    triplet = data.build_triplet_coef_dataset(
        input_fits_path=cfg.input_fits_meta,
        sky_near_decomp_fits_path=cfg.coef_fits('sky1'),
        sky_far_decomp_fits_path=cfg.coef_fits('sky2'),
        sci_decomp_fits_path=cfg.coef_fits('sci'),
        context_columns=list(cfg.context_columns),
        return_chi2=True,
    )
    filtered = data.apply_triplet_filters(
        triplet,
        thin_every_n=1, chi2_qmax=90.0, chi2_min=0.0, chi2_max=10.0,
        hard_coef_bounds={'feo': (0.0, 36.01), 'atom_k': (0.0, 10.01)},
        kappa=8.0, kappa_iter=3, oh_kappa=6.0, oh_kappa_iter=3,
        exclude_field_regions=[data.LMC_EXCLUSION, data.SMC_EXCLUSION],
        reversal_decomp_fits={
            'near': cfg.decomp_fits('sky1'), 'far': cfg.decomp_fits('sky2'), 'sci': cfg.decomp_fits('sci')},
        reversal_wave=native_wave(cfg.input_fits_for_basis),
        reversal_min_component_frac=0.05, reversal_min_separation=0.0,
        colour_excess_input_fits=cfg.input_fits_flux,
        colour_excess_max=data.SCI_COLOUR_EXCESS_MAX,
        diffuse_zeroed_frac=data.DIFFUSE_ZEROED_FRAC,
        show_plots=False,
    )
    data._augment_triplet_with_ecliptic(filtered, force=True, meta_fits_path=cfg.input_fits_meta)
    data._augment_triplet_with_physics_priors(filtered, force=True)
    data._augment_triplet_with_moon_model(filtered, cfg.decomp_prefix, force=True)
    return filtered


def reorder_ctx(filtered, names):
    have = [str(n) for n in filtered['ctx_names']]
    missing = [n for n in names if n not in have]
    if missing:
        raise ValueError(f'context features missing for the ensemble: {missing}')
    idx = [have.index(n) for n in names]
    return {arm: np.asarray(filtered[f'ctx_{arm}'], np.float64)[:, idx] for arm in ('near', 'far', 'sci')}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--corpus-dir', required=True)
    p.add_argument('--ensemble', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()

    cfg = DataConfig(decomp_data_root=str(a.corpus_dir), decomp_stem=STEM)
    filtered = build_triplet(cfg)
    names = [str(n) for n in filtered['coef_names']]
    ensemble = serialization.load_ensemble(a.ensemble)
    if names != list(ensemble['coef_names']):
        raise ValueError('corpus coefficient names/order differ from the ensemble')
    ctx = reorder_ctx(filtered, list(ensemble['ctx_names']))

    coef = {arm: np.asarray(filtered[f'coef_{arm}'], np.float64) for arm in ('near', 'far', 'sci')}
    pred = trainer.predict_sci_coefficients_default(
        ensemble, coef['near'], coef['far'], ctx['near'], ctx['far'], ctx['sci'])

    oi = [names.index(n) for n in OI_NAMES]
    tr, va, te = split_indices_by_moon_phase(
        filtered['obstime_mjd'], moon_phase_deg_from_ctx(filtered), seed=42)
    row = np.asarray(filtered['row_index'], np.int64)
    with fits.open(cfg.input_fits_meta) as h:
        meta = h['META'].data
        tile = np.asarray(meta['tile_id'])[row]
        expnum = np.asarray(meta['expnum'])[row]

    out = dict(
        oi_names=np.array(OI_NAMES), ctx_names=np.array(ensemble['ctx_names']),
        row_index=row, expnum=expnum, tile_id=tile,
        obstime_mjd=np.asarray(filtered['obstime_mjd'], np.float64),
        train_idx=tr, val_idx=va, test_idx=te,
        deployed_pred=np.asarray(pred, np.float64)[:, oi],
    )
    for arm in ('near', 'far', 'sci'):
        out[f'coef_{arm}'] = coef[arm][:, oi]
        out[f'err_{arm}'] = np.asarray(filtered[f'coef_err_{arm}'], np.float64)[:, oi]
        out[f'ctx_{arm}'] = ctx[arm]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f'rows {row.size}: train {tr.size}, val {va.size}, test {te.size} -> {a.out}')


if __name__ == '__main__':
    main()
