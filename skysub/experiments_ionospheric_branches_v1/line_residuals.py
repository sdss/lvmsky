"""O I / N I sky-line residuals after subtracting each method's predicted sky from held-out Science medians.

For every held-out corpus row the predicted Science sky (all 388 coefficients, rendered with that row's own
Science LSF, O2 vector and telluric basis, exactly as the notebook's full-spectrum verification does) is
subtracted from the observed Science median spectrum.  Per line, on the native grid:

    leftover % = 100 * sum_core(resid - b_resid) / sum_core(sky - b_sky)

with core |lambda - lambda_0| <= 3 A and b the median of the 6-12 A side bands.  Negative: the sky line was
over-subtracted.  ``sci_fit`` (the Science row's own decomposition coefficients) is the floor every
prediction is measured against: it carries the decomposition's own miss and any science-field emission.
"""

import numpy as np
import pandas as pd
from astropy.io import fits

from mlp_predictor import cell_parallel, serialization, trainer
from mlp_predictor.data import make_reconstruction_decomposer, make_telluric_row_lookup
from sky_decomp.result_io import load_lsf_surface_cubes, lsf_surface_state_from_cubes

LINES = {'[O I] 5577': 5577.34, '[N I] 5199': 5199.0, '[O I] 6300': 6300.30, '[O I] 6364': 6363.78,
         'O I 7774': 7774.2, 'O I 8446': 8446.4}
OI_COEFS = ('ATOM_Or', 'ATOM_Orc_OI0777', 'ATOM_Orc_OI0845', 'ATOM_Og')
CORE_A, SIDE_A = 3.0, (6.0, 12.0)


def predict_all(path, triplet, rows):
    """Ensemble-mean Science coefficients (all 388) for ``rows``."""
    ens = serialization.load_ensemble(path, device='cpu')
    have = [str(n) for n in triplet['ctx_names']]
    ctx = [have.index(n) for n in ens['ctx_names']]
    arm = {a: np.asarray(triplet[f'coef_{a}'], np.float32)[rows] for a in ('near', 'far')}
    cx = {a: np.asarray(triplet[f'ctx_{a}'], np.float32)[rows][:, ctx] for a in ('near', 'far', 'sci')}
    return np.asarray(trainer.predict_sci_coefficients_default(
        ens, coef_near_phys=arm['near'], coef_far_phys=arm['far'], ctx_near_phys=cx['near'],
        ctx_far_phys=cx['far'], ctx_sci_phys=cx['sci']), np.float64)


def with_dedicated(base, dedicated, triplet):
    """``base`` with the four O I coefficients replaced by the dedicated networks' (as SkyML's oi_mlp run)."""
    from experiments_ionospheric_split_heads_v1.compare import LINES as CMP_LINES
    out, names = base.copy(), [str(n) for n in triplet['coef_names']]
    for j, line in enumerate(CMP_LINES):
        if line in OI_COEFS:
            ok = np.isfinite(dedicated[:, j])
            out[ok, names.index(line)] = dedicated[ok, j]
    return out


def _leftover(wave, obs, sky):
    out = {}
    for key, w0 in LINES.items():
        d = np.abs(wave - w0)
        core, side = d <= CORE_A, (d > SIDE_A[0]) & (d <= SIDE_A[1])
        line = np.nansum(sky[core] - np.nanmedian(sky[side]))
        resid = obs - sky
        out[key] = 100 * np.nansum(resid[core] - np.nanmedian(resid[side])) / line if line > 0 else np.nan
    return out


def line_residuals(triplet, rows, coef_sets, *, corpus_fits, sci_decomp_fits, base_dir, n_moon_knots,
                   n_zodi_knots, factor, n_workers=8):
    """Leftover % per row, line and method; a long DataFrame (row position, method, line, leftover)."""
    corpus_rows = np.asarray(triplet['row_index'], int)[rows]
    with fits.open(corpus_fits, memmap=True) as h:
        obs = np.asarray(h['FLUX_SCI'].data[corpus_rows], np.float64)
        wave = np.asarray(h['WAVE'].data, np.float64)
    with fits.open(sci_decomp_fits, memmap=True) as h:
        o2 = np.asarray(h['VECTOR_O2'].data[corpus_rows], np.float64) if 'VECTOR_O2' in h else None
    lsf_cubes = load_lsf_surface_cubes(str(sci_decomp_fits))
    telluric_for = make_telluric_row_lookup(corpus_fits, verbose=False)

    def work(i):
        model = make_reconstruction_decomposer(wave, n_spline_knots=n_moon_knots, base_dir=base_dir,
                                               n_zodi_spline_knots=n_zodi_knots,
                                               telluric=telluric_for('sci', int(corpus_rows[i])))
        model._set_lsf_state(lsf_surface_state_from_cubes(lsf_cubes, int(corpus_rows[i])))
        mats = model._assemble_refined_matrices()
        if o2 is not None and np.isfinite(o2[i]).any() and np.nansum(np.abs(o2[i])) > 0:
            mats['o2'] = o2[i][None, :]
        res = {}
        for method, coef in coef_sets.items():
            c = model._components_from_coef(np.asarray(coef[i], float), mats)
            sky = sum(np.asarray(c[k], float) for k in ('oh', 'moon', 'zodi', 'diffuse', 'atom', 'orc', 'o2'))
            res[method] = _leftover(wave, obs[i], sky / factor)
        return res

    per_row = cell_parallel.map_indexed(work, range(len(rows)), n_workers=n_workers, verbose=False)
    return pd.DataFrame([dict(row=i, method=m, line=k, leftover=v)
                         for i, r in enumerate(per_row) for m, d in r.items() for k, v in d.items()])


def summary(df, methods):
    """Median leftover % (signed: systematic over/under-subtraction) and median |leftover| (typical size)."""
    g = df.groupby(['method', 'line'])['leftover']
    tab = pd.concat({'median %': g.median().unstack(), 'median |%|': g.apply(lambda x: np.nanmedian(np.abs(x)))
                     .unstack()}, axis=1)
    order = list(LINES)
    return tab.reindex(methods).reindex(columns=pd.MultiIndex.from_product([['median %', 'median |%|'], order]))
