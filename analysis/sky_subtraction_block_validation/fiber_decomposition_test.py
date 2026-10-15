"""Nine reference Science fibers per exposure: in-fiber Sky Decomposition versus the global block LSF (R and Z).

Variants, all on the same data, fit pixels and weights:
  fiber        stored SolarBlue fit (its own R/Z LSF and coefficients), re-rendered
  block_same   block R/Z LSF with the stored coefficients
  fiber_refit  coefficients re-estimated by the production QP with the fiber's own LSF fixed
  block_refit  coefficients re-estimated by the production QP with the block R/Z LSF fixed
B keeps the fiber's own fitted LSF in every variant.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ[key] = '1'
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from pathlib import Path
import glob
from functools import lru_cache
import json
import multiprocessing as mp
import sys
import time
import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
WORK = Path('/data/sas/sdsswork/users/u6058164/lvm_sky_subtraction_block_validation_drp132_work')
SOLARBLUE = Path('/data/sas/sdsswork/users/u6058164/lvm_solarblue_all_fibers_drp132_moonalt60_fli095')
EXTRA = Path('/data/sas/sdsswork/users/u6058164/lvm_lsf_epoch_coverage_drp132/additional_clean_fits')
EXPOSURES = (4108, 16819, 40620, 44808, 60919)
FIT = ('R', 'Z')
CHANNELS = {'R': (5787., 7454.), 'Z': (7454., 9801.)}
CHANNEL_INDEX = {'B': 0, 'R': 1, 'Z': 2}
VARIANTS = ('fiber', 'block_same', 'fiber_refit', 'block_refit')
GROUPS = {'OH core (<1 Å)': (0,), 'OH wing (1–2.5 Å)': (4,), 'Other line cores (<1 Å)': (1, 2, 3), 'Other line wings (1–2.5 Å)': (5, 6, 7),
          'Outside LSF support (2.5–4 Å)': (8,), 'Continuum (>4 Å from lines)': (9, 10)}


def reference_file(expnum):
    for root in (SOLARBLUE, EXTRA):
        for path in sorted(glob.glob(str(root / 'reference_models/*.fits'))):
            with fits.open(path, memmap=True) as h:
                rows = np.flatnonzero(h['INPUT_META'].data['expnum'] == expnum)
            if len(rows):
                return path, rows
    raise ValueError(f'No reference model for exposure {expnum}')


@lru_cache(maxsize=None)
def gaia_status(expnum):
    """{fiberid: (clean_1pct, gaia_ratio)} from the Gaia-verified fiber table; the reference-model INPUT_META carries no Gaia values."""
    table = pd.read_parquet('/data/sas/sdsswork/users/u6058164/lvm_lsf_epoch_coverage_drp132/fiber_exposure_channel.parquet',
                            columns=['expnum', 'fiberid', 'telescope', 'channel', 'gaia_ratio', 'gaia_clean_1pct'])
    table = table[(table.expnum == expnum) & (table.telescope == 'Sci') & (table.channel == 'Z')]
    return {int(f): (bool(c), float(r)) for f, c, r in zip(table.fiberid, table.gaia_clean_1pct, table.gaia_ratio)}


def selected_record(expnum):
    for root in (SOLARBLUE, EXTRA):
        table = pd.read_csv(root / 'selected_exposures.csv')
        if (table.expnum == expnum).any():
            return table[table.expnum == expnum].iloc[0].to_dict()
    raise ValueError(f'Exposure {expnum} is not in a selected_exposures.csv')


def fiber_state(base, h, row):
    state = deepcopy(base)
    config = {str(c['channel']).strip(): c for c in h['LSF_CONFIG'].data}
    for ch in 'BRZ':
        state.coefficients[ch] = np.asarray(h['LSF_COEF'].data[row, CHANNEL_INDEX[ch]], float)
        state.knot_vectors[ch] = np.asarray(h['LSF_KNOTS'].data[row, CHANNEL_INDEX[ch], :int(config[ch]['n_knots'])], float)
        state.degrees[ch] = int(config[ch]['degree'])
    assert abs(float(config['R']['offset_half_width_angstrom']) - state.config['offset_half_width_angstrom']) < 1e-12
    return state


def block_state(fiber, block):
    state = deepcopy(fiber)
    for ch in FIT:
        state.coefficients[ch] = block.coefficients[ch].copy()
        state.knot_vectors[ch] = block.knot_vectors[ch].copy()
        state.degrees[ch] = block.degrees[ch]
        state.metrics[ch] = dict(status='Solved', model='sky_block_centered', quality='measured')
    return state


def render_model(dec, coef, state, t_o2):
    from skysub.sky_subtraction.rendering import render
    dec.t_o2 = float(t_o2)
    total, components, matrices = render(dec, coef, state)
    return total, components, matrices


def categories(dec, coef, matrices, wave, total, sigma):
    from residual_budget import classify
    slices = dec._component_slices(matrices)
    line_coef = np.concatenate([coef[slices[k]] for k in dec._LINE_KEYS])
    line_flux = dec._line_weights() * line_coef[dec._line_group]
    family = np.empty(len(line_flux), dtype=object)
    for key, group in dec._group_slices.items():
        family[(dec._line_group >= group.start) & (dec._line_group < group.stop)] = key
    return classify(wave, dec._line_wave, line_flux, family, total, sigma)[0]


def analyze_fiber(task):
    """One fiber of one exposure; returns its metric rows and spectra."""
    expnum, k = task
    from skysub import decompose_parallel as dp
    from skysub.sky_subtraction_block.reconstruction import SkyBlockReconstructor
    started = time.perf_counter()
    manifest = {r['expnum']: r for r in json.loads((HERE / 'run_manifest.json').read_text())}
    dp.decompose_in_process(str(WORK / f'stack_{expnum}.fits'), rows=(), kinds=('sci',))
    path, rows = reference_file(expnum)
    wave = np.asarray(dp._WORKER_WAVE, float)
    metric_rows, category_rows = [], []
    with SkyBlockReconstructor(manifest[expnum]['output'], cframe=manifest[expnum]['input']) as r, fits.open(path, memmap=True) as h:
        row = rows[k]
        meta = h['INPUT_META'].data
        names = h['COEF'].columns.names
        f = int(meta['fiberid'][row])
        sp = int(meta['spectrographid'][row])
        stored_coef = np.array([h['COEF'].data[row][n] for n in names], float)
        t_o2 = float(h['META'].data['t_o2'][row])
        stored_best = np.asarray(h['BESTFIT_LSF'].data[row], float)
        fstate = fiber_state(r.global_state(f), h, row)
        bstate = block_state(fstate, r.global_state(f))
        flux = np.asarray(dp._WORKER_FLUX['sci'][k], float) * dp._WORKER_FACTOR
        ivar = dp._fit_ivar_row('sci', k, flux)
        ivar[dp._science_line_mask_for_row(k)] = 0.0
        sigma = np.divide(1, np.sqrt(ivar), out=np.full(len(wave), np.nan), where=ivar > 0)
        dec = dp._telluric_decomposer('sci', k)
        models, comps, coefs, t_values = {}, {}, {}, {}
        models['fiber'], comps['fiber'], matrices_fiber = render_model(dec, stored_coef, fstate, t_o2)
        coefs['fiber'], t_values['fiber'] = stored_coef, t_o2
        check = dict(expnum=expnum, fiberid=f, max_abs_render_minus_stored_bestfit=float(np.nanmax(abs(models['fiber'] - stored_best))),
                     median_abs_stored_bestfit=float(np.nanmedian(abs(stored_best))))
        models['block_same'], _, _ = render_model(dec, stored_coef, bstate, t_o2)
        coefs['block_same'], t_values['block_same'] = stored_coef, t_o2
        for name, state in (('fiber_refit', fstate), ('block_refit', bstate)):
            result, _ = dp._fit_worker_row('sci', k, flux, ivar.copy(), fixed_lsf_state=state)
            if str(result.fit_status) != 'Solved':
                raise RuntimeError(f'{expnum} fiber {f} {name}: {result.fit_status}')
            coefs[name], t_values[name] = np.asarray(result.coef, float), float(result.t_o2)
            models[name], comps_name, _ = render_model(dec, coefs[name], state, t_values[name])
            np.testing.assert_allclose(models[name], np.asarray(result.bestfit_lsf, float), rtol=2e-5, atol=2e-5)
            comps[name] = comps_name
            print(expnum, f, name, 'fit', round(time.perf_counter() - started, 1), 's', flush=True)
        dec.t_o2 = t_o2
        cat = categories(dec, stored_coef, matrices_fiber, wave, models['fiber'], sigma)
        for ch, (lo, hi) in CHANNELS.items():
            inside = (wave >= lo) & (wave < hi)
            edge = wave[inside]
            good = inside & (ivar > 0) & (wave >= edge[0] + 5) & (wave <= edge[-1] - 5)
            for variant in VARIANTS:
                z = (flux - models[variant]) / sigma
                dz = (models[variant] - models['fiber']) / sigma
                metric_rows.append(dict(expnum=expnum, fiberid=f, spectrograph=sp, channel=ch, variant=variant, n_pixels=int(good.sum()),
                                        chi2_per_pixel=float(np.mean(z[good] ** 2)), within_1sigma=float(np.mean(abs(z[good]) <= 1)),
                                        mean_z=float(np.mean(z[good])), rms_flux_residual=float(np.sqrt(np.mean((flux - models[variant])[good] ** 2))),
                                        rms_model_difference_sigma=float(np.sqrt(np.mean(dz[good] ** 2))), t_o2=t_values[variant]))
                for group, indices in GROUPS.items():
                    take = good & np.isin(cat, indices)
                    if take.any():
                        category_rows.append(dict(expnum=expnum, fiberid=f, spectrograph=sp, channel=ch, variant=variant, category=group,
                                                  n_pixels=int(take.sum()), chi2_sum=float(np.sum(z[take] ** 2)),
                                                  chi2_per_pixel=float(np.mean(z[take] ** 2)),
                                                  model_difference_sigma=float(np.sqrt(np.mean(dz[take] ** 2)))))
        arrays = dict(flux=flux.astype('f4'), sigma=sigma.astype('f4'), fit_pixels=ivar > 0, fiberid=f, spectrograph=sp, wave=wave)
        for v in VARIANTS:
            arrays[f'model_{v}'] = models[v].astype('f4')
            arrays[f'coef_{v}'] = coefs[v]
        for v in ('fiber', 'block_refit'):
            for key in ('oh', 'moon', 'zodi', 'diffuse', 'atom', 'orc', 'o2'):
                arrays[f'comp_{v}_{key}'] = np.asarray(comps[v][key], 'f4')
        arrays['coef_names'] = np.array(names)
    print(expnum, 'fiber', f, 'done', round(time.perf_counter() - started, 1), 's', flush=True)
    return expnum, k, metric_rows, category_rows, check, arrays


def prepare(expnum):
    from analysis.solarblue_all_fibers import run_all_fibers as rf
    WORK.mkdir(parents=True, exist_ok=True)
    stack = WORK / f'stack_{expnum}.fits'
    if not stack.is_file():
        rf.build_exposure_input(selected_record(expnum), stack)
    return expnum


if __name__ == '__main__':
    exposures = [int(x) for x in sys.argv[1:]] or list(EXPOSURES)
    for e in exposures:
        prepare(e)
    results = {e: {} for e in exposures}
    tasks = [(e, k) for e in exposures for k in range(9)]
    started = time.perf_counter()
    with ProcessPoolExecutor(min(len(tasks), 64), mp_context=mp.get_context('fork')) as pool:
        for e, k, *rest in pool.map(analyze_fiber, tasks):
            results[e][k] = rest
    for e in exposures:
        parts = [results[e][k] for k in range(9)]
        pd.DataFrame([x for p in parts for x in p[0]]).to_csv(HERE / f'fiber_decomp_metrics_{e}.csv', index=False)
        pd.DataFrame([x for p in parts for x in p[1]]).to_csv(HERE / f'fiber_decomp_categories_{e}.csv', index=False)
        pd.DataFrame([p[2] for p in parts]).to_csv(HERE / f'fiber_decomp_render_check_{e}.csv', index=False)
        merged = {key: np.stack([p[3][key] for p in parts]) for key in parts[0][3] if key not in ('wave', 'coef_names')}
        np.savez_compressed(WORK / f'fiber_decomposition_{e}.npz', wave=parts[0][3]['wave'], coef_names=parts[0][3]['coef_names'], **merged)
    print('Finished', round(time.perf_counter() - started, 1), 's', flush=True)
