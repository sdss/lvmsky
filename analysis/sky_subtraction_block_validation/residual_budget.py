"""Where does the Sky-fiber residual come from? chi2 budget by pixel category, telluric window and model S/N."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
CHANNELS = {'B': (3600., 5787.), 'R': (5787., 7454.), 'Z': (7454., 9801.)}
FAMILY_LABEL = {'oh': 'OH', 'atom': 'Atomic', 'orc': 'Oxygen recombination', 'o2': 'O₂ airglow'}
CATEGORIES = [f'{FAMILY_LABEL[k]} core (<1 Å)' for k in FAMILY_LABEL] + \
    [f'{FAMILY_LABEL[k]} wing (1–2.5 Å)' for k in FAMILY_LABEL] + \
    ['Outside LSF support (2.5–4 Å)', 'Solar structure (>4 Å from lines)', 'Smooth continuum (>4 Å from lines)']
TELLURIC = {'O₂ B absorption (6860–6930 Å)': (6860., 6930.), 'H₂O (7160–7340 Å)': (7160., 7340.),
            'O₂ A absorption (7590–7700 Å)': (7590., 7700.), 'H₂O (8120–8350 Å)': (8120., 8350.),
            'H₂O (8900–9200 Å)': (8900., 9200.), 'H₂O (9300–9801 Å)': (9300., 9801.)}
STRONG_FRACTION = 0.01
SMOOTH_PIXELS = 41
PARTS = HERE / 'budget_parts'
WORK = Path('/data/sas/sdsswork/users/u6058164/lvm_sky_subtraction_block_validation_drp132_work')


def classify(wave, line_wave, line_flux, line_family, total=None, sigma=None):
    """Exclusive pixel category from distance to, and family of, the nearest strong model line."""
    strong = line_flux >= STRONG_FRACTION * line_flux.max()
    order = np.argsort(line_wave[strong])
    lw = line_wave[strong][order]
    fam = np.asarray(line_family)[strong][order]
    position = np.clip(np.searchsorted(lw, wave), 1, len(lw) - 1)
    left, right = lw[position - 1], lw[position]
    nearest = np.where(abs(wave - left) <= abs(right - wave), position - 1, position)
    distance = abs(wave - lw[nearest])
    category = np.full(len(wave), len(CATEGORIES) - 1)
    names = list(FAMILY_LABEL)
    family_index = np.array([names.index(f) for f in fam])[nearest]
    category[distance >= 2.5] = len(CATEGORIES) - 3
    category[distance >= 4.] = len(CATEGORIES) - 1
    category = np.where(distance < 1., family_index, category)
    category = np.where((distance >= 1.) & (distance < 2.5), 4 + family_index, category)
    if total is not None:
        structure = np.abs(total - lowpass(total)) > sigma
        category = np.where((category == len(CATEGORIES) - 1) & structure, len(CATEGORIES) - 2, category)
    return category, distance


def lowpass(values):
    """NaN-aware running median over SMOOTH_PIXELS native pixels (about 20 Å)."""
    return pd.Series(values).rolling(SMOOTH_PIXELS, center=True, min_periods=SMOOTH_PIXELS // 3).median().to_numpy()


def budget(z, take, categories):
    rows = []
    total = float(np.sum(z[take] ** 2)) if take.any() else np.nan
    for index, name in enumerate(CATEGORIES):
        selected = take & (categories == index)
        if not selected.any():
            continue
        chi = float(np.sum(z[selected] ** 2))
        rows.append(dict(category=name, n_pixels=int(selected.sum()), chi2_sum=chi,
                         chi2_per_pixel=chi / selected.sum(), mean_z=float(z[selected].mean()),
                         chi2_share=chi / total, pixel_share=float(selected.sum() / take.sum())))
    return rows


def analyze(row):
    from skysub.sky_subtraction_block.reconstruction import SkyBlockReconstructor
    from skysub.sky_subtraction.pipeline import render, FACTOR
    e = row['expnum']
    print('Residual budget', e, flush=True)
    category_rows, telluric_rows, scale_rows = [], [], []
    stored = dict(fiberid=[], telescope=[], spectrograph=[], z_low=[], flux_low=[], moon=[], zodi=[], diffuse=[], model=[])
    with SkyBlockReconstructor(row['output'], cframe=row['input']) as r:
        r_wave = np.asarray(r.wave, float)
        fit_mask = r.h['MED_FIT_MASK'].data.astype(bool)
        for fit in r.h['SKY_FIBER_COEF'].data:
            f, i, label, sp = int(fit['fiberid']), int(fit['input_row']), str(fit['telescope']), int(fit['spectrograph'])
            model = r.model(label)
            coef = np.asarray(fit['coef'], float)
            total, comps, matrices = render(model, coef, r.global_state(f))
            total = total / FACTOR
            slices = model._component_slices(matrices)
            line_coef = np.concatenate([coef[slices[k]] for k in model._LINE_KEYS])
            line_flux = model._line_weights() * line_coef[model._line_group]
            family = np.empty(len(line_flux), dtype=object)
            for key, group in model._group_slices.items():
                family[(model._line_group >= group.start) & (model._line_group < group.stop)] = key
            observed = r.observed[i]
            sigma = np.divide(1, np.sqrt(np.maximum(r.ivar[i], 0)), out=np.full(len(r.wave), np.nan), where=r.ivar[i] > 0)
            mask = r.mask[i].copy()
            mask[fit_mask] |= 1
            valid = (mask == 0) & np.isfinite(observed) & np.isfinite(sigma) & (sigma > 0)
            z = np.where(valid, (observed - total) / sigma, np.nan)
            z_low = lowpass(z)
            flux_low = lowpass(np.where(valid, (observed - total) * FACTOR, np.nan))
            for key, value in (('fiberid', f), ('telescope', label), ('spectrograph', sp), ('z_low', z_low.astype('f4')),
                               ('flux_low', flux_low.astype('f4')), ('moon', comps['moon'].astype('f4')),
                               ('zodi', comps['zodi'].astype('f4')), ('diffuse', comps['diffuse'].astype('f4')),
                               ('model', (total * FACTOR).astype('f4'))):
                stored[key].append(value)
            categories, _ = classify(r.wave, model._line_wave, line_flux, family, total, sigma)
            ident = dict(expnum=e, fiberid=f, telescope=label, spectrograph=sp)
            for ch, (lo, hi) in CHANNELS.items():
                inside = (r.wave >= lo) & (r.wave < hi)
                wave = r.wave[inside]
                good = valid & inside & (r.wave >= wave[0] + 5) & (r.wave <= wave[-1] - 5)
                high = (z - z_low)[good]
                high = high[np.isfinite(high)]
                scale_rows.append(dict(**ident, channel=ch, n_pixels=int(good.sum()), rms_z=float(np.sqrt(np.nanmean(z[good] ** 2))),
                                       high_frequency_rms_z=float(1.4826 * np.median(np.abs(high - np.median(high)))),
                                       low_frequency_rms_z=float(np.sqrt(np.nanmean(z_low[good] ** 2))),
                                       low_frequency_chi2_share=float(np.nansum(z_low[good] ** 2) / np.nansum(z[good] ** 2))))
                for item in budget(z, good, categories):
                    category_rows.append(dict(**ident, channel=ch, **item))
            for name, (lo, hi) in TELLURIC.items():
                ch = 'R' if hi <= 7454 else 'Z'
                inside = (r.wave >= lo) & (r.wave < hi) & valid
                if inside.sum() > 10:
                    telluric_rows.append(dict(**ident, channel=ch, window=name, n_pixels=int(inside.sum()),
                                              chi2_per_pixel=float(np.mean(z[inside] ** 2)), mean_z=float(np.mean(z[inside])),
                                              model_to_noise=float(np.median(np.abs(total[inside]) / sigma[inside]))))
    WORK.mkdir(parents=True, exist_ok=True)
    PARTS.mkdir(exist_ok=True)
    np.savez_compressed(WORK / f'residual_arrays_{e}.npz', wave=r_wave, **{k: np.asarray(v) for k, v in stored.items()})
    pd.DataFrame(scale_rows).to_csv(PARTS / f'budget_scale_{e}.csv', index=False)
    pd.DataFrame(category_rows).to_csv(PARTS / f'budget_category_{e}.csv', index=False)
    pd.DataFrame(telluric_rows).to_csv(PARTS / f'budget_telluric_{e}.csv', index=False)
    return e


def assemble(rows):
    sets = {}
    for name in ('category', 'telluric', 'scale'):
        sets[name] = pd.concat([pd.read_csv(PARTS / f'budget_{name}_{row["expnum"]}.csv') for row in rows], ignore_index=True)
        sets[name].to_csv(HERE / f'residual_{name}_metrics.csv', index=False)
    summary = sets['category'].groupby(['expnum', 'channel', 'category']).agg(
        fibers=('fiberid', 'nunique'), n_pixels=('n_pixels', 'sum'), chi2_sum=('chi2_sum', 'sum'),
        median_chi2_per_pixel=('chi2_per_pixel', 'median'), median_mean_z=('mean_z', 'median')).reset_index()
    summary['chi2_per_pixel'] = summary.chi2_sum / summary.n_pixels
    summary['chi2_share'] = summary.chi2_sum / summary.groupby(['expnum', 'channel']).chi2_sum.transform('sum')
    summary['pixel_share'] = summary.n_pixels / summary.groupby(['expnum', 'channel']).n_pixels.transform('sum')
    summary.to_csv(HERE / 'residual_budget_summary.csv', index=False)
    return summary


if __name__ == '__main__':
    rows = json.loads((HERE / 'run_manifest.json').read_text())
    if '--assemble-only' not in sys.argv:
        from concurrent.futures import ProcessPoolExecutor
        import multiprocessing as mp
        with ProcessPoolExecutor(4, mp_context=mp.get_context('fork')) as pool:
            list(pool.map(analyze, rows))
    summary = assemble(rows)
    pd.set_option('display.width', 250)
    print(summary[summary.expnum.isin([3469, 16819, 60919])].round(3).to_string(index=False))
