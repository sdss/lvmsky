"""R/Z comparison of DRP, per-fiber in-situ LSF and block LSF through the same ML sky; B keeps the survey LSF."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ[key] = '1'
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from pathlib import Path
import html
import json
import multiprocessing as mp
import sys
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
CORPUS = Path('/data/sas/sdsswork/users/u6058164/lvm_lsf_epoch_coverage_drp132')
WORK = Path('/data/sas/sdsswork/users/u6058164/lvm_sky_subtraction_block_validation_drp132_work')
PARAMETERS = [Path('/data/sas/sdsswork/users/u6058164/lvm_solarblue_all_fibers_drp132_moonalt60_fli095/parameters'),
              CORPUS / 'additional_clean_fits/parameters']
FACTOR = 1e14
EXPOSURES = (4108, 16819, 40620, 44808, 60919)
CHANNELS = {'B': (3600., 5787.), 'R': (5787., 7454.), 'Z': (7454., 9801.)}
CHANNEL_INDEX = {'B': 0, 'R': 1, 'Z': 2}
FIT = ('R', 'Z')
_CONTEXT = {}


def parameter_files():
    return [p for root in PARAMETERS for p in sorted(root.glob('solarblue-shard-*.fits'))]


def individual_profiles(expnum):
    """{fiberid: {channel: (coefficients, knots, degree)}} from the per-exposure SolarBlue Science fits."""
    profiles = {}
    for path in parameter_files():
        with fits.open(path, memmap=True) as h:
            meta = h['INPUT_META'].data
            rows = np.flatnonzero(meta['expnum'] == expnum)
            if not len(rows):
                continue
            config = {str(c['channel']).strip(): c for c in h['LSF_CONFIG'].data}
            coefficient, knots = h['LSF_COEF'].data, h['LSF_KNOTS'].data
            for i in rows:
                profiles[int(meta['fiberid'][i])] = {
                    ch: (np.asarray(coefficient[i, CHANNEL_INDEX[ch]], float), np.asarray(knots[i, CHANNEL_INDEX[ch], :int(config[ch]['n_knots'])], float),
                         int(config[ch]['degree'])) for ch in FIT}
            return profiles
    raise ValueError(f'No SolarBlue parameter FITS contains exposure {expnum}')


def accepted_channels(expnum):
    table = pd.read_parquet(CORPUS / 'fiber_exposure_channel.parquet', columns=['expnum', 'fiberid', 'channel', 'telescope', 'accepted'])
    table = table[(table.expnum == expnum) & (table.telescope == 'Sci') & table.accepted & table.channel.isin(FIT)]
    return {(int(f), c) for f, c in zip(table.fiberid, table.channel)}


def individual_state(base, profile, accepted, fiberid):
    state = deepcopy(base)
    for ch in FIT:
        if (fiberid, ch) in accepted:
            coefficient, knots, degree = profile[ch]
            state.coefficients[ch] = coefficient.copy()
            state.knot_vectors[ch] = knots.copy()
            state.degrees[ch] = degree
            state.metrics[ch] = dict(status='Solved', model='individual_insitu', quality='measured')
    return state


def _render(index):
    from skysub.sky_subtraction.rendering import sky_parts
    c = _CONTEXT
    r = c['reconstructor']
    row = c['science'][index]
    fiberid = int(c['fiberid'][index])
    base = r.state('', index)
    state = individual_state(base, c['profiles'][fiberid], c['accepted'], fiberid)
    sky = sky_parts(r.model_sci, r.coef, r.oh_coef, state, r.arms, r.settings).total
    return int(row), (r.observed[row] - sky).astype('f4'), [(fiberid, ch) in c['accepted'] for ch in FIT]


def compute_flux(row, jobs=8):
    """FLUX = observed - ML sky rendered with each fiber's own accepted R/Z profile; everything else unchanged."""
    from skysub.sky_subtraction_block.reconstruction import SkyBlockReconstructor
    e = row['expnum']
    with SkyBlockReconstructor(row['output'], cframe=row['input']) as r:
        science = np.asarray(r.h['LSF_FIBERS'].data['input_row'], int)
        fiberid = np.asarray(r.h['LSF_FIBERS'].data['fiberid'], int)
        profiles = individual_profiles(e)
        _CONTEXT.update(reconstructor=r, science=science, fiberid=fiberid, profiles=profiles, accepted=accepted_channels(e))
        shape = r.h['FLUX'].data.shape
        flux = np.full(shape, np.nan, 'f4')
        accepted = np.zeros((shape[0], len(FIT)), bool)
        with ProcessPoolExecutor(jobs, mp_context=mp.get_context('fork')) as pool:
            for row_index, values, ok in pool.map(_render, range(len(science)), chunksize=8):
                flux[row_index] = values
                accepted[row_index] = ok
        _CONTEXT.clear()
    np.savez_compressed(WORK / f'individual_flux_{e}.npz', flux=flux, accepted=accepted)
    return e


def flux_figure(row, r_wave, slit, science, drp, block, individual, accepted):
    e = row['expnum']
    columns = (('DRP SFrame', drp), ('Per-fiber in-situ LSF', individual), ('Block LSF', block))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#888888')
    fig, axes = plt.subplots(3, 3, figsize=(18, 11), layout='constrained')
    rows = []
    for j, (ch, (lo, hi)) in enumerate(CHANNELS.items()):
        take = (r_wave >= lo) & (r_wave < hi)
        valid = np.isfinite(drp[science][:, take]) & np.isfinite(block[science][:, take]) & np.isfinite(individual[science][:, take])
        if ch in FIT:
            valid &= accepted[science][:, FIT.index(ch)][:, None]
        pooled = np.concatenate([values[science][:, take][valid] for _, values in columns])
        limit = float(np.nanpercentile(np.abs(pooled), 98))
        for k, (name, values) in enumerate(columns):
            image = np.full((len(slit), take.sum()), np.nan)
            image[science] = np.where(valid, values[science][:, take], np.nan)
            ax = axes[j, k]
            im = ax.imshow(image, aspect='auto', origin='lower', extent=[r_wave[take][0], r_wave[take][-1], slit['fiberid'].min() - .5, slit['fiberid'].max() + .5],
                           cmap=cmap, vmin=-limit, vmax=limit, interpolation='nearest')
            ax.set_title(f'{ch} · {name}')
            ax.set_ylabel('Physical fiber ID')
            ax.set_xlabel('Air wavelength (Å)')
            for sp in (1, 2):
                ax.axhline(np.max(slit['fiberid'][slit['spectrographid'] == sp]) + .5, color='black', lw=.6)
            if k == 2:
                fig.colorbar(im, ax=axes[j, :], label='Flux (10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹)', shrink=.85)
            rows.append(dict(expnum=e, channel=ch, variant=name, shared_abs_limit=limit, valid_pixels=int(valid.sum()), n_fibers=int(valid.any(axis=1).sum()),
                             median_abs_flux=float(np.median(np.abs(values[science][:, take][valid]))),
                             p90_abs_flux=float(np.percentile(np.abs(values[science][:, take][valid]), 90)),
                             median_flux=float(np.median(values[science][:, take][valid]))))
    fig.suptitle(f'Exposure {e}: remaining source + sky residual (same ML sky, B survey LSF; R/Z from the stated source)', fontsize=16)
    fig.savefig(HERE / f'individual_flux_maps_{e}.png', dpi=100)
    plt.close(fig)
    return rows


def figures(row):
    e = row['expnum']
    print('Figures', e, flush=True)
    from skysub.sky_subtraction_block.reconstruction import SkyBlockReconstructor
    flux = np.load(WORK / f'individual_flux_{e}.npz')
    with SkyBlockReconstructor(row['output'], cframe=row['input']) as r, fits.open(row['drp'], memmap=True) as drp:
        slit = r.h['SLITMAP'].data
        science = np.asarray(r.h['LSF_FIBERS'].data['input_row'], int)
        np.testing.assert_array_equal(drp['WAVE'].data, r.wave)
        scale = u.Unit(drp['FLUX'].header['BUNIT']).to(u.erg / u.s / u.cm**2 / u.AA)
        reference = np.full(r.h['FLUX'].data.shape, np.nan)
        usable = (drp['MASK'].data == 0) & (drp['IVAR'].data > 0) & (r.h['MASK'].data == 0) & (r.h['IVAR'].data > 0)
        reference[science] = np.where(usable[science], np.asarray(drp['FLUX'].data[science], float) * scale * FACTOR, np.nan)
        block = np.where(usable, np.asarray(r.h['FLUX'].data, float) * FACTOR, np.nan)
        individual = np.where(usable, flux['flux'].astype(float) * FACTOR, np.nan)
        flux_rows = flux_figure(row, r.wave, slit, science, reference, block, individual, flux['accepted'])
    return flux_rows


def section(flux):
    f = flux.copy()
    f['variant'] = f.variant.map({'DRP SFrame': 'DRP', 'Per-fiber in-situ LSF': 'Per-fiber', 'Block LSF': 'Block'})
    table = '<table><tr><th>Exposure</th><th>Channel</th><th>Fibers</th><th>Pixels</th><th>Median |FLUX|: DRP</th><th>Per-fiber in-situ</th><th>Block</th><th>90th percentile |FLUX|: DRP</th><th>Per-fiber in-situ</th><th>Block</th></tr>'
    for (e, ch), g in f.groupby(['expnum', 'channel']):
        g = g.set_index('variant')
        table += f'<tr><td>{e}</td><td>{ch}</td><td>{int(g.n_fibers.iloc[0])}</td><td>{int(g.valid_pixels.iloc[0])}</td>' + ''.join(f'<td>{g.loc[v, "median_abs_flux"]:.3f}</td>' for v in ('DRP', 'Per-fiber', 'Block')) \
            + ''.join(f'<td>{g.loc[v, "p90_abs_flux"]:.3f}</td>' for v in ('DRP', 'Per-fiber', 'Block')) + '</tr>'
    table += '</table>'
    body = '''<div id="individual-comparison"><h2>R and Z: remaining source + sky residual with each LSF</h2>
<p>The block LSF is fitted and compared in R and Z only. B is not fitted: every Science fiber keeps its survey LSF in all variants below, so B differs between columns only through the other changes in the subtraction. Each variant uses the same ML component coefficients, Science-median OH correction and B continuum transfer; only the R/Z LSF changes. <b>DRP</b>: standard DRP SFrame. <b>Per-fiber in-situ</b>: each fiber's own R/Z surface from the SolarBlue fit of that exposure's Science spectrum (accepted channels only). <b>Block</b>: the sky-only surface evaluated at that fiber. Fiber, channel and pixel selection are identical in the three columns (common valid pixels, accepted per-fiber R/Z fits).</p>
<h3>Remaining source + sky residual</h3>
<p>FLUX = observed − subtracted sky, in 10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹. Source emission is included, so the values are not a zero-target residual; they compare variants on identical pixels. Smaller values mean less residual only when the source contribution is the same, which holds here because the pixels and fibers are identical.</p>''' + table
    for e in sorted(flux.expnum.unique()):
        body += f'<h3>Exposure {e}</h3><img src="individual_flux_maps_{e}.png" alt="Flux maps of DRP, per-fiber in-situ and block subtraction for exposure {e}"><p>Rows are B, R and Z; each row has one symmetric scale shared by the three columns (pooled 98th absolute percentile; larger values saturate). Gray: masked, non-Science, or without an accepted per-fiber R/Z fit (the Gaia-cleanliness criterion leaves only 293 fibers in 16819 and 15 in 44808, so those two maps are mostly gray in R and Z). Horizontal lines separate spectrographs.</p>'
    body += '<p><a href="individual_flux_summary.csv">FLUX summary</a></p></div>'
    return body


if __name__ == '__main__':
    manifest = {r['expnum']: r for r in json.loads((HERE / 'run_manifest.json').read_text())}
    rows = [manifest[e] for e in EXPOSURES]
    WORK.mkdir(exist_ok=True, parents=True)
    if '--figures-only' not in sys.argv:
        for row in rows:
            print('Rendering', row['expnum'], flush=True)
            compute_flux(row)
    flux_rows = []
    for row in rows:
        flux_rows += figures(row)
    flux = pd.DataFrame(flux_rows)
    flux.to_csv(HERE / 'individual_flux_summary.csv', index=False)
    (HERE / 'individual_comparison_section.html').write_text(section(flux))
    print(flux.pivot_table(index=['expnum', 'channel'], columns='variant', values='median_abs_flux').round(3).to_string())
