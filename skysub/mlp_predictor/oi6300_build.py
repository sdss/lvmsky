"""Build, validate and save the Science [O I] 6300 model (``oi6300.OI6300Model``) from a decomposed median stack.

One command from a decomposed median stack to the saved model; run from ``skysub/`` (the model class and the
method are described in ``oi6300.py``)::

    python -m mlp_predictor.oi6300_build --decomp-root ROOT --stem STEM --out DIR
        [--sci-lines [oh o2 5577 7774 8446]] [--no-sky-field] [--members 5] [--name NAME] [--no-cv]
        [--rebuild-context] [--workers 4] [--corpus FILE]

Without options it builds the default model: Science lines OH, O2, 5577, 7774, 8446 and the night field.  It first
builds what is missing of the ionospheric context of the stack (``iono_corpus.build_context``: ionosphere, ROTI and
thermosphere tables, then the corpus), so a new DRP version needs this one command after the decomposition.

Input files: ``ROOT/STEM.fits`` (median stack: WAVE, FLUX_*, META), ``ROOT/STEM_meta_only.fits``, the decompositions
``ROOT/STEM_decomp_{sci,sky1}<suffix>.fits`` (RESID) and ``ROOT/STEM_{sci,sky1,sky2}_meta_coef<suffix>.fits`` (COEF).

Steps, one function each; call them from a notebook to change any of them:

    iono_corpus.build_context -> context tables and corpus of the stack (skipped where they exist)
    build_table     -> table of the ``oi6300.inputs`` columns and per-row quantities of the selection
    row_sets        -> clean rows, LMC+SMC, application domain, gas-poor, training rows
    control_lines   -> Science/SkyNear OH and O2 ratios (inputs ``'oh'``, ``'o2'``)
    night_field     -> night field of the sky telescopes (inputs of ``sky_field=True``)
    cross_validate  -> test-night predictions: the expected accuracy of the saved model
    main            -> all of the above, the fit on all training rows, ``<name>.pkl``, ``<name>.json`` (provenance and
                       scores) and ``<name>_cv.csv`` (the test-night prediction of every training row)

Spectra are read on the native wavelength grid; only pixel windows are cut out of it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

from . import data, iono_corpus, oi6300, sky_field
from .config import PipelineConfig

SKYSUB = Path(__file__).resolve().parents[1]
ARMS = {'sci': ('sci', 'sci', 'FLUX_SCI'), 'near': ('sky_near', 'sky1', 'FLUX_SKY_NEAR'),
        'far': ('sky_far', 'sky2', 'FLUX_SKY_FAR')}  # arm -> META prefix, decomposition tag, flux HDU
POINTINGS = ('sci', 'sky_near', 'sky_far', 'skye', 'skyw')  # a pointing at ra = dec = 0 has broken coordinates
HALPHA, HALPHA_WINDOW = 6562.8, (6560.3, 6565.3)  # [A]
GAS_POOR_MAX_HALPHA, GAS_POOR_MIN_ABS_B = 2.0, 10.0  # Halpha(Science) - Halpha(SkyNear) [R]; |b| [deg]
N_OH, O2_COLUMN = 40, 'O2_b01'  # control lines: the N_OH OH groups brightest in SkyNear, and O2 b(0-1)
DEFAULT_SCI_LINES = ('oh', 'o2', '5577', '7774', '8446')  # the default model, with the night field
N_FOLDS, SEP_BINS = 5, ((0, 2), (2, 4), (4, 6), (6, 10), (10, 90))  # runs by nights; Science-SkyNear separation [deg]


def sigma(r):
    """Robust scatter, 1.4826 MAD."""
    return float(1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r))))


def _context(dcfg, obstime_mjd):
    """SkyML observing context of each telescope on the corpus rows: arm -> (n_rows, n_ctx), and the ctx names."""
    t = data.build_triplet_coef_dataset(
        input_fits_path=dcfg.input_fits_meta, sky_near_decomp_fits_path=dcfg.coef_fits('sky1'),
        sky_far_decomp_fits_path=dcfg.coef_fits('sky2'), sci_decomp_fits_path=dcfg.coef_fits('sci'),
        context_columns=list(dcfg.context_columns))
    data._augment_triplet_with_ecliptic(t, meta_fits_path=dcfg.input_fits_meta)
    data._augment_triplet_with_physics_priors(t)
    order = np.argsort(t['obstime_mjd'])
    k = iono_corpus._match(np.asarray(t['obstime_mjd'], float)[order], obstime_mjd)
    hit, row = k >= 0, order[np.maximum(k, 0)]
    print(f'observing context: {(~hit).sum()} of {hit.size} corpus rows have no match (context NaN)', flush=True)
    ctx = {a: np.where(hit[:, None], np.asarray(t[f'ctx_{a}'], float)[row], np.nan) for a in ARMS}
    return ctx, list(t['ctx_names'])


def _window(hdu, rows, wave, lo, hi):
    """Native-grid pixels of ``lo < wave < hi`` for file rows ``rows`` of an image HDU (read in file order)."""
    sl = slice(*np.searchsorted(wave, [lo, hi]))
    o = np.argsort(rows)
    out = np.empty((rows.size, sl.stop - sl.start))
    out[o] = np.asarray(hdu.data[rows[o], sl], float)
    return wave[sl], out


def build_table(dcfg, corpus):
    """Model table and selection quantities for every corpus row.

    Returns ``tab`` (the ionosphere corpus plus ``<arm>_ctx_<name>`` and ``sep_<a>_<b>`` [deg]; see
    ``oi6300.inputs``) and ``aux``: ``ctx_names``, ``file_row`` (row of the stack and decomposition files),
    ``valid_coords``, ``lmc``, ``smc``, ``gal_b`` [deg], ``n_fibres`` and ``log_ew`` per arm, and
    ``halpha_excess`` = Halpha(Science) - Halpha(SkyNear) [R] in the decomposition residuals.
    """
    ic, _ = iono_corpus.load(corpus)
    tab = dict(ic)
    meta = Table.read(dcfg.input_fits_meta, hdu='META')
    expnum = np.asarray(meta['expnum'])
    if np.unique(expnum).size != expnum.size or np.unique(ic['expnum']).size != ic['expnum'].size:
        raise ValueError('an exposure number is repeated in the median stack or in the corpus')
    row_of = {e: i for i, e in enumerate(expnum)}
    file_row = np.array([row_of[e] for e in ic['expnum']])
    col = lambda c: np.asarray(meta[c], float)[file_row]  # noqa: E731
    pos = {a: SkyCoord(col(f'{p}_ra'), col(f'{p}_dec'), unit='deg') for a, (p, _, _) in ARMS.items()}
    for a, b in (('sci', 'near'), ('sci', 'far'), ('near', 'far')):
        tab[f'sep_{a}_{b}'] = pos[a].separation(pos[b]).deg
    atom = {}
    for a, (_, tag, _) in ARMS.items():  # the corpus must come from these decompositions: every line it holds
        with fits.open(dcfg.coef_fits(tag), memmap=True) as h:
            atom[tag] = np.asarray(h['COEF'].data['ATOM_Or'], float)
            if atom[tag].size != expnum.size:
                raise ValueError(f'{dcfg.coef_fits(tag)} does not hold one row per exposure of the stack')
            for column, line in iono_corpus.LINES.items():
                c, v = np.asarray(h['COEF'].data[column], float)[file_row], ic[f'{a}_coef_{line}']
                ok = np.isfinite(v)  # the corpus leaves a few rows empty
                if not np.allclose(c[ok], v[ok], rtol=1e-5):
                    raise ValueError(f'{corpus} does not hold the {line} coefficients of {dcfg.coef_fits(tag)}: '
                                     'rebuild the corpus')
    ctx, ctx_names = _context(dcfg, ic['obstime_mjd'])
    for a in ARMS:
        tab |= {f'{a}_ctx_{c}': ctx[a][:, i] for i, c in enumerate(ctx_names)}
    in_region = lambda reg: data._angular_separation_deg_vec(  # noqa: E731
        col('sci_ra'), col('sci_dec'), reg['ra_deg'], reg['dec_deg']) <= reg['radius_deg']
    aux = {'ctx_names': ctx_names, 'file_row': file_row, 'gal_b': pos['sci'].galactic.b.deg,
           'valid_coords': ~np.any([(col(f'{p}_ra') == 0) & (col(f'{p}_dec') == 0) for p in POINTINGS], axis=0),
           'lmc': in_region(data.LMC_EXCLUSION), 'smc': in_region(data.SMC_EXCLUSION),
           'n_fibres': {a: col(f'fibers_{p}_used') for a, (p, _, _) in ARMS.items()}}
    with fits.open(dcfg.input_fits_flux, memmap=True) as stack, np.errstate(divide='ignore', invalid='ignore'):
        if not np.array_equal(stack['META'].data['expnum'], expnum):
            raise ValueError('the median stack and its meta_only file list different exposures')
        wave = np.asarray(stack['WAVE'].data, float)
        aux['log_ew'] = {a: np.log10(oi6300.equivalent_width(*_window(stack[h], file_row, wave, 6200, 6400),
                                                             ic[f'{a}_coef_6300'])) for a, (_, _, h) in ARMS.items()}
        flux = {a: _window(stack[h], file_row, wave, *HALPHA_WINDOW)[1] for a, (_, _, h) in ARMS.items()}
    # Halpha left in the decomposition residual (it has no Halpha template); SkyNear removes the geocoronal part.
    # RESID is in units of 1 / dcfg.factor erg s-1 cm-2 A-1; the sum over pixels times the pixel width is a flux.
    photons = HALPHA * 1e-8 / oi6300.HC * 4 * np.pi * 1e-6 / oi6300.OMEGA  # R per erg s-1 cm-2 per fibre
    to_rayleigh = photons / dcfg.factor * float(np.median(np.diff(wave)))
    ha = {}
    for a, (_, tag, _) in ARMS.items():
        # The decomposition files carry no exposure number.  Their rows and grid are those of the stack if
        # residual + best fit gives back its flux, and the coefficient file is the same fit if ATOM_Or agrees.
        with fits.open(dcfg.decomp_fits(tag), memmap=True) as h:
            ok = h['RESID'].data.shape == (expnum.size, wave.size) and np.array_equal(
                np.asarray(h['COEF'].data['ATOM_Or'], float), atom[tag], equal_nan=True)
            resid = _window(h['RESID'], file_row, wave, *HALPHA_WINDOW)[1] if ok else None
            total = resid + _window(h['BESTFIT'], file_row, wave, *HALPHA_WINDOW)[1] if ok else None
            fitted = np.isfinite(total)  # a failed fit leaves its row empty
            if not ok or not np.allclose(total[fitted], dcfg.factor * flux[a][fitted], rtol=1e-6):
                raise ValueError(f'{dcfg.decomp_fits(tag)} does not match the rows or the grid of the stack, or the '
                                 f'fit of {dcfg.coef_fits(tag)}')
            ha[a] = to_rayleigh * resid.sum(1)
    aux['halpha_excess'] = ha['sci'] - ha['near']
    return tab, aux


def row_sets(tab, aux):
    """Boolean row masks: ``clean``, ``lmcsmc``, ``domain``, ``gaspoor`` and ``train`` = ``lmcsmc | gaspoor``."""
    clean = oi6300.clean_rows(tab, aux['n_fibres'], aux['log_ew']) & aux['valid_coords']
    clouds = aux['lmc'] | aux['smc']  # their own [O I] is redshifted off the atmospheric template
    gaspoor = (clean & ~clouds & (aux['halpha_excess'] < GAS_POOR_MAX_HALPHA)
               & (np.abs(aux['gal_b']) > GAS_POOR_MIN_ABS_B))
    return {'clean': clean, 'lmcsmc': clean & clouds, 'domain': clean & ~clouds, 'gaspoor': gaspoor,
            'train': (clean & clouds) | gaspoor}


def control_lines(dcfg, tab, aux, rows, oh_groups=None):
    """Add ``sci_near_log_oh`` and ``sci_near_log_o2`` to ``tab``; returns the OH groups used.

    ``oh_groups=None`` takes the ``N_OH`` groups with the largest median SkyNear coefficient on ``rows``; pass the
    ``oh_groups`` of a saved model to build the inputs it was trained with.
    """
    coef = {a: fits.getdata(dcfg.coef_fits(ARMS[a][1]), 'COEF') for a in ('sci', 'near')}
    grab = lambda a, cols: np.column_stack([np.asarray(coef[a][c], float)[aux['file_row']] for c in cols])  # noqa: E731
    if oh_groups is None:
        oh = [n for n in coef['near'].names if n.startswith('OH_')]
        oh_groups = [oh[i] for i in np.argsort(-np.nanmedian(grab('near', oh)[rows], axis=0))[:N_OH]]
    v = tab['sci_vanrhijn_95km'], tab['near_vanrhijn_95km']
    tab['sci_near_log_oh'] = oi6300.control_ratio(grab('sci', oh_groups), grab('near', oh_groups), *v)
    tab['sci_near_log_o2'] = oi6300.control_ratio(grab('sci', [O2_COLUMN]), grab('near', [O2_COLUMN]), *v)
    return tuple(oh_groups)


def night_field(tab, aux, hyper=None, fit_rows=None):
    """Add the ``oi6300.SKY_FIELD`` columns to ``tab``; returns the kernel parameters used.

    The field of each night rests on every exposure whose SkyNear and SkyFar medians hold enough fibres
    (``oi6300.CLEAN_MIN_FIBRES``) and whose pointings are valid.  ``hyper=None`` learns the kernel parameters from
    the sky telescopes of ``tab`` (``sky_field.fit``; only from the rows ``fit_rows`` if given); pass the
    ``field_hyper`` of a saved model to build the inputs it was trained with.
    """
    usable = np.all([aux['n_fibres'][a] >= oi6300.CLEAN_MIN_FIBRES for a in ('near', 'far')], axis=0)
    usable &= aux['valid_coords']  # a pointing at ra = dec = 0 puts the sample at a wrong pierce point
    hyper = hyper or sky_field.fit(tab, usable if fit_rows is None else usable & fit_rows)
    tab |= sky_field.night_field(tab, usable, hyper)
    return hyper


def night_folds(tab, n=N_FOLDS):
    """Group index of every row: the corpus nights are split at random into ``n`` groups."""
    nights = np.unique(tab['night_mjd'])
    group = dict(zip(np.random.default_rng(0).permutation(nights), np.arange(nights.size) % n, strict=True))
    return np.array([group[x] for x in tab['night_mjd']])


def cross_validate(tab, rows, domain, members=1, prepare=None, **model_kw):
    """Test-night prediction of log10 I_sci for ``rows`` (NaN elsewhere), and whether a guarded line was not used.

    Each of the ``N_FOLDS`` runs trains the whole model (with ``members`` > 1 the whole ensemble) without one group
    of nights and predicts that group, so no row is predicted by a model that saw its night.  Nothing else is
    learned from the test nights either: the coverage cut reads the ``domain`` rows of the other nights, and
    ``prepare(seen)``, if given, rebuilds the inputs that are learned from the corpus (OH groups, kernel of the
    night field) from the rows ``seen`` of the other nights.  ``prepare`` leaves ``tab`` as built for the last run.
    """
    fold, rows = night_folds(tab), np.flatnonzero(rows)
    pred, dropped = np.full(fold.size, np.nan), np.zeros(fold.size, bool)
    for k in range(N_FOLDS):
        if prepare is not None:
            prepare(fold != k)
        tr, te = rows[fold[rows] != k], rows[fold[rows] == k]
        m = oi6300.OI6300Model(**model_kw).fit(tab, tr, domain=np.asarray(domain) & (fold != k), members=members)
        pred[te] = m.predict(tab)[0][te]
        dropped[te] = ~m.lines_used(tab)[te].all(1)
    return pred, dropped


def scores(tab, pred, sets):
    """Test-night scatter [dex] of log10 I_obs - log10 I_pred on the training rows and its subsets."""
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.log10(oi6300.K['6300'] * np.asarray(tab['sci_coef_6300'], float)) - pred
    sep = tab['sep_sci_near']
    groups = {'all': sets['train'], 'LMC+SMC': sets['lmcsmc'], 'gas-poor': sets['gaspoor']}
    groups |= {f'separation {lo}-{hi} deg': sets['train'] & (sep >= lo) & (sep < hi) for lo, hi in SEP_BINS}
    rms = lambda x: float(np.sqrt(np.mean(x**2)))  # noqa: E731
    return {k: {'rows': int(g.sum()), 'sigma': round(sigma(r[g]), 5), 'rms': round(rms(r[g]), 5),
                'bias': round(float(np.median(r[g])), 5), 'abs_p95': round(float(np.percentile(np.abs(r[g]), 95)), 5)}
            for k, g in groups.items()}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()[:16]


def _provenance(dcfg, corpus):
    """Identity of the input files, of the code and of the settings behind a model."""
    files = [dcfg.input_fits_flux, dcfg.input_fits_meta, corpus, *(dcfg.coef_fits(t) for _, t, _ in ARMS.values()),
             *(dcfg.decomp_fits(t) for _, t, _ in ARMS.values())]
    stat = lambda p: Path(p).resolve().stat()  # noqa: E731
    inputs = {Path(p).name: {'bytes': stat(p).st_size} | (  # the multi-GB spectral files: size and time only
        {'sha256_16': _sha256(p)} if stat(p).st_size < 1e9 else {'mtime': int(stat(p).st_mtime)}) for p in files}
    here = Path(__file__).parent
    git = subprocess.run(['git', '-C', str(here), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
    import astropy
    import scipy
    import sklearn

    return {'input_files': inputs, 'git_head': git, 'inputs_version': oi6300.INPUTS_VERSION,
            'code_sha256_16': {f: _sha256(here / f) for f in ('oi6300.py', 'oi6300_build.py', 'sky_field.py')},
            'settings': {'K_rayleigh_per_coefficient': oi6300.K, 'hgb': oi6300.HGB,
                         'clean': [oi6300.CLEAN_MIN_FIBRES, oi6300.CLEAN_MIN_LOG_EW],
                         'gas_poor': [GAS_POOR_MAX_HALPHA, GAS_POOR_MIN_ABS_B], 'n_oh': N_OH, 'folds': N_FOLDS,
                         'roti_floor': oi6300.ROTI_FLOOR},
            'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
                         'scikit-learn': sklearn.__version__, 'astropy': astropy.__version__}}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--decomp-root', required=True, help='directory of the median stack and its decompositions')
    p.add_argument('--stem', required=True, help='median-stack file stem (file name without .fits)')
    p.add_argument('--corpus', help='use this ionosphere corpus .npz and build no context')
    p.add_argument('--rebuild-context', action='store_true', help='rebuild the context tables and the corpus')
    p.add_argument('--workers', type=int, default=4, help='parallel downloads of GNSS station-days')
    p.add_argument('--out', required=True, help='output directory')
    p.add_argument('--name', help='output file stem (default: oi6300_model[_<lines>]_ens<members>)')
    p.add_argument('--sci-lines', nargs='*', default=list(DEFAULT_SCI_LINES),
                   choices=oi6300.SCI_CONTROL + oi6300.SCI_GUARDED,
                   help='Science-spectrum airglow lines used as inputs (default: %(default)s; no value: none)')
    p.add_argument('--sky-field', action=argparse.BooleanOptionalAction, default=True,
                   help='night field of the sky telescopes among the inputs (default: on)')
    p.add_argument('--guard', type=float, default=3.0, help='drop a guarded Science line above this many sigma')
    p.add_argument('--members', type=int, default=5, help='ensemble size (1 = single model)')
    p.add_argument('--no-cv', action='store_true', help='skip the cross-validation (no accuracy in the summary)')
    args = p.parse_args(argv)
    if args.members < 1 or len(set(args.sci_lines)) != len(args.sci_lines):
        p.error('--members must be >= 1 and --sci-lines must not repeat a line')

    dcfg = PipelineConfig().data
    dcfg.decomp_data_root, dcfg.decomp_stem = args.decomp_root, args.stem
    corpus = Path(args.corpus) if args.corpus else iono_corpus.build_context(
        args.decomp_root, args.stem, args.workers, args.rebuild_context)['corpus']
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    name = args.name or '_'.join(
        ['oi6300_model', *args.sci_lines, *['field'] * args.sky_field, f'ens{args.members}'])

    tab, aux = build_table(dcfg, corpus)
    sets = row_sets(tab, aux)
    controls = bool(set(args.sci_lines) & set(oi6300.SCI_CONTROL))

    def prepare(seen=True):
        """Inputs learned from the corpus, from the rows ``seen`` only: OH groups, kernel of the night field."""
        oh = control_lines(dcfg, tab, aux, sets['train'] & seen) if controls else ()
        return oh, (night_field(tab, aux, fit_rows=None if seen is True else seen) if args.sky_field else None)

    for k, v in sets.items():
        print(f'{k:8s} {v.sum():6d} rows, {np.unique(tab["night_mjd"][v]).size:4d} nights', flush=True)
    kw = dict(ctx_names=aux['ctx_names'], sci_lines=tuple(args.sci_lines), guard=args.guard, sky_field=args.sky_field)
    summary = {'name': name, 'created': time.strftime('%Y-%m-%d %H:%M:%S%z'), 'sci_lines': args.sci_lines,
               'sky_field': args.sky_field, 'guard_sigma': args.guard, 'members': args.members,
               'decomp_root': args.decomp_root, 'stem': args.stem,
               'decomp_suffix': dcfg.decomp_suffix, 'corpus': corpus.name, 'corpus_sha256_16': _sha256(corpus),
               'rows': {k: int(v.sum()) for k, v in sets.items()},
               'rows_dropped_for_zero_coordinates': int((~aux['valid_coords']).sum())}
    pred = np.full(sets['train'].size, np.nan)
    if not args.no_cv:
        pred, dropped = cross_validate(tab, sets['train'], sets['domain'], members=args.members, prepare=prepare, **kw)
        tr = sets['train']
        summary['test_nights'] = scores(tab, pred, sets) | {
            'rows with a dropped line': round(float(dropped[tr].mean()), 4)}
        for k, v in summary['test_nights'].items():
            print(f'test nights, {k}: {v}', flush=True)
        y = np.log10(oi6300.K['6300'] * np.asarray(tab['sci_coef_6300'], float)[tr])
        np.savetxt(out / f'{name}_cv.csv', np.column_stack(
            [tab['expnum'][tr], tab['night_mjd'][tr], night_folds(tab)[tr], y, pred[tr], dropped[tr]]),
            fmt=['%d', '%d', '%d', '%.5f', '%.5f', '%d'], delimiter=',', comments='',
            header='expnum,night_mjd,fold,log10_observed_R,log10_predicted_R,line_dropped')

    oh_groups, field_hyper = prepare()  # the saved model learns from every night
    m = oi6300.OI6300Model(**kw).fit(tab, sets['train'], domain=sets['domain'], members=args.members)
    m.oh_groups, m.field_hyper = oh_groups, field_hyper
    summary |= {'oh_groups': list(oh_groups), 'field_hyper': field_hyper, 'alpha': float(m.alpha),
                'inputs_kept': int(m.keep.size), 'inputs_total': len(m.names),
                'inputs_dropped': [n for i, n in enumerate(m.names) if i not in set(m.keep)],
                'line_level_sigma_dex': {n: round(float(s), 5) for n, (_, s) in m.level.items()}}
    summary |= _provenance(dcfg, corpus)
    # Per corpus exposure, for SkyML and for training the sky-transfer MLP: the night field, and the prediction a
    # model that did not see the night gives (the test-night run for the training rows, this model for the rest).
    cv = sets['train'] & np.isfinite(pred)
    m.corpus = {'expnum': np.asarray(tab['expnum']), 'obstime_mjd': np.asarray(tab['obstime_mjd'], float),
                'log10_pred': np.where(cv, pred, m.predict(tab)[0]), 'pred_from_cv': cv,
                **{c: np.asarray(tab[c], float) for c in oi6300.SKY_FIELD if args.sky_field}}
    summary['corpus_predictions'] = {'rows': int(cv.size), 'finite': int(np.isfinite(m.corpus['log10_pred']).sum()),
                                     'from_test_nights': int(cv.sum())}
    m.meta = summary
    m.save(out / f'{name}.pkl.part')  # the pair appears only when the build is complete
    (out / f'{name}.pkl.part').replace(out / f'{name}.pkl')
    (out / f'{name}.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(f'saved {out / name}.pkl ({(out / f"{name}.pkl").stat().st_size / 1e6:.0f} MB) and {name}.json')


if __name__ == '__main__':
    main()
