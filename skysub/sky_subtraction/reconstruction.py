"""Render SkyML spectra and fit diagnostics from compact, recorded parameters; no fitting or ML."""

import json
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from functools import cached_property, lru_cache
from pathlib import Path

import numpy as np
from astropy.io import fits
from lvmdrp.core.fluxcal import TelluricCalculator
from scipy.interpolate import BSpline

from skysub.sky_decomp.lsf_reference import load_survey_ascii_lsf
from skysub.sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig
from skysub.sky_decomp.result_io import lsf_surface_state_from_cubes
from skysub.sky_decomp.solar_lsf import SkyDecompPalaceCorrSolarBlueLSF

from .fiber_lsf import sky_designs
from .io import BUNDLED_CALIBRATION, read_cframe
from .product import read_json_hdu, sha256
from .rendering import (
    FACTOR,
    FAMILIES,
    SkyArm,
    channel_slices,
    exposure_line_masks,
    fit_pixels,
    fit_variance,
    huber_factor,
    render,
    sky_parts,
)
from .settings import Settings

TEMPLATES = Path(__file__).resolve().parents[1] / 'sky_decomp/data'
FIT_KEYS = ('FIT_MODEL', 'FIT_SKY', 'FIT_LINES', 'FIT_SOLAR', 'FIT_CONT', 'FIT_IVAR')
DIAGNOSTIC_JOBS = 8


@lru_cache(maxsize=16)
def _input_hash(path, size, mtime):
    return sha256(path)


class SkyMLReconstructor:
    """One compact product plus its original CFrame; no fitting or ML inference."""

    def __init__(self, path, cframe=None):
        self.h = fits.open(path, memmap=True)
        try:
            if self.h[0].header.get('SKYMFMT') != 'COEFFICIENTS':
                raise ValueError('Expected a compact coefficient-based SkyML SFrame')
            self.provenance = read_json_hdu(self.h['SKYML_META'])
            self.settings = Settings.from_record(self.provenance.get('settings'))
            self.source = Path(cframe or self.provenance['input'])
            self._observed_input = None
            self.wave = np.asarray(self.h['WAVE'].data, float)
            self.context = read_json_hdu(self.h['RENDER_META'])
            self.labels = list(self.h['MED_LABELS'].data['telescope'])
            self.models = {}
            self.medians = {label: self._median(i, label) for i, label in enumerate(self.labels)}
            self.model_sci = self.model('Sci')
            ml = self.h['ML_COEF'].data
            self.coef = np.asarray(ml['coef'], float)
            self.oh_coef = np.asarray(ml['oh_coef'], float)
            if list(ml['name']) != list(self.model_sci.design_names):
                raise ValueError('Stored coefficient names differ from renderer schema')
            science = self.labels.index('Sci')
            _, _, matrices = render(self.model_sci, self.coef, self.state('MED_', science))
            self.designs = sky_designs(self.model_sci, self.oh_coef, self.state('INIT_', science), matrices)
            self.arms = [self._sky_arm(i, label) for i, label in enumerate(self.labels) if label != 'Sci']
            velocities = self.provenance.get('residual_correction', {}).get('line_protection_velocities_km_s', [])
            self.lines = exposure_line_masks(self.wave, self.state('MED_', science), velocities, self.settings)
        except BaseException:
            self.h.close()
            raise

    def _median(self, index, label):
        total, components, _ = render(self.model(label), self.h['MED_COEF'].data[index], self.state('MED_', index))
        return total / FACTOR, {k: v / FACTOR for k, v in components.items()}

    def _sky_arm(self, index, label):
        total, components = self.medians[label]
        weight = np.asarray(self.h['MED_IVAR'].data[index], float)
        variance = np.divide(1, weight, out=np.full_like(weight, np.inf), where=weight > 0)
        return SkyArm(self.h['MED_FLUX'].data[index], total, components['moon'] + components['zodi'], variance)

    def _read_observed(self):
        if self._observed_input is None:
            stat = self.source.stat()
            if _input_hash(str(self.source), stat.st_size, stat.st_mtime_ns) != self.provenance['input_sha256']:
                raise ValueError('CFrame hash differs from recorded input')
            frame = read_cframe(self.source)
            if not np.array_equal(frame.slit['fiberid'], self.h['SLITMAP'].data['fiberid']):
                raise ValueError('CFrame SLITMAP fiber IDs differ from the product')
            if not np.array_equal(frame.wave, self.wave):
                raise ValueError('CFrame wavelength grid differs from the product')
            self._observed_input = (frame.flux, frame.ivar, frame.mask)
        return self._observed_input

    @property
    def observed(self):
        return self._read_observed()[0]

    @property
    def ivar(self):
        return self._read_observed()[1]

    @property
    def mask(self):
        return self._read_observed()[2]

    def close(self):
        self.h.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def state(self, prefix, index):
        cubes = tuple(self.h[prefix + name].data for name in ('LSF_COEF', 'LSF_KNOTS', 'LSF_META'))
        return lsf_surface_state_from_cubes(cubes, index)

    def model(self, label):
        if label not in self.models:
            i = next(i for i, m in enumerate(self.context['models']) if m['telescope'] == label)
            meta = self.context['models'][i]
            kw = dict(self.context['settings'])
            kw['config'] = LSFSurfaceIterativeConfig(**kw['config'])
            # Templates ship with this codebase; the stored path is provenance only.
            kw['base_dir'] = TEMPLATES
            model = SkyDecompPalaceCorrSolarBlueLSF(
                self.wave,
                telluric_calculator=TelluricCalculator(),
                pwv_mm=meta['pwv_mm'],
                source_airmass=meta['source_airmass'],
                drp_transmission=self.h['RENDER_TRANSMISSION'].data[i],
                **kw,
            )
            model.t_o2 = meta['t_o2']
            self.models[label] = model
        return self.models[label]

    @cached_property
    def calibration(self):
        """Survey LSF directory of the run; products before it was recorded used the bundled one."""
        directory = Path(self.provenance.get('calibration') or BUNDLED_CALIBRATION)
        recorded = self.provenance.get('calibration_manifest')
        if recorded is not None and json.loads((directory / 'manifest.json').read_text()) != recorded:
            raise ValueError(f'Survey LSF calibration {directory} differs from the one recorded in the product')
        return directory

    def _row(self, fiberid):
        match = np.flatnonzero(self.h['SLITMAP'].data['fiberid'] == fiberid)
        if len(match) != 1:
            raise ValueError('Unknown physical fiber ID')
        return int(match[0])

    def sky(self, fiberid, survey=False):
        """Adopted sky and its parts; needs coefficients and model assets, not the CFrame."""
        row = self._row(fiberid)
        match = np.flatnonzero(self.h['LSF_FIBERS'].data['input_row'] == row)
        if len(match) != 1:
            raise ValueError('Fiber was not processed')
        if survey:
            state = load_survey_ascii_lsf(self.calibration, int(fiberid), allow_warning=True)
        else:
            state = self.state('', int(match[0]))
        parts = sky_parts(self.model_sci, self.coef, self.oh_coef, state, self.arms, self.settings, self.lines)
        out = dict(SKY_ML=parts.ml, SKY_OH=parts.oh, SKY_RES=parts.residual)
        total = parts.total
        if not survey:
            qa = self.h['LSF_QA'].data
            for record in qa[qa['input_row'] == row]:
                if record['reason'] == 'unusable_survey_fallback':
                    channel = str(record['channel'])
                    total[channel_slices(self.wave, state, (channel,))[channel]] = np.nan
            if not np.allclose(total * FACTOR, self.h['SKY'].data[row] * FACTOR, rtol=2e-5, atol=2e-5, equal_nan=True):
                raise ValueError('Reconstructed SKY differs from recorded SKY; check renderer/model assets')
        out['SKY'] = total
        basis = np.vstack([parts.matrices[k] for k in FAMILIES])
        out['SKY_EPI'] = (np.asarray(self.h['ML_COEF'].data['coef_std'], float) ** 2) @ (basis**2) / FACTOR**2
        return out

    def fiber(self, fiberid, iteration=None):
        """Observed-data diagnostics; optionally reconstruct a saved fit iteration."""
        out = self.sky(fiberid)
        row = self._row(fiberid)
        observed = self.observed[row]
        out['OBSFLUX'] = observed
        out['FLUX'] = observed - out['SKY']
        for name in FIT_KEYS:
            out[name] = np.full(len(self.wave), np.nan)
        out['FIT_MASK'] = np.ones(len(self.wave), dtype='u1')
        if 'FIT_LSF_HISTORY' not in self.h:
            return out
        history = self.h['FIT_LSF_HISTORY'].data
        ivar = np.where(self.mask[row] == 0, self.ivar[row], 0.0)
        fit_mask = self.h['LSF_FIT_MASK'].data.astype(bool)
        for channel, (take, design) in self.designs.items():
            selected = history[(history['input_row'] == row) & (history['channel'] == channel)]
            if iteration is not None:
                selected = selected[selected['iteration'] == iteration]
            if len(selected):
                self._fit_diagnostics(out, selected[-1], take, design, observed[take], ivar[take], fit_mask[take])
        return out

    def _fit_diagnostics(self, out, record, take, design, observed, ivar, fit_mask):
        """Spectral model and robust weights of one recorded R/Z fit iteration."""

        def vector(key):
            return np.asarray(record[key][: int(record[key + '_n'])], float)

        wave = self.wave[take]
        knots = vector('continuum_knots')
        continuum = BSpline(knots, vector('continuum_coef'), 3)(wave)
        variance_continuum = BSpline(knots, vector('variance_continuum_coef'), 3)(wave)
        lsf = np.asarray(record['coefficients'], float)
        lines = record['line_amplitude'] * design['lines'] @ lsf
        solar = record['solar_amplitude'] * design['solar'] @ lsf + design['fixed']
        model = lines + solar + continuum
        good = fit_pixels(observed, ivar, fit_mask)
        variance = fit_variance(ivar, good, variance_continuum / FACTOR, self.settings)
        z = (observed - model / FACTOR) * np.sqrt(np.divide(1, variance, out=np.zeros_like(variance), where=good))
        for key, value in (
            ('FIT_MODEL', model),
            ('FIT_SKY', lines + solar),
            ('FIT_LINES', lines),
            ('FIT_SOLAR', solar),
            ('FIT_CONT', continuum),
        ):
            out[key][take] = value / FACTOR
        out['FIT_IVAR'][take] = np.where(good, huber_factor(z, record['huber_scale'], self.settings) / variance, 0.0)
        out['FIT_MASK'][take] = ~good


_RENDERER = None


def _initialize_renderer(renderer):
    global _RENDERER
    _RENDERER = renderer


def _render_fiber(fiberid):
    return _RENDERER.fiber(int(fiberid))


def diagnostic_hdul(path, fiberids=None, survey_baseline=False):
    """Transient report view: derived spectra exist in memory, never in the FITS."""
    with SkyMLReconstructor(path) as r:
        product = fits.HDUList([x.copy() for x in r.h])
        shape = product['FLUX'].data.shape
        keys = ['SKY_ML', 'SKY_OH', 'SKY_RES', 'SKY_EPI', *FIT_KEYS]
        arrays = {key: np.full(shape, np.nan, dtype='f4') for key in keys}
        if survey_baseline:
            arrays['SKY_SURVEY'] = np.full(shape, np.nan, dtype='f4')
        arrays['FIT_MASK'] = np.ones(shape, dtype='u1')
        arrays['OBSFLUX'] = r.observed.astype('f4')
        ids = product['LSF_FIBERS'].data['fiberid'] if fiberids is None else fiberids
        pool = None
        if len(ids) > 2 * DIAGNOSTIC_JOBS:
            pool = ProcessPoolExecutor(
                DIAGNOSTIC_JOBS, mp_context=mp.get_context('fork'), initializer=_initialize_renderer, initargs=(r,)
            )
        try:
            spectra = pool.map(_render_fiber, ids, chunksize=8) if pool else map(r.fiber, ids)
            for fiberid, values in zip(ids, spectra, strict=True):
                row = r._row(fiberid)
                if survey_baseline:
                    values['SKY_SURVEY'] = r.sky(int(fiberid), survey=True)['SKY']
                for key in arrays:
                    if key != 'OBSFLUX':
                        arrays[key][row] = values[key]
        finally:
            if pool:
                pool.shutdown(cancel_futures=True)
        for key, array in arrays.items():
            product.append(fits.ImageHDU(array, name=key))
        models = np.asarray([r.medians[label][0] for label in r.labels])
        components = np.asarray([[r.medians[label][1][name] for name in FAMILIES] for label in r.labels])
        product.append(fits.ImageHDU(models, name='MED_MODEL'))
        product.append(fits.ImageHDU(components, name='MED_COMPONENTS'))
        product.append(fits.ImageHDU(product['MED_FLUX'].data - models, name='MED_RESID'))
        return product
