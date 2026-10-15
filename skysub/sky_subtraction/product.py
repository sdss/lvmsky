"""Assembly and atomic publication of the compact SkyML SFrame."""

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

from skysub.sky_decomp.result_io import build_lsf_hdus

from .io import FLUX_UNIT, flux_scale
from .rendering import FAMILIES, fwhm_grid

FORMAT_VERSION = 3
# DRP MASK is nonzero for unusable pixels; SkyML sets this value where its sky is invalid.
INVALID_SKY_MASK = 1
DRP_SKY_EXTENSIONS = ('SKY_EAST', 'SKY_EAST_IVAR', 'SKY_WEST', 'SKY_WEST_IVAR')
HISTORY_VECTORS = ('continuum_knots', 'continuum_coef', 'variance_continuum_coef')


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _json(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def json_hdu(name, value, **options):
    payload = json.dumps(value, default=_json, **options).encode()
    return fits.ImageHDU(np.frombuffer(payload, dtype='u1'), name=name)


def read_json_hdu(hdu):
    return json.loads(bytes(hdu.data).decode())


def history_table(records):
    """Pad explicit-length spline vectors so R/Z histories share one FITS table."""
    widths = {key: max(len(r[key]) for r in records) for key in HISTORY_VECTORS}
    rows = []
    for record in records:
        row = dict(record)
        for key in HISTORY_VECTORS:
            values = np.asarray(row[key], float)
            row[key + '_n'] = len(values)
            row[key] = np.pad(values, (0, widths[key] - len(values)), constant_values=np.nan)
        rows.append(row)
    return Table(rows=rows)


@dataclass
class MLCoefficients:
    names: list
    coef: np.ndarray  # common Science-IFU prediction; float32 ensemble output, stored as float32
    coef_std: np.ndarray  # ensemble spread, stored as float32
    oh_coef: np.ndarray  # after the median OH band correction


@dataclass
class ProductInputs:
    sky: np.ndarray  # (n_fiber, n_wave) subtracted sky; zero outside ``rows``
    rows: np.ndarray  # processed Science rows
    adopted: list  # LSF states per processed row
    candidates: list
    qa: list  # LSF_QA records
    history: list  # FIT_LSF_HISTORY records
    photometry: Table  # GAIA, including selection columns
    sources: Table | None
    medians: dict  # telescope -> io.Median
    median_fits: object  # decomp_adapter.MedianFits
    initial_states: list  # survey-mean median LSF states, in ``medians`` order
    stack_meta: Table
    ml: MLCoefficients
    fit_mask: np.ndarray  # science-line mask of the individual fiber fits


def _frame_hdus(input_path, inputs):
    with fits.open(input_path, memmap=False) as h:
        product = fits.HDUList([x.copy() for x in h if x.name not in DRP_SKY_EXTENSIONS])
    scale = flux_scale(product['FLUX'].header)
    product['FLUX'].data *= scale
    product['IVAR'].data /= scale**2
    product['FLUX'].header['BUNIT'] = str(FLUX_UNIT)
    product['IVAR'].header['BUNIT'] = str(1 / FLUX_UNIT**2)
    sky, rows = inputs.sky, inputs.rows
    processed = np.isin(np.arange(len(sky)), rows)
    product['FLUX'].data[rows] -= sky[rows]
    invalid = processed[:, None] & ~np.isfinite(sky)
    product['MASK'].data[invalid] |= INVALID_SKY_MASK
    product['IVAR'].data[invalid] = 0
    wave = np.asarray(product['WAVE'].data, float)
    for row, state in zip(rows, inputs.adopted, strict=True):
        product['LSF'].data[row] = fwhm_grid(state, wave)
    product['FLUX'].header['COMMENT'] = 'Science flux minus SKY; other rows retain observed CFrame flux'
    header = product[0].header
    header['SKYMETH'] = 'skyml'
    header['SKYMLVER'] = FORMAT_VERSION
    header['IVARERR'] = 'OBS_ONLY'
    header['COMMENT'] = 'IVAR retains observed variance; sky fitting correlations are not modeled'
    header['SKYMFMT'] = 'COEFFICIENTS'
    header['COMMENT'] = 'SKY_EPI is reconstructed ML ensemble variance, not full sky uncertainty'
    image = fits.ImageHDU(sky, name='SKY')
    image.header['BUNIT'] = str(FLUX_UNIT)
    product.append(image)
    slit = Table(product['SLITMAP'].data)
    slit['skyml_applied'] = processed
    product['SLITMAP'] = fits.BinTableHDU(slit, name='SLITMAP')
    return product, np.asarray(slit['fiberid'])


def _lsf_hdus(prefix, states):
    hdus = build_lsf_hdus(states)
    for hdu in hdus:
        hdu.name = prefix + hdu.name
    return hdus


def _model_hdus(inputs, fiberids):
    labels = list(inputs.medians)
    fits_by_label = inputs.median_fits.results
    ml = inputs.ml
    hdus = [fits.BinTableHDU(inputs.photometry, name='GAIA')]
    if inputs.sources is not None:
        hdus.append(fits.BinTableHDU(inputs.sources, name='GAIA_SOURCES'))
    hdus += [
        fits.BinTableHDU(Table({'telescope': labels}), name='MED_LABELS'),
        fits.BinTableHDU(inputs.stack_meta, name='MED_META'),
        fits.BinTableHDU(
            Table(
                rows=[
                    dict(
                        telescope=k,
                        fit_status=fits_by_label[k].fit_status,
                        reduced_chi2=fits_by_label[k].reduced_chi2,
                        t_o2=fits_by_label[k].t_o2,
                    )
                    for k in labels
                ]
            ),
            name='MED_QA',
        ),
    ]
    for name, field in (('FLUX', 'flux'), ('IVAR', 'ivar'), ('COUNT', 'count')):
        hdus.append(fits.ImageHDU(np.array([getattr(inputs.medians[k], field) for k in labels]), name='MED_' + name))
    hdus.append(
        fits.BinTableHDU(Table({'input_row': inputs.rows, 'fiberid': fiberids[inputs.rows]}), name='LSF_FIBERS')
    )
    hdus += _lsf_hdus('', inputs.adopted)
    hdus += _lsf_hdus('CAND_', inputs.candidates)
    hdus += _lsf_hdus('MED_', [fits_by_label[k].lsf_state for k in labels])
    hdus += [
        fits.BinTableHDU(Table(rows=inputs.qa), name='LSF_QA'),
        fits.ImageHDU(inputs.median_fits.fit_mask, name='MED_FIT_MASK'),
        fits.ImageHDU(np.array([fits_by_label[k].coef for k in labels]), name='MED_COEF'),
        fits.BinTableHDU(Table({'family': FAMILIES}), name='COMPONENT_NAMES'),
        fits.BinTableHDU(
            Table(
                {
                    'name': ml.names,
                    'coef': ml.coef.astype('f4'),
                    'coef_std': ml.coef_std.astype('f4'),
                    'oh_coef': ml.oh_coef,
                }
            ),
            name='ML_COEF',
        ),
    ]
    if inputs.history:
        hdus.append(fits.BinTableHDU(history_table(inputs.history), name='FIT_LSF_HISTORY'))
    hdus += [
        json_hdu('RENDER_META', inputs.median_fits.render_context, sort_keys=True),
        fits.ImageHDU(inputs.median_fits.transmission, name='RENDER_TRANSMISSION'),
        fits.ImageHDU(np.asarray(inputs.fit_mask, dtype='u1'), name='LSF_FIT_MASK'),
    ]
    hdus += _lsf_hdus('INIT_', inputs.initial_states)
    return hdus


def publish(product, output, mode):
    """Write beside ``output`` and hard-link: atomic, and refuses a destination created meanwhile."""
    output = Path(output)
    fd, temporary = tempfile.mkstemp(prefix='.skyml-', suffix='.fits', dir=output.parent)
    os.close(fd)
    try:
        product.writeto(temporary, overwrite=True, checksum=True)
        os.chmod(temporary, mode)
        os.link(temporary, output)
    finally:
        os.unlink(temporary)


def write_product(input_path, output, inputs, provenance):
    """Preserve CFrame rows; processed Science rows carry the subtraction and model validity."""
    product, fiberids = _frame_hdus(input_path, inputs)
    try:
        product.extend(_model_hdus(inputs, fiberids))
        product.append(json_hdu('SKYML_META', provenance, indent=2, allow_nan=True))
        publish(product, output, Path(input_path).stat().st_mode & 0o777)
    finally:
        product.close()
