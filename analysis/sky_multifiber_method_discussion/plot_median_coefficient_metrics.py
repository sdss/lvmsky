"""Plot exact coefficient contractions for saved SkyE/SkyW median LSFs."""
from pathlib import Path
import csv
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from skysub.sky_decomp.result_io import lsf_surface_state_from_cubes
from skysub.sky_decomp.lsf_spline2d import (
    LSF_SPLINE2D_REPRESENTATION, _offset_knots, evaluate_bspline_basis,
    mspline_basis, mspline_bin_integrals,
)


def basis_moments(n_basis, half_width, order=4):
    """Integrate cubic basis times u**r exactly, piece by piece, for r <= 3."""
    knots = _offset_knots(n_basis, half_width)
    z, w = leggauss(order)
    mid, half = (knots[:-1] + knots[1:]) / 2, np.diff(knots) / 2
    u = (mid[:, None] + half[:, None] * z).ravel()
    weights = (half[:, None] * w).ravel()
    basis = mspline_basis(u, n_basis=n_basis, half_width=half_width)
    return np.stack([(weights * u**r) @ basis for r in range(4)])


def coefficient_metrics(state, channel, wavelength):
    assert state.legacy_kernel_representation == LSF_SPLINE2D_REPRESENTATION
    h = state.coefficients[channel]
    n, width = h.shape[0], state.config['offset_half_width_angstrom']
    basis = evaluate_bspline_basis(wavelength, state.knot_vectors[channel], state.degrees[channel])
    moments = basis_moments(n, width)
    np.testing.assert_allclose(moments, basis_moments(n, width, 8), atol=2e-14, rtol=2e-14)
    # Contract the offset basis before evaluating the wavelength basis.
    raw = basis @ (moments @ h).T
    np.testing.assert_allclose(raw[:, 0], 1, atol=2e-8, rtol=0)
    assert np.min(h) >= -1e-10 and np.isfinite(raw).all()
    mu, second, third = raw[:, 1], raw[:, 2], raw[:, 3]
    variance = second - mu**2
    assert np.all(variance > 0)
    central_third = third - 3 * mu * second + 2 * mu**3
    metrics = dict(centroid_A=mu, raw_second_A2=second, variance_A2=variance,
                   sigma_A=np.sqrt(variance), raw_third_A3=third,
                   central_third_A3=central_third, skewness=central_third / variance**1.5)
    for radius in (1.5, 2.0):
        assert radius < width
        intervals = mspline_bin_integrals(np.array([-width, radius]), np.array([-radius, width]),
                                         n_basis=n, half_width=width)
        tail = basis @ (intervals.sum(axis=0) @ h)
        assert np.all((tail >= -1e-10) & (tail <= 1 + 1e-10))
        metrics[f'tail_outside_{radius:g}A'] = tail
    # Independent dense-profile integration at three wavelengths checks contractions.
    u = np.linspace(-width, width, 20001)
    profiles = (basis[[0, len(basis)//2, -1]] @ h.T) @ mspline_basis(u, n_basis=n, half_width=width).T
    dense = np.stack([np.trapz(profiles * u**r, u, axis=1) for r in range(4)], axis=1)
    np.testing.assert_allclose(dense, raw[[0, len(basis)//2, -1]], atol=3e-10, rtol=1e-8)
    return metrics


def main():
    out = Path(__file__).resolve().parent
    manifest = json.loads((ROOT/'analysis/sky_subtraction_validation/run_manifest.json').read_text())
    source = next(r for r in manifest if r['expnum'] == 59750)['output']
    records, summaries = [], []
    with fits.open(source, memmap=True) as hdus:
        wave = np.asarray(hdus['WAVE'].data, float)
        labels = list(hdus['MED_LABELS'].data['telescope'])
        cubes = tuple(hdus['MED_'+name].data for name in ('LSF_COEF','LSF_KNOTS','LSF_META'))
        mjd = hdus[0].header.get('MJD')
        for telescope in ('SkyE', 'SkyW'):
            state = lsf_surface_state_from_cubes(cubes, labels.index(telescope))
            fig, axes = plt.subplots(5, 3, figsize=(14.5, 12), sharex='col', sharey='row')
            for col, channel in enumerate(('B','R','Z')):
                knots, degree = state.knot_vectors[channel], state.degrees[channel]
                x = wave[(wave >= knots[degree]) & (wave <= knots[-degree-1])]
                m = coefficient_metrics(state, channel, x)
                for j, wavelength in enumerate(x):
                    records.append(dict(exposure=59750, telescope=telescope, channel=channel,
                                        wavelength_A=wavelength, **{k: v[j] for k,v in m.items()}))
                axes[0,col].plot(x,m['centroid_A'],color='#22689b',label='First moment about zero')
                axes[0,col].axhline(0,color='#888888',lw=.6)
                axes[1,col].plot(x,m['raw_second_A2'],color='#22689b',label='Second moment about zero')
                axes[1,col].plot(x,m['variance_A2'],color='#d76c22',ls='--',label='Central variance')
                axes[2,col].plot(x,100*m['tail_outside_1.5A'],color='#22689b',label='Outside ±1.5 Å')
                axes[2,col].plot(x,100*m['tail_outside_2A'],color='#d76c22',ls='--',label='Outside ±2.0 Å')
                axes[3,col].plot(x,m['raw_third_A3'],color='#22689b',label='Third moment about zero')
                axes[3,col].plot(x,m['central_third_A3'],color='#d76c22',ls='--',label='Central third moment')
                axes[4,col].plot(x,m['skewness'],color='#22689b',label='Central standardized skewness')
                for row in (3,4): axes[row,col].axhline(0,color='#888888',lw=.6)
                for ax in axes[:,col]:
                    ax.grid(alpha=.18);ax.tick_params(labelsize=9);ax.margins(x=.01)
                axes[0,col].set_title(channel+(' · fixed survey reference' if channel=='B' else ' · fitted median LSF'),pad=12,fontsize=12)
                axes[4,col].set_xlabel('Air wavelength (Å)',fontsize=11)
                summaries.append(dict(telescope=telescope,channel=channel,
                    ranges={k:[float(v.min()),float(v.max())] for k,v in m.items()}))
            ylabels=['Centroid (Å)\nlinear in coefficients','Second moment (Å²)\nraw: linear; central: nonlinear',
                     'Tail mass (%)\nfixed origin; linear', 'Third moment (Å³)\nraw: linear; central: nonlinear',
                     'Skewness\ncentral; nonlinear']
            for row,label in enumerate(ylabels):
                axes[row,0].set_ylabel(label,fontsize=10)
                if row in (1,2,3):
                    axes[row,1].legend(loc='upper center',bbox_to_anchor=(.5,-.17),ncol=2,fontsize=9,frameon=False)
            fig.suptitle(f'Exposure 59750 · MJD {mjd} · {telescope} median LSF\nMetrics from saved coefficient contractions',fontsize=15,y=.99)
            fig.text(.5,.025,'Original fitted profiles; no recentering. Full stored support: ±2.5 Å. Tail thresholds are illustrative, not fitted constraints.',ha='center',fontsize=10)
            fig.subplots_adjust(left=.12,right=.985,top=.90,bottom=.105,hspace=.62,wspace=.16)
            stem=out/f'59750_{telescope}_median_lsf_coefficient_metrics'
            fig.savefig(stem.with_suffix('.png'),dpi=150)
            fig.savefig(stem.with_suffix('.pdf'))
            plt.close(fig)
    with (out/'59750_median_lsf_coefficient_metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    (out/'59750_median_lsf_coefficient_metrics.json').write_text(json.dumps(dict(
        source=source,exposure=59750,mjd=mjd,normalization_tolerance=2e-8,
        method='Exact piecewise Gauss integration of cubic M-splines through order 3; analytic tail integrals; coefficient contractions before wavelength evaluation.',
        checks='4-point versus 8-point Gauss quadrature; independent dense-profile moment integration at three wavelengths per channel; positivity, normalization and tail bounds.',
        profiles_recentered=False,tail_origin='Stored offset zero',summary=summaries),indent=2))
    print(json.dumps(summaries,indent=2))


if __name__ == '__main__':
    main()
