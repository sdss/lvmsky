# %%
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord


SAS_BASE_DIR = os.environ.get('SAS_BASE_DIR', '/Users/ik52/obs/sas')
f_b = os.path.join(SAS_BASE_DIR, "sdsswork/lvm/spectro/redux/1.1.1/1027XX/1027983/60581/lvmFrame-b-00026500.fits")
f_r = os.path.join(SAS_BASE_DIR, "sdsswork/lvm/spectro/redux/1.1.1/1027XX/1027983/60581/lvmFrame-r-00026500.fits")
f_z = os.path.join(SAS_BASE_DIR, "sdsswork/lvm/spectro/redux/1.1.1/1027XX/1027983/60581/lvmFrame-z-00026500.fits")
fits.info(f_b)

N_STD = 12

%matplotlib widget

plt.close('all')

def create_sky(wave, slitmap, flux, wave_trace, skytype='SkyE'):
    msk = slitmap['telescope'] == skytype
    sky_spec = flux[msk]
    sky_mask = mask[msk]
    sky_wave_trace = wave_trace[msk]
    sky_waves = np.zeros_like(sky_spec)

    sky_out_2d = np.zeros((sky_spec.shape[0], wave.size))
    x = np.arange(sky_spec.shape[1])
    for i in range(sky_spec.shape[0]):
        current_waves = np.polyval(sky_wave_trace[i]['COEFF'][::-1], x)
        sky_out_2d[i] = np.interp(wave, current_waves, sky_spec[i])

    _, sky_out_1d, _ = sigma_clipped_stats(sky_out_2d, axis=0)

    return sky_out_1d


def get_nearby_sky(hdr, ra, dec):
    coord_current_spec = SkyCoord(ra=ra, dec=dec, unit='deg')
    coord_SkyE = SkyCoord(ra=hdr['TESKYERA'], dec=hdr['TESKYEDE'], unit='deg')
    coord_SkyW = SkyCoord(ra=hdr['TESKYWRA'], dec=hdr['TESKYWDE'], unit='deg')

    separation_Spec_SkyE = coord_current_spec.separation(coord_SkyE).deg
    separation_Spec_SkyW = coord_current_spec.separation(coord_SkyW).deg
    nearby_sky = "SkyE" if separation_Spec_SkyE < separation_Spec_SkyW else "SkyW"

    return nearby_sky

# Collect data from all stars and all channels and plot them
stack_stars = []
plt.close()
for istar in range(N_STD):
# for istar in [1, 3]:
    current_data = {}
    fig, ax = plt.subplots(figsize=(10, 4))
    for file, suffix in [(f_b, 'b'), (f_r, 'r'), (f_z, 'z')]:
        slitmap = Table.read(file, 'SLITMAP')
        hdr = fits.getheader(file)
        flux = fits.getdata(file, 'FLUX')
        mask = fits.getdata(file, 'MASK')
        ivar = fits.getdata(file, 'IVAR')
        wave_trace = fits.getdata(file, 'WAVE_TRACE')
        lsf_trace = fits.getdata(file, 'LSF_TRACE')
        exptime_sky = hdr['EXPTIME']

        # read infor from header
        std_info = {
            hdr[f"STD{i}FIB"].strip(): {
                'fib': hdr[f"STD{i}FIB"],
                'id': hdr[f"STD{i}ID"],
                'exp': hdr[f"STD{i}EXP"],
                'ra': hdr[f"STD{i}RA"],
                'de': hdr[f"STD{i}DE"],
                't0': hdr[f"STD{i}T0"],
                't1': hdr[f"STD{i}T1"],
                'acq': hdr[f"STD{i}ACQ"],
            } for i in range(1, N_STD+1)
        }

        msk_spec = slitmap['telescope'] == 'Spec'
        msk_skye = slitmap['telescope'] == 'SkyE'
        msk_skyw = slitmap['telescope'] == 'SkyW'

        idx_stars_in_Spec = np.argwhere( np.sum(np.isfinite( flux[msk_spec]), axis=1) != 0).flatten()

        orig_ifulabel = slitmap[msk_spec][idx_stars_in_Spec[istar]]['orig_ifulabel']
        exptime_std = std_info[orig_ifulabel]['exp']

        w_coef = wave_trace['COEFF'][msk_spec][idx_stars_in_Spec[istar]][::-1]
        lsf_coef = lsf_trace['COEFF'][msk_spec][idx_stars_in_Spec[istar]][::-1]
        wave = np.polyval(w_coef, np.arange(flux.shape[1]))
        lsf = np.polyval(lsf_coef, np.arange(flux.shape[1]))
        spec = flux[msk_spec][idx_stars_in_Spec[istar]] / exptime_std
        espec = 1.0 / np.sqrt(ivar[msk_spec][idx_stars_in_Spec[istar]]) / exptime_std
        cur_mask = mask[msk_spec][idx_stars_in_Spec[istar]]
        bad = cur_mask != 0
        nearby_sky = get_nearby_sky(hdr, std_info[orig_ifulabel]['ra'], std_info[orig_ifulabel]['de'])
        spec_sky = create_sky(wave, slitmap, flux, wave_trace, skytype=nearby_sky) / exptime_sky

        current_data[suffix] = dict(
            wave=wave,
            spec=spec-spec_sky,
            espec=espec,
            mask=cur_mask,
            lsf=lsf,
            sky=spec_sky,
            info_dict_all=std_info,
            info=std_info[orig_ifulabel],
            hdr=hdr,
        )

        # plotting
        ax.plot(wave, spec)
        ax.plot(wave, spec_sky)
        ax.plot(wave, espec)
        ax.plot(wave, spec - spec_sky)

        ax.set_ylabel("counts / s")
        ax.set_title(f"Star {orig_ifulabel} {std_info[orig_ifulabel]['id']} {std_info[orig_ifulabel]['ra']} {std_info[orig_ifulabel]['de']} {nearby_sky}")
    plt.tight_layout()
    plt.show()

    stack_stars.append(current_data)



# %%
fits.info(file)


# %% [markdown]
# # Stellar fitting based on FBS
# %%
from time import process_time
from time import perf_counter
from astropy.io import fits
import numpy as np
from astropy.table import Table, vstack
from scipy import linalg
from scipy.spatial import Delaunay
from tqdm import tqdm


def read_stellar_templates(filenames, verbose=True, use_model="ap"):
    """
    Read stellar templates from BOSZ grid.
    
    Excludes duplicated models for the same parameters, keeping only the
    specified model type (default is 'ap' for ATLAS9) in the output grid.
    """
    t0 = process_time()
    l_temps = []
    l_pars = []
    for idx, filename in enumerate(filenames):
        if verbose:
            print(filename)

        temps_current, hdr = fits.getdata(filename, header=True)
        pars_current = Table.read(filename, 'PARAMS')
        pars_current['idx'] = np.arange(len(pars_current))

        # Find duplicates and select only one model type
        mask = np.ones(len(pars_current), dtype=bool)

        pars_grouped = pars_current.group_by(['teff', 'logg', 'metal', 'alpha'])
        for grp in pars_grouped.groups:
            if len(grp) > 1:
                msk = grp['model_type'] != use_model
                mask[grp[msk]['idx']] = False

        l_temps.append(temps_current[mask])
        l_pars.append(pars_current[mask])

    waves = fits.getdata(filenames[0], 'WAVE')

    if verbose:
        print("Concatenating grid templates...")
    temps = np.vstack(l_temps)
    pars = vstack(l_pars)

    if verbose:
        print("Creating Delaunay triangulation...")

    tri = Delaunay(np.vstack((np.log10(pars['teff']),
                              pars['logg'],
                              pars['metal'],
                              pars['alpha'])).T)
    if verbose:
        print('Elapsed time in templates preparing: %.3f s' %
              (process_time() - t0))

    return waves, temps, pars, tri


STLIB_DIR = "/Users/ik52/sci/work/binary_stars/data/spec_libs/BOSZ/"

stellar_wave, stellar_templates, stellar_pars, stellar_tri = read_stellar_templates([
    STLIB_DIR + "bosz2024_r20000_m+0.50_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m+0.25_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m+0.00_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m-0.25_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m-0.50_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m-0.75_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m-1.00_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m-1.25_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m-1.50_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m-1.75_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
    STLIB_DIR + "bosz2024_r20000_m-2.00_3000_11000_Teff5000_8500_carbon0_vturb1.fits",
])

plt.close('all')
plt.plot(stellar_pars['teff'], stellar_pars['logg'], '.')
plt.xlabel('Teff')
plt.ylabel('logg')
plt.title('BOSZ Grid')
plt.show()



# %%
stellar_pars



# %%
def interp_rbfsimplex(data, tri, points):
    # find closes N+1 vertices of the N-dim triangle including the point
    npoints = points.shape[0]
    ndim = points.shape[1]
    out_spectra = np.zeros((data.shape[1], npoints))
    for j in range(npoints):
        point = points[j]
        current_simplex = tri.find_simplex(point)
        if current_simplex == -1:
            dst = linalg.norm(tri.points - point, axis=1)
            idx_n = dst.argsort()[:ndim+1]
        else:
            # these are incedes of vertices for simplex included xyz point
            idx_n = tri.simplices[current_simplex]
        simplex_vertices = tri.points[idx_n]

        distances = linalg.norm(simplex_vertices - point, axis=1)

        idx_zero = np.argwhere(distances == 0.0)

        if idx_zero.size == 0:
            weights = 1.0/distances**2 / np.sum(1.0/distances**2)
            # This stupid syntax provides best perfomrans (need re-check again!)
            if ndim == 5:
                out_spectra[:, j] = \
                    weights[0]*data[idx_n[0]] + \
                    weights[1]*data[idx_n[1]] + \
                    weights[2]*data[idx_n[2]] + \
                    weights[3]*data[idx_n[3]] + \
                    weights[4]*data[idx_n[4]] + \
                    weights[5]*data[idx_n[5]]
            elif ndim == 4:
                out_spectra[:, j] = \
                    weights[0]*data[idx_n[0]] + \
                    weights[1]*data[idx_n[1]] + \
                    weights[2]*data[idx_n[2]] + \
                    weights[3]*data[idx_n[3]] + \
                    weights[4]*data[idx_n[4]]
            elif ndim == 3:
                out_spectra[:, j] = \
                    weights[0]*data[idx_n[0]] + \
                    weights[1]*data[idx_n[1]] + \
                    weights[2]*data[idx_n[2]] + \
                    weights[3]*data[idx_n[3]]
            elif ndim == 2:
                out_spectra[:, j] = \
                    weights[0]*data[idx_n[0]] + \
                    weights[1]*data[idx_n[1]] + \
                    weights[2]*data[idx_n[2]]
            elif ndim == 1:
                raise ValueError("ndim = 1 is not supported")
            else:
                raise ValueError("ndim > 3 is not supported")
        else:
            # if point exactly equal to on the vertices. In this case just
            # return model spectrum from the vertex.
            out_spectra[:, j] = np.copy(data[idx_n[idx_zero[0][0]]])

    return out_spectra


msk = (stellar_pars['teff'] == 5500) & (stellar_pars['logg'] == 4) & (stellar_pars['metal'] == -0.0) & (stellar_pars['alpha'] == 0.0)
spec_ref = stellar_templates[msk].flatten()
spec_int = interp_rbfsimplex(stellar_templates, stellar_tri, np.array([[np.log10(5600), 4.0, 0, 0]]))[:, 0]

plt.close()
plt.plot(stellar_wave, spec_int, alpha=0.7, lw=0.7, label='interpolated')
plt.plot(stellar_wave, spec_ref, alpha=0.7, lw=0.7, label='reference')
plt.plot(stellar_wave, spec_ref-spec_int, lw=0.7, label="ref-int")
plt.legend()
plt.show()



# %%
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix

%load_ext line_profiler

def edges_from_centers(centers):
    centers = np.asarray(centers, dtype=float)
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5*(centers[1:] + centers[:-1])
    edges[0]  = centers[0]  - 0.5*(centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5*(centers[-1] - centers[-2])
    return edges


def fluxconserve_rebin(output_wave, input_wave, input_flux, normalize=True):
    """
    Flux-conserving rebin using linear inteprolation of the commulative function.
    """
    edges_out = edges_from_centers(output_wave)
    edges_in = edges_from_centers(input_wave)
    cdf_in = np.r_[0.0, np.cumsum(input_flux)]
    cdf_out = np.interp(edges_out, edges_in, cdf_in)

    if normalize:
        # in this case output spectrum will be on the same level of the input
        x_edges_in = np.arange(len(edges_in))
        x_edges_out = np.interp(edges_out, edges_in, x_edges_in)
        norm_factors = np.diff(x_edges_out)
    else:
        # spectrum level will be higer because of integration over pixels
        norm_factors = np.diff(edges_out)

    output = np.diff(cdf_out) / norm_factors

    return output


def convolve_and_rebin(input_wave, input_flux, output_wave, lsf_fwhm, pad_pix=30):
    """
    Function to convolve spectrum using sparse matrix multiplication, then flux conservative rebinning.

    Assume here that lsf is FWHM in Angstroms
    https://sdss-wiki.atlassian.net/wiki/spaces/LVM/pages/627048457/LVM+DRP+User+Guide#lvmFrame
    """

    # calculate range of pixels in the original spectrum covered by the output spectrum
    idx_left = np.argmax(input_wave - output_wave[0] > 0)
    idx_right = np.argmin(input_wave - output_wave[-1] < 0)
    pix_range = idx_left, idx_right
    pix_range_extended = idx_left - pad_pix, idx_right + pad_pix

    wave_sub = input_wave[pix_range[0]:pix_range[1]]
    wave_sub_extended = input_wave[pix_range_extended[0]:pix_range_extended[1]]

    # here assume that minimal delta lambda is on the blue part of the spectrum
    min_delta_ang = wave_sub_extended[1] - wave_sub_extended[0]

    # interpolate lsf to the input wave (stellar models wavelength)
    # interpolating outside of the range will return edge values
    lsf_sigma = np.interp(input_wave[pix_range[0]:pix_range[1]], output_wave, lsf_fwhm / 2.355)

    max_sigma = np.max(lsf_sigma)
    npix_half = 5 * np.ceil(max_sigma / min_delta_ang).astype(int)
    npix = 2 * npix_half + 1 # full size of the kernel

    x_px = np.linspace(-npix_half, npix_half, npix).astype(float) # X axis in pixels
    s_px = lsf_sigma / min_delta_ang # gaussian sigma in pixels

    yy = x_px[None, :] / s_px[:, None]
    kernels = np.exp(-0.5 * yy**2) / np.sqrt(2 * np.pi) / s_px[:, None]

    x_size = kernels.shape[0] # number of pixels in the NON-extenede range
    x_size_extended = wave_sub_extended.size

    rows2d = np.broadcast_to(np.arange(x_size)[:, None], (x_size, npix))
    cols2d = (pix_range[0] - pix_range_extended[0]) + rows2d + x_px.astype(int)

    data = kernels.ravel()
    rows = ((pix_range[0] - pix_range_extended[0]) + np.arange(x_size)[:, None] + x_px.astype(int)).ravel()
    cols = np.arange(0, x_size * npix + 1, npix)
    matrix_csr = csr_matrix((data, rows, cols), shape=(x_size, x_size_extended))

    conv_flux = matrix_csr @ input_flux[pix_range_extended[0]:pix_range_extended[1]]

    conv_rebin_flux = fluxconserve_rebin(output_wave, wave_sub, conv_flux)

    return conv_rebin_flux


spec_model = interp_rbfsimplex(stellar_templates, stellar_tri, np.array([[np.log10(5500), 4.3, -0.0, 0.00]]))[:, 0] / 1.2e4

spec_model_rebinned = fluxconserve_rebin(wave, stellar_wave, spec_model)

# %lprun -f convolve_and_rebin spec_conv = convolve_and_rebin(stellar_wave, spec_model, wave, lsf)

spec_conv = convolve_and_rebin(stellar_wave, spec_model, wave, lsf)

plt.close('all')
plt.step(wave, spec-spec_sky, where='mid', alpha=0.7, lw=0.7, color='C0')
plt.plot(wave, spec-spec_sky, '.', ms=5, alpha=0.7, lw=0.7, color='C0')
plt.step(stellar_wave, spec_model, where='mid', alpha=0.7, lw=0.7, color='C1')
plt.plot(stellar_wave, spec_model, '.', ms=5, alpha=0.7, lw=0.7, color='C1')

plt.step(wave, spec_conv, where='mid', alpha=0.7, lw=0.7, color='C3', label='conv')
plt.plot(wave, spec_conv, '.', ms=5, alpha=0.7, lw=0.7, color='C3')

plt.legend()
plt.show()



# %%
# Sanity check of the convolve and rebin function
step_h = 0.11
step_l = 0.56
w_high = np.arange(-50, 50, step_h)
w_low = np.arange(-30, 30, step_l)
s1 = 2.0
lsf = 3.5
s2 = np.sqrt(s1**2 + lsf**2)
g = np.exp(-0.5*(w_high/s1)**2) / np.sqrt(2 * np.pi) / s1
g_conv_true = np.exp(-0.5*(w_low/s2)**2) / np.sqrt(2 * np.pi) / s2
g_conv_true_high = np.exp(-0.5*(w_high/s2)**2) / np.sqrt(2 * np.pi) / s2
print(w_high.shape, g.shape, w_low.shape)
g_conv = convolve_and_rebin(w_high, g, w_low, np.full_like(w_low, lsf*2.355), pad_pix=5)
plt.close()
plt.step(w_high, g, where='mid', label='orig')
plt.step(w_low, g_conv_true, where='mid', label='true conv')
plt.step(w_high, g_conv_true_high, where='mid', label='true conv high')
plt.step(w_low, g_conv, where='mid', label='conv')

plt.legend()
plt.show()

np.sum(g_conv_true*step_l), np.sum(g_conv_true_high*step_h), np.sum(g_conv*step_l)

# %% [markdown]
# # Telluric fitting


# %%
###############################################################################
###############################################################################

# %%
import lmfit
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

def vac2air(vac):
    """
    Convert vacuum wavelengths [Angstrom] to air wavelengths [Angstrom].
    Ciddor (1996) formula.
    https://www.sdss4.org/dr16/spectro/spectro_basics/#Conversionbetweenvacuumandairwavelengths
    """
    vac = np.asarray(vac, float)
    return vac / ( 1.0 +  5.792105e-2/(238.0185 - (1e4/vac)**2) + 1.67917e-3/( 57.362 - (1e4/vac)**2))


def prepare_palace_data(palace_dir, wrange=[3500, 10000]):
    t_palace = Table.read(f"{palace_dir}/src/palace/data/palace_cont.fits")
    palace_wave = np.asarray(t_palace['lam'] * 1e4)
    palace_wave = vac2air(palace_wave)

    palace_mask = (palace_wave >= wrange[0]) & (palace_wave <= wrange[1])
    palace_wave = palace_wave[palace_mask]
    palace_trans = np.asarray(t_palace['trans'][palace_mask])
    palace_fH2O = np.asarray(t_palace['fH2O'][palace_mask])

    return palace_wave, palace_trans, palace_fH2O


def calculate_mpoly(spec, espec, model, mask, mdegree, sigma_threshold=3.0, niter=0):
    xp = np.linspace(-1.0, 1.0, spec.size)
    xpoly = np.sin(xp * np.pi/2)
    ratio = spec / model
    good = np.isfinite(ratio) & mask
    weights = 1.0 / espec[good]

    def fit_and_get_outliers(w):
        coef = np.polynomial.legendre.legfit(xpoly[good], ratio[good], mdegree, w=w)
        fitted = np.polynomial.legendre.legval(xpoly[good], coef)
        resid = ratio[good] - fitted
        w_resid = resid * w
        mad = np.median(np.abs(w_resid - np.median(w_resid)))
        return coef, np.abs(w_resid - np.median(w_resid)) > sigma_threshold * 1.4826 * mad
    
    # Apply sigma clipping if requested
    if niter > 0:
        for _ in range(niter):
            coef, outliers = fit_and_get_outliers(weights)
            if not np.any(outliers): break
            weights[outliers] = 0.0
    
    final_coef, _ = fit_and_get_outliers(weights)
    return np.polynomial.legendre.legval(xpoly, final_coef)


def residual_3spectra(p, arr_wave, arr_spec, arr_espec, arr_lsf, arr_mask, arr_sky,
                      data=None, mode="residual"):

    # getting stellar template
    stellar_parameters = np.array([[np.log10(p['teff']), p['logg'], p['metal'], p['alpha']]])
    stellar_model = interp_rbfsimplex(data['stellar_templates'], data['stellar_tri'], stellar_parameters)[:, 0]
    stellar_wave = data['stellar_wave']

    # getting transmission curve (See PALACE paper)
    cosz = np.cos(np.deg2rad(p['za']))
    XZ = 1.0 / (cosz + 0.025 * np.exp(-11.0 * cosz))
    rpwv = p['pwv'] / 2.5 - 1
    transmission = np.power(data['palace_trans'], (1.0 + rpwv * data['palace_fH2O']) * XZ)
    transmission_wave = data['palace_wave']

    # iterate over channels
    channels_to_fit = data['channels_to_fit']
    n_channels = len(channels_to_fit)

    (output_resid, output_resid_masked, output_model,
     output_trans, output_mpoly, output_lsf) = ([], [], [], [], [], [])

    for ich, suffix in enumerate(channels_to_fit):
        lsf_corr2 = p[f'lsf_{suffix}']
        wave = arr_wave[ich]
        spec = arr_spec[ich]
        espec = arr_espec[ich]
        lsf = arr_lsf[ich]
        mask = arr_mask[ich]
        sky = arr_sky[ich]

        # lsf_corrected = np.sqrt( np.clip(lsf**2 + lsf_corr2, 0.1, None) )
        lsf_corrected = np.full_like(lsf, lsf_corr2)

        wave_restframe = wave / (1 + p['vel'] / 299792.45)
        stellar_model_convolved = convolve_and_rebin(stellar_wave, stellar_model, wave_restframe, lsf_corrected)
        transmission_rebined = convolve_and_rebin(transmission_wave, transmission, wave, lsf_corrected)

        stellar_model_convolved *= transmission_rebined
        mpoly = calculate_mpoly(spec, espec, stellar_model_convolved, mask, data['mdegree'],
                                sigma_threshold=3.0, niter=3)
        model = stellar_model_convolved * mpoly

        output_resid.append(spec - model)
        output_resid_masked.append((spec - model)[mask] / espec[mask])
        output_model.append(model)
        output_trans.append(transmission_rebined)
        output_mpoly.append(mpoly)
        output_lsf.append(lsf_corrected)

    if mode == "residual":
        return np.concatenate(output_resid_masked)
    elif mode == "model":
        return output_model, output_trans, output_mpoly, output_lsf


def iter_cb(p, iter, resid, *args, **kws):
    if iter % 1 == 0:
        rc2 = np.sum(resid**2) / (len(resid) - len(p))
        print(f"{iter:5d} rchi2={rc2:.3f} Teff={p['teff'].value:.2f} "
              f"logg={p['logg'].value:.3f} [M/H]={p['metal'].value:.3f} "
              f"[a/Fe]={p['alpha'].value:.3f} v={p['vel'].value:.2f} "
              f"PWV={p['pwv'].value:.3f} ZA={p['za'].value:.2f} "
              f"lsf={p['lsf_b'].value:.2f}, {p['lsf_r'].value:.2f}, {p['lsf_z'].value:.2f}")


def fit_star_3spectra(star_data, channels_to_fit=['b', 'r', 'z'], method="leastsq",
                      mdegree=28, verbose=False):
    t_wall0 = perf_counter()
    t_cpu0 = process_time()
    p = lmfit.Parameters()
    p.add("vel", value=0, min=-500, max=500, vary=True)
    p.add("teff", value=6666, min=5000, max=8500)
    p.add("logg", value=4.34, min=1.0, max=5.0)
    p.add("metal", value=0.12, min=-2.0, max=0.5)
    p.add("alpha", value=0.15, min=-0.25, max=0.5)
    p.add("pwv", value=1.0, min=0.1, max=60.0)
    p.add("za", value=30.0, min=0.0, max=80.0)
    p.add("lsf_b", value=1.7, min=0.3, max=3.0, vary=True)
    p.add("lsf_r", value=1.7, min=0.3, max=3.0, vary=True)
    p.add("lsf_z", value=1.7, min=0.3, max=3.0, vary=True)

    palace_wave, palace_trans, palace_fH2O = prepare_palace_data("/Users/ik52/progs/palace/PALACE")

    mask_bright_sky_lines = spec_sky > np.median(spec_sky) + np.std(spec_sky)

    extra_data = dict(
        stellar_wave=stellar_wave,
        stellar_templates=stellar_templates / 1.2e4,
        stellar_tri=stellar_tri,
        mdegree=mdegree,
        palace_wave=palace_wave,
        palace_trans=palace_trans,
        palace_fH2O=palace_fH2O,
        channels_to_fit=channels_to_fit,
    )
    arr_wave = [star_data[suffix]['wave'] for suffix in channels_to_fit]
    arr_spec = [star_data[suffix]['spec'] for suffix in channels_to_fit]
    arr_espec = [star_data[suffix]['espec'] for suffix in channels_to_fit]
    arr_lsf = [star_data[suffix]['lsf'] for suffix in channels_to_fit]
    arr_sky = [star_data[suffix]['sky'] for suffix in channels_to_fit]
    arr_mask = [star_data[suffix]['mask'] for suffix in channels_to_fit]
    arr_mask = [m == 0 for m in arr_mask]
    arr_mask[0] &= (arr_wave[0] > 3600) & (arr_wave[0] < 5800)
    arr_mask[1] &= (arr_wave[1] > 5780) & (arr_wave[1] < 7570)
    arr_mask[2] &= (arr_wave[2] > 7520) & (arr_wave[2] < 9800)


    args = (arr_wave, arr_spec, arr_espec, arr_lsf, arr_mask, arr_sky)
    out = lmfit.minimize(residual_3spectra, p, args=args, kws={'data': extra_data}, nan_policy='omit',
                         method=method, iter_cb=iter_cb if verbose else None)
    print(lmfit.fit_report(out))

    model, trans, mpoly, lsf = residual_3spectra(out.params, *args, data=extra_data, mode='model')

    # Timing
    runtime_wall_s = perf_counter() - t_wall0
    runtime_cpu_s = process_time() - t_cpu0
    extra_data['runtime_wall_s'] = runtime_wall_s
    extra_data['runtime_cpu_s'] = runtime_cpu_s
    extra_data['runtime_s'] = runtime_wall_s
    # Record optimizer details
    extra_data['method'] = method
    try:
        extra_data['nfev'] = int(getattr(out, 'nfev', None))
    except Exception:
        extra_data['nfev'] = None

    return out, model, trans, mpoly, lsf, args, extra_data


def masked_regions(wave, mask):
    """
    Return list of [start, end] wavelength intervals where mask == False.
    """
    bad = ~mask
    edges = np.diff(bad.astype(int), prepend=0, append=0)
    starts, ends = np.where(edges == 1)[0], np.where(edges == -1)[0]
    return list(zip(wave[starts], wave[ends-1]))


def get_zenith_angle(info):
    # Las Campanas Observatory location
    las_campanas_location = EarthLocation(lat=-29.0089*u.deg, lon=-70.6920*u.deg, height=2281*u.m)

    # Midpoint time
    tmid = Time([info['t0'], info['t1']], format='isot', scale='utc').mean()

    # Target coordinates
    target = SkyCoord(ra=info['ra']*u.deg, dec=info['de']*u.deg, frame='icrs')

    # Transform to AltAz and compute zenith angle
    altaz = target.transform_to(AltAz(obstime=tmid, location=las_campanas_location))
    return (90*u.deg - altaz.alt).to(u.deg).value


def plot_results(star_data, out, model, trans, mpoly, lsf, arr_wave, arr_spec, arr_espec, arr_lsf, arr_mask, arr_sky, data=None):
    channels_to_fit = data['channels_to_fit']
    n_channels = len(channels_to_fit)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    colors = ['k', 'purple', 'darkblue']
    for ich, suffix in enumerate(channels_to_fit):
        for mask_region in masked_regions(arr_wave[ich], arr_mask[ich]):
            ax1.axvspan(mask_region[0], mask_region[1], color=colors[ich], alpha=0.1)
        ax1.plot(arr_wave[ich], arr_spec[ich], color=colors[ich])
        ax1.plot(arr_wave[ich], arr_sky[ich], color='pink', lw=0.5)
        ax1.plot(arr_wave[ich], model[ich], color='orange')

        ax1.plot(arr_wave[ich], model[ich] / trans[ich], color='red')
        ax1.plot(arr_wave[ich], mpoly[ich] / np.nanmedian(mpoly[ich]) * np.nanmedian(arr_spec[ich]),
                 color='grey', lw=0.5)
        ax1.axhline(0, ls='--', lw=0.5)

        arr_wave = [star_data[suffix]['wave'] for suffix in channels_to_fit]
        ax1.plot(arr_wave[ich], arr_spec[ich] - model[ich], color=colors[ich], lw=0.5)
        ax1.plot(arr_wave[ich], arr_espec[ich], color=colors[ich], lw=0.3)
        ax1.plot(arr_wave[ich], -arr_espec[ich], color=colors[ich], lw=0.3)
        ax2.plot(arr_wave[ich], trans[ich], color=colors[ich], alpha=0.7, lw=0.5, label='Transmission' if ich == 0 else None)
        ax2.plot(arr_wave[ich], mpoly[ich] / np.nanmax(mpoly), ls=':', color=colors[ich], label='P-Cont. normalizes' if ich == 0 else None)
        ax2.plot(arr_wave[ich], arr_lsf[ich], color=colors[ich])
        ax2.plot(arr_wave[ich], lsf[ich], '--', color=colors[ich], label='LSF corrected' if ich == 0 else None)
        for ax in ax1, ax2:
            ax.set_ylabel("counts / s")

    obs_za = get_zenith_angle(star_data['b']['info'])
    runtime_total = data.get('runtime_s', None)
    method_used = data.get('method', None)

    msg = (f"V={out.params['vel'].value:.2f} km/s "
           f"Teff={out.params['teff'].value:.2f} K "
           f"logg={out.params['logg'].value:.3f} "
           f"[M/H]={out.params['metal'].value:.3f} "
           f"[a/Fe]={out.params['alpha'].value:.3f} "
           f"PWV={out.params['pwv'].value:.3f} mm "
           f"ZA={out.params['za'].value:.1f}°, true={obs_za:.1f}° "
           f"chi2_red={out.redchi:.3f} "
           f"Nfev={out.nfev} "
           f"{method_used} "
           f"Runtime={runtime_total:.2f}s"
           )

    # ax1.text(0.01, 0.99, msg, transform=ax1.transAxes, fontsize=10, ha='left', va='top')
    fig.suptitle(f"{star_data['b']['info']['fib']} GAIA {star_data['b']['info']['id']} "
                 f"{star_data['b']['info']['ra']} {star_data['b']['info']['de']}")
    ax1.set_title(msg, fontsize=10)

    _, _, rms = sigma_clipped_stats(np.concatenate(arr_spec) - np.concatenate(model))
    print(rms)
    ax1.set_ylim(-5*rms, np.nanpercentile(arr_spec, 99)+5*rms)

    plt.legend()
    plt.tight_layout()
    plt.show()



def plot_results_plotly(star_data, out, model, trans, mpoly, lsf, arr_wave,
                        arr_spec, arr_espec, arr_lsf, arr_mask, arr_sky,
                        data=None, output_html_path=None, odir="./"):
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    channels_to_fit = data['channels_to_fit']
    n_channels = len(channels_to_fit)

    hdr = stack_stars[7]['b']['hdr']

    fig = make_subplots(rows=2, cols=1, shared_xaxes=False, vertical_spacing=0.08)

    color_map = ['black', 'purple', 'darkblue']

    # Upper panel: spectra and models
    for ich, suffix in enumerate(channels_to_fit):
        # Accept either per-channel arrays (expected) or a scalar; broadcast scalars
        try:
            lsf_corr_i = lsf[ich]
        except Exception:
            lsf_corr_i = lsf
        if np.isscalar(lsf_corr_i) or (np.asarray(lsf_corr_i).ndim == 0):
            lsf_corr_i = np.full_like(arr_wave[ich], float(lsf_corr_i))

        for x0, x1 in masked_regions(arr_wave[ich], arr_mask[ich]):
            fig.add_vrect(x0=float(x0), x1=float(x1), fillcolor=color_map[ich], opacity=0.1, line_width=0, row=1, col=1)

        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=arr_spec[ich], mode='lines', line=dict(color=color_map[ich]), 
                      name=f'{suffix}: Obs. spectrum', showlegend=(ich == 0)), row=1, col=1)

        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=arr_sky[ich], mode='lines', line=dict(color='pink', width=1), 
                      name='Sky spectrum', showlegend=(ich == 0)), row=1, col=1)

        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=model[ich], mode='lines', line=dict(color='orange'), 
                      name='Model', showlegend=(ich == 0)), row=1, col=1)

        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=(model[ich] / np.where(trans[ich] == 0, np.nan, trans[ich])), mode='lines', 
                      line=dict(color='red'), name='Model (no transmission)', showlegend=(ich == 0)), row=1, col=1)

        scale = np.nanmedian(arr_spec[ich])
        denom = np.nanmedian(mpoly[ich]) if np.isfinite(np.nanmedian(mpoly[ich])) and np.nanmedian(mpoly[ich]) != 0 else 1.0
        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=(mpoly[ich] / denom * scale), mode='lines', 
                      line=dict(color='grey', width=1, dash='dot'), name='Normalization continuum', 
                      showlegend=(ich == 0)), row=1, col=1)

        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=(arr_spec[ich] - model[ich]), mode='lines', 
                      line=dict(color=color_map[ich], width=1), name=f'{suffix}: residual',
                      showlegend=False), row=1, col=1)

        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=arr_espec[ich], mode='lines', 
                      line=dict(color=color_map[ich], width=0.5), name=f'{suffix}: error',
                      showlegend=False), row=1, col=1)
        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=-arr_espec[ich], mode='lines', 
                      line=dict(color=color_map[ich], width=0.5), name=f'{suffix}: error',
                      showlegend=False), row=1, col=1)

        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=trans[ich], mode='lines', 
                      line=dict(color=color_map[ich], width=1), name='Bestfit transmission', 
                      showlegend=(ich == 0), legend="legend2"), row=2, col=1)

        norm_mpoly = np.nanmedian(mpoly[ich])
        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=(mpoly[ich] / norm_mpoly), mode='lines', 
                      line=dict(color=color_map[ich], width=1, dash='dot'), name='Normalization continuum', 
                      showlegend=(ich == 0), legend="legend2"), row=2, col=1
        )

        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=arr_lsf[ich], mode='lines', 
                      line=dict(color=color_map[ich], width=1, dash='longdashdot'), name='DRP LSF', 
                      showlegend=(ich == 0), legend="legend2"), row=2, col=1
        )

        fig.add_trace(
            go.Scatter(x=arr_wave[ich], y=lsf_corr_i, mode='lines', 
                      line=dict(color=color_map[ich], width=1, dash='dash'), name='Bestfit LSF', 
                      showlegend=(ich == 0), legend="legend2"), row=2, col=1
        )
    fig.update_layout(
        title=None,
        # First legend (for upper plot)
        legend=dict(
            orientation='h', yanchor='bottom', font=dict(size=11), xanchor='center',
            y=0.99, x=0.5, bgcolor='rgba(255, 255, 255, 0)',
        ),
        # Second legend (for lower plot)
        legend2=dict(
            orientation='h', yanchor='top', font=dict(size=11), xanchor='center',
            y=0.495, x=0.5, bgcolor='rgba(255, 255, 255, 0)',
        )
    )

    # Y-limits for upper panel based on RMS and 99th percentile
    _, _, rms = sigma_clipped_stats(np.concatenate(arr_spec) - np.concatenate(model))
    y_top_min = -5 * rms
    y_top_max = np.nanpercentile(np.concatenate(arr_spec), 99) + 5 * rms
    if np.isfinite(y_top_min) and np.isfinite(y_top_max):
        fig.update_yaxes(range=[y_top_min, y_top_max], row=1, col=1)

    fig.update_yaxes(title_text="counts / s", row=1, col=1)
    fig.update_xaxes(title_text="Wavelengths (A)", row=2, col=1)

    fig.update_yaxes(range=[0, 2.5], row=2, col=1)

    fig.update_layout(height=800, width=1200, margin=dict(l=60, r=20, t=50, b=100))

    fig.add_annotation(
        text=f"EXPNUM {hdr['EXPOSURE']} {star_data['b']['info']['fib']} GAIA {star_data['b']['info']['id']} {star_data['b']['info']['ra']} {star_data['b']['info']['de']}",
        xref="paper", yref="paper", showarrow=False, font=dict(size=13, weight=600),
        x=0.5, y=1.08, align="center"
    )
    obs_za = get_zenith_angle(star_data['b']['info'])
    runtime_total = data.get('runtime_s', None)
    method_used = data.get('method', None)
    msg = (f"V={out.params['vel'].value:.2f} km/s "
        f"Teff={out.params['teff'].value:.2f} K "
        f"logg={out.params['logg'].value:.3f} "
        f"[M/H]={out.params['metal'].value:.3f} "
        f"[a/Fe]={out.params['alpha'].value:.3f} "
        f"PWV={out.params['pwv'].value:.3f} mm "
        f"ZA={out.params['za'].value:.1f}°, true={obs_za:.1f}° "
        f"&#967;<sup>2</sup>={out.redchi:.3f} "
        f"Nfev={out.nfev} "
        f"{method_used} "
        f"T={runtime_total:.2f}s"
        )
    fig.add_annotation(
        text=msg,
        xref="paper", yref="paper", showarrow=False, font=dict(size=12),
        x=0.5, y=1.05, align="center"
    )

    # Save to HTML
    if output_html_path is None:
        output_html_path = f"{odir}/telluric_fit_{hdr['EXPOSURE']:08g}_{star_data['b']['info']['fib']}.html"
    pio.write_html(fig, file=output_html_path, auto_open=False, include_plotlyjs='cdn')

    return output_html_path

# %%
# run one spectrum
istar = 7
out, model, trans, mpoly, lsf, args, extra_data = \
    fit_star_3spectra(stack_stars[istar], method="least_squares", mdegree=38, verbose=True)


# %%
# plot one spectrum
plot_results_plotly(stack_stars[istar], out, model, trans, mpoly, lsf, *args, data=extra_data, odir="./figs")

# %%
# run all spectra
for istar in range(N_STD):
    out, model, trans, mpoly, lsf, args, extra_data = \
        fit_star_3spectra(stack_stars[istar], method="powell", mdegree=28, verbose=True)
    plot_results_plotly(stack_stars[istar], out, model, trans, mpoly, lsf, *args, data=extra_data, odir="./figs")
