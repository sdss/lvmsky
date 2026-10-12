# LSF references for Sky Decomposition

## Sample and individual LSF recovery

The line-spread function (LSF) is the response to a monochromatic input at wavelength $\lambda$. Displacement $\Delta\lambda$ is output minus input wavelength. B/R/Z are the blue, red and near-infrared spectral channels. Sci denotes the science telescope; SkyE and SkyW denote the east and west sky telescopes. SP1–SP3 identify the three spectrographs. MJD means Modified Julian Date. The analyzed sample contains 558 local DRP 1.3.2 exposures selected with Moon elevation >60°, lunar illumination >0.95 and Sun elevation <−18°; it is the analyzed bright-Moon sample, not every survey observation.

SolarBlue fits each physical fiber independently to the observed FLUX+SKY spectrum, with its masks and variance. Continuum amplitudes, emission-line amplitudes and the LSF are refined alternately. The B channel combines sky emission with solar Fraunhofer absorption carried by the Moon and zodiacal-light continua; R/Z use sky emission lines. The high-resolution solar source and emission-line model are integrated over the native detector bins. In B, the LSF-dependent model is

$$
y_B=A_{\rm line}\theta+C(\lambda)(G\theta)+m_{\rm diffuse}+m_{\rm background},
$$

where $G$ is the precomputed solar convolution tensor and $C$ the fitted smooth solar-continuum envelope. For fixed spectral amplitudes the LSF step is a constrained quadratic least-squares problem. The displacement second-difference penalty uses fraction 0.01 of the mean data-Hessian diagonal; wavelength roughness and fallback-prior fractions are zero. Spectral parameters are refitted after updating the LSF. Refinement stops when the maximum pointwise relative changes in W50/W90 and total chi-square are <1%, and centroid changes are <0.015 Å, or after 30 cycles. Convergence flags and original source settings are retained in the FITS provenance.

Each channel uses the continuous density

$$
L_i(\Delta\lambda,\lambda)=\sum_{k=0}^{10}\sum_{b=0}^{3}\theta_{i,kb}M_k(\Delta\lambda)B_b(\lambda),\qquad
\theta_{i,kb}\geq0,\quad\sum_k\theta_{i,kb}=1.
$$

The 11 cubic M-splines have unit integral and displacement support ±2.5 Å; four cubic wavelength B-splines form a partition of unity on their stored interval. Coefficients are dimensionless; density is in Å⁻¹. Each coefficient column increases toward component $k=5$ and decreases thereafter. Exact wavelength knots, bounds and degrees are saved. The displacement knots are reconstructed as $t_j=-2.5+5j/14$ Å for $j=0,\ldots,14$. Each local cubic B-spline $\beta_k$ on $t_k,\ldots,t_{k+4}$ is zero outside its support and gives $M_k=4\beta_k/(t_{k+4}-t_k)$. This recipe specifies the continuous density independently of the installed model code. No centroid-zero or fixed-width condition is imposed.

## Accepted measurements

The verified aperture-source Gaia G-band flux divided by synthetic G-band flux of the current observed FLUX+SKY spectrum must be <0.1. The separate <0.01 flag is retained; unknown contamination is not accepted. Gaia is taken only from `$SAS_BASE_DIR/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2/gaia/fibers/`. Table checksums and current SFrame size/timestamp are checked before photometry-cache reuse. Synthetic G uses current FLUX+SKY with verified passband coverage; coordinates and passband metadata are checked by the parent analysis.

Technical acceptance requires an available, converged channel, positive valid-pixel count, main status `Solved`, `fibstatus=0`, fitted valid fraction ≥0.8, finite chi-square and no recorded LSF failure. Supplemental fits retain their original <0.01 input restriction. Accepted channels total 2,556,278 across 1,871 physical fibers; Gaia-table availability and completed fitting are distinct counts.

## Epochs and lunations

Camera calibration IDs are read from each SFrame's calibration-file provenance. They define the actual camera epoch; nominal epoch definitions describe the parent timeline. The portable selector trusts the actual camera IDs in target provenance rather than reassigning epochs from a date-only timeline. Each exposure timestamp $t_i$ is assigned to the nearest computed full-Moon timestamp $T_j$:

$$
\ell(i)=\operatorname*{arg\,min}_j|t_i-T_j|.
$$

The stored full-Moon UTC timestamp labels that observed lunation. A lunation crossing a camera-epoch boundary is split. This product uses the parent report's epoch–lunation groups, with no additional illumination or waxing/waning subdivisions.

| Scope | Contributing population |
|---|---|
| Survey | All accepted measurements of the same physical fiber and channel; descriptive pooling across epochs |
| Epoch | Same fiber and channel within one actual camera epoch |
| Epoch–lunation | Same fiber, channel, actual epoch and full-Moon label |

{{SUMMARY_TABLE}}

For a target observation, restrict candidate references to its physical fiber, channel and actual camera epoch. Choose the nearest available epoch–lunation reference by its recorded weighted observation MJD; the full-Moon label identifies the group, while first/last MJD record actual sampling. The supplied loader performs this selection; no reference is silently borrowed from another fiber or epoch. Target times can lie outside the sampled first/last MJD while retaining the same camera calibration ID. This nearest-neighbor transfer assumes stability within that epoch; its predictive accuracy away from sampled nights has not been independently established. The survey-wide file is descriptive rather than the default calibration.

## Averaging

Within each group, verify identical bases and flatten each 11×4 matrix as $x_{ij}$. Repair archival float32 unit-sum rounding at float64 precision. Compute

$$
c_j=\operatorname{median}_i x_{ij},\qquad
s_j=\max\left(1.4826\operatorname{median}_i|x_{ij}-c_j|,10^{-6}\right),\qquad
d_i=\sqrt{\frac1{44}\sum_j\left(\frac{x_{ij}-c_j}{s_j}\right)^2}.
$$

One scalar weight applies to the whole surface:

$$
r_i=\min(1,2.5/d_i),\qquad q_i=1/n_{h(i)},\qquad
\alpha_i=\frac{q_ir_i}{\sum_\ell q_\ell r_\ell},\qquad
\bar\theta_{kb}=\sum_i\alpha_i\theta_{i,kb}.
$$

Here $n_{h(i)}$ counts exposures in the same observed night. Set $r_i=1$ when $d_i=0$ or the group has fewer than three exposures. Nights have equal mass before robust downweighting. Survey, epoch and epoch–lunation averages use the same formula independently on their respective populations. These are empirical robustness weights, not inverse formal variances. Positive weights summing to one preserve nonnegativity, normalization and the common shape inequalities. No centroid alignment is applied, so centroid variation can broaden a pooled profile.

## Figures and statistics

All figures are 2D. Maps of FWHM, W50, W90, centroid and asymmetry use wavelength on X and one available physical fiber per ordered row on Y. Fiber IDs are row labels: spacing is categorical, missing fibers are not filled, and no interpolation between fibers is used. Maps evaluate the whole-sample mean coefficient solutions, separately for every telescope/spectrograph group and channel.

For profile examples, choose the working fibers nearest the lower, middle and upper blue-channel slit positions in each group. This geometry-based choice does not favor stable fits. Their density maps use the same survey solution. Every wavelength-curve figure has a legend: **black = whole-sample mean LSF; colored = epoch–lunation mean LSF of the same fiber**. The colorbar gives that channel's weighted observation MJD.

FWHM is full width at half maximum. $Q_p$ is the pth percentile of the normalized profile; W50 = $Q_{75}-Q_{25}$, W90 = $Q_{95}-Q_{05}$, and asymmetry = $(Q_{95}+Q_{05}-2Q_{50})/\mathrm{W90}$. Centroid is $\int\Delta\lambda\,L(\Delta\lambda,\lambda)\,d\Delta\lambda$. Quantiles and widths use the established 1201-point displacement evaluator. The plots cover the stored wavelength-basis intervals. Curve differences include fit noise and recovered instrumental differences, not a calibrated prediction error.

`COEF_SCATTER` contains coefficient standard deviations describing weighted exposure dispersion. `COEF_CONDERR` contains coefficient standard errors estimating conditional mean error using fixed weights and observed-night groups; it assumes independent nights and does not propagate uncertainty in the weights. One exposure gives undefined scatter; one night gives undefined conditional error, recorded as NaN. Neither field establishes predictive transfer accuracy. Writing $v_i=\theta_i-\bar\theta$, their saved diagonals are

$$
\mathrm{COEF\_SCATTER}_{kb}^{2}=\frac{\sum_i\alpha_i v_{i,kb}^{2}}{1-\sum_i\alpha_i^{2}},\qquad
\mathrm{COEF\_CONDERR}_{kb}^{2}=\frac{H}{H-1}\sum_{h=1}^{H}\left(\sum_{i\in h}\alpha_i v_{i,kb}\right)^2,
$$

where $H$ is the number of observed nights. Full contributor weights are retained in the separate calculation archive.

Survey-profile characterization at the existing channel anchors is tabulated as median [16th, 84th percentile] across the 1,871 physical fibers. These intervals describe variation between fiber mean profiles, not uncertainty on an individual profile.

{{CHANNEL_METRICS}}

## Files and use

The publication directory contains this README, the HTML report, `figures/`, the local formula-rendering asset, the portable `load_lsf.py`, and compact coefficient files under `survey/`, `epoch/` and `lunation/`. Full calculation records are kept separately in the sibling `lsf_reference_audit/`; they are unnecessary for using the LSF.

| File location | Content |
|---|---|
| `survey/lsf_survey.fits` | All analyzed epochs pooled; 5,613 fiber/channel surfaces; about 7.4 MiB |
| `epoch/lsf_<epoch>.fits` | One actual camera epoch per file; 18 files |
| `lunation/lsf_<epoch>_<full-Moon-UTC>.fits` | One epoch–lunation group per file; 20 files |
| Adjacent `.README.md` | Short calculation identification and a link to this common method |

Every FITS has the same schema:

| Extension | Content |
|---|---|
| `REFERENCES` | Fiber, telescope, spectrograph, channel, epoch, period, exposure/night counts, dates and statistics |
| `LSF_COEF` | Float64 array `(reference, 11 displacement components, 4 wavelength basis functions)` |
| `LSF_KNOTS`, `LSF_CONFIG`, `WAVE` | Exact wavelength knots, degrees, support, configuration and native grid |
| `COEF_SCATTER`, `COEF_CONDERR` | Dimensionless coefficient standard deviations and conditional standard errors, each `(reference, 11, 4)` |
| `SOURCES`, `SOURCE_CONFIG` | Original input provenance and source fit settings |

`reference_id` is local to the file. Retrieve a row by physical fiber and channel in `REFERENCES`; that row's array position indexes `LSF_COEF` and `LSF_KNOTS`. The first matrix axis is displacement component; the second is wavelength basis. Native-bin convolution uses the existing continuous-spline implementation.

Copy this entire `lsf_calibration/` directory. Python needs NumPy, Astropy, and a compatible installed `skysub.sky_decomp` with continuous-spline SolarBlue and `fixed_lsf_state` support. The included `load_lsf.py` is the portable loader; it reads only local compact FITS. Original SAS paths in `SOURCES` are provenance, not runtime dependencies.

Run from the copied directory, or add it to Python's import path:

```python
from load_lsf import load_lsf_from_sframe
from skysub.sky_decomp.solar_lsf import SkyDecompPalaceCorrSolarBlueLSF

state = load_lsf_from_sframe("lvmSFrame-00044724.fits", fiberid=323)
model = SkyDecompPalaceCorrSolarBlueLSF(wave, **target_model_options)
result = model.fit(flux, ivar, fixed_lsf_state=state)
```

To explicitly use the mean over the entire analyzed survey for that same fiber and all three channels:

```python
state = load_lsf_from_sframe("lvmSFrame-00044724.fits", fiberid=323, scope="survey")
result = model.fit(flux, ivar, fixed_lsf_state=state)
```

`scope="survey"` reads `survey/lsf_survey.fits` and pools the analyzed calibration epochs. It needs the target `SLITMAP` and `WAVE`; timestamps and camera-epoch provenance are not required in this mode. It validates the physical fiber and exact native grid. `scope="lunation"` is the default and selects the nearest measured lunation inside each actual camera epoch. Both modes load saved coefficients and perform zero LSF optimization.

`wave`, `flux`, `ivar` and `target_model_options` are the normal target inputs for Sky Decomposition. In the default lunation mode, the helper reads `OBSTIME`, `SLITMAP` and the three camera calibration filenames from the target SFrame, then chooses the closest available weighted observation MJD within each camera's epoch. B/R/Z can come from different files. Selection works for arbitrary exposure timestamps **within the fiber/channel/epoch coverage of this calibration sample**; an uncovered epoch or fiber raises a clear error. It never substitutes a different epoch or the survey pool. The supplied data cover 18 actual epoch IDs, not every possible calibration epoch of the survey.

With explicit provenance, use `select_reference_lsf(directory, fiberid, {"B": epoch_B, "R": epoch_R, "Z": epoch_Z}, observation_mjd)`. To select an exact file, use `load_reference_lsf(file, fiberid)`; a B/R/Z file mapping is also accepted. `state.metrics[channel]` records the chosen file, reference ID, calibration epoch and exposure/night counts.

Loading copies coefficients exactly and performs zero LSF optimization. The target native wavelength grid and displacement support must match; other spectral parameters are still fitted. The high-level loader selects survey pooling explicitly with `scope="survey"`; direct `load_reference_lsf` calls require `allow_epoch_pooling=True`. Historical native-tap/median-stack inputs require the compatible individual-target SolarBlue path rather than substituting these matrices into the old kernel format. The high-level state loader uses installed model code; the FITS definition above also permits an independent implementation of the density.
