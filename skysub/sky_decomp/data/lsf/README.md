# Quality-screened per-fiber survey LSF

The B/R/Z ASCII tables contain a separate normalized continuous LSF for each of 1871 measured physical Sci/Sky fibers in DRP 1.3.2. Standard-star and unmeasured/hardware-bad positions are not populated by interpolating another fiber.

## Files and representation

- `survey_lsf_B.csv`, `survey_lsf_R.csv`, `survey_lsf_Z.csv`: one physical fiber per row, coefficient matrix `coef_ii_jj`, night coefficient scatter, coverage and repeatability.
- `wavelength_knots.csv`: wavelength B-spline knots in Angstrom, indexed by channel.
- `coverage.csv`: per-fiber/channel coverage, quality flags and metrics.
- `manifest.json`: schema, cubic M-spline displacement basis, wavelength degrees, native grid identity, checksums and input provenance.

For displacement d and wavelength w, L(d,w) = sum_ij coef_ij M_i(d) B_j(w). M_i are unit-integral cubic M-splines with 11 components and uniform displacement knots over [-2.5,+2.5] Angstrom. B_j are the stored wavelength B-splines. The coefficient matrix is nonnegative, each wavelength column sums to one, and central-component inequalities are preserved. CSV floats use 17 significant digits. This is the full wavelength-dependent profile, not just Gaussian FWHM.

## Selection and averaging

The authoritative Gaia source is `$SAS_BASE_DIR/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2/`; `gaia/fibers/` is its resolved internal loose-file layout. Verified Gaia G / current observed FLUX+SKY synthetic G <0.1, with the <0.01 flag retained, and confirmed OBJECT/physical-scale/spectrum guards. Archived channel-solve reduced chi-square is <2 for Blue and <3 for R/Z. Shape screening uses the five channel anchors and both wavelength endpoints: FWHM 0.9-3.0 Angstrom, W90/W50 1.8-3.5, |barycenter| <=0.25 Angstrom, |asymmetry| <=0.25, coefficient maximum <=0.85 and fewer than two nearly flat coefficient columns (range <0.025). For normalized LSF area percentiles Qp, W50 = Q75-Q25, W90 = Q95-Q05, and asymmetry = (Q95+Q05-2*Q50)/W90. The anchor wavelengths are B: 3800,4250,4700,5200,5650; R: 5900,6250,6600,7000,7350; Z: 7600,8000,8500,9000,9600 Angstrom. These are explicit conservative operational screens, not instrument-calibrated confidence limits.

After screening, every retained observed MJD night gets equal weight; exposures within a night get equal weight. Coefficient distances across epochs never downweight a clean night. Each physical fiber and channel is measured separately, with no cross-fiber smoothing, borrowing or centroid alignment. The survey reference represents the retained observed-night mixture, rather than all survey exposure times or a single calibration epoch.

All Blue fibers are measured. Z fibers 84/85 carry `elevated_chi2`: their diagnostic tables use chi-square <5 because no <3 archival solution exists. Strict loading refuses these Z rows; `allow_warning=True` is an explicit diagnostic opt-in. Fixed-Blue loading requests B only and uses its measured per-fiber reference. R/Z initialization then uses the target DRP profiles and is refined on target sky lines.

## Precision and epoch variation

`night_scatter_ii_jj` is sample scatter of nightly coefficient means, not uncertainty of the reference mean. `night_fwhm_scatter` is 1.4826*MAD of nightly FWHM at the channel anchor. `loo_night_width_error_pct` is the median absolute percentage difference between each nightly mean's FWHM and a reference constructed from all other retained nights, across the five channel anchors. Selection remains based on archived fits. This checks recovered-fit repeatability, not held-out spectral residual prediction or absolute hardware accuracy. `historical_width_error_pct` compares the historical full-sample reference directly to the selected nightly means; it is a descriptive departure, not an old-reference leave-one-night-out estimate.

`epoch_width_range_pct` is the largest, across the five anchors, P90-P10 range of equal-night epoch-profile FWHM divided by survey-profile FWHM. `epoch_dependent` marks >3%; it is a descriptive flag. A fixed survey profile cannot reproduce a genuine calibration change such as fiber 973 in every epoch. `n_eff` is 1/sum(exposure_weight**2), not an independent-night count. `n_gaia_clean_1pct` retains the number of selected exposures also passing the separate <0.01 Gaia criterion.

## Loading

```python
from skysub.sky_decomp.lsf_reference import load_survey_ascii_lsf
state = load_survey_ascii_lsf("skysub/sky_decomp/data/lsf", 1812, channels=("B",))
```

The existing `--blue-lsf reference` mode defaults to this bundled directory and fixes each physical fiber's B profile; R/Z are still refined on the target spectrum. Existing FITS calibration directories remain supported through `--lsf-calibration`. The real-spectrum check is reproducible with `python -m analysis.lsf_epoch_coverage.check_fixed_blue_ascii`.

Regenerate with `python -m analysis.lsf_epoch_coverage.build_survey_ascii` in the documented lvmdrp311 environment; publish the figures with `python -m analysis.lsf_epoch_coverage.publish_survey_ascii`.

[Published validation report](http://lumen.abudhabi.nyu.edu/sas/sdsswork/users/u6058164/lvm_lsf_survey_ascii_drp132/survey_lsf_ascii.html).
