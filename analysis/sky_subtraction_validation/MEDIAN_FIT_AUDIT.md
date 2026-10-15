# Median-fit investigation

## Measured fit quality

RMS = unweighted √mean[(observed median − fitted model)²] on valid fitted pixels. Flux units: 10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹. Normalized RMS = √mean[(data−model)² MED_IVAR]. σ = MED_IVAR⁻¹ᐟ² is the approximate input error of a masked median, propagated from CFrame IVAR under independent Gaussian fibers sharing a signal. It measures input noise; intrinsic spatial variation, covariance and sky-model uncertainty are separate.

| Exposure | Telescope | RMS B / R / Z | Normalized RMS B / R / Z |
|---|---|---|---|
| 16819 | Sci | 0.932 / 0.555 / 0.837 | 16.1 / 26 / 43.8 |
| 16819 | SkyE | 1.36 / 0.645 / 0.842 | 8.63 / 11.6 / 20.3 |
| 16819 | SkyW | 0.956 / 0.57 / 0.861 | 6.5 / 10.4 / 17.3 |
| 59750 | Sci | 0.587 / 0.336 / 0.291 | 17.7 / 22.4 / 21.3 |
| 59750 | SkyE | 0.522 / 0.263 / 0.288 | 3.79 / 4.46 / 5.11 |
| 59750 | SkyW | 0.798 / 0.788 / 0.564 | 4.27 / 10.5 / 11.7 |

These fits have substantial model mismatch above input noise. `Solved` identifies numerical optimization status; it does not identify a noise-limited fit or accurate LSF.

## Reproduction and isolation

A repeated Science-median fit for 16819 uses the recorded input medians, masks and initial survey profiles. Its coefficient-rendered spectrum agrees with the solver's own best-fit vector to 2.85×10⁻¹⁴ **in display units**. The repeat differs from the saved model by at most 0.0144 display units; the large displayed residuals persist. Thus reconstruction from compact parameters does not explain the large residuals in this repeated exposure.

Local Z-line checks fit the same modeled emission template plus amplitude, baseline and slope in ±4 Å windows, first at zero shift and then with a free displacement within ±0.5 Å. They use input median IVAR. Positive shift means moving the template to longer wavelengths. Select up to 12 strongest modeled Z peaks separated by ≥30 Å; show the three regions with largest saved-fit RMS. The JSON records all probes and all three telescopes.

For Science 16819, shifts change sign over short wavelength intervals (approximately −0.10 Å near 9376 Å and +0.15 Å near 9476.5 Å), and reduce local mismatch considerably. Near 9376 Å, zero-shift local RMS 3.78 becomes 0.634 after freeing the shift. SkyE/SkyW yield shifts −0.102/−0.098 Å at the same feature; the similar displacement across telescopes is additional evidence of a common spectral mismatch. This establishes a position-sensitive component of the residual. It does not establish whether the cause is wavelength calibration, template wavelengths, line blends, or inadequate LSF freedom. A four-basis smooth LSF surface cannot describe arbitrary short-scale displacements. Broader or less restrictive fits require independent validation before adopting them.

## Final subtraction is a separate measurement

The Science-median decomposition is used for the OH correction and LSF comparison; the final sky coefficients come from the two sky-telescope fits and the trained ML predictor. The report therefore shows the fit RMS separately from the matched-fiber DRP → ML → OH → final subtraction RMS.

The existing whole-exposure line-scale metric gives 16819 DRP → SKYML values B 1.09 → 0.94, R 4.77 → 3.37, Z 7.96 → 5.21 in units of observed per-fiber input σ. On the same metric, 3469 R/Z and 59750 Z are worse than DRP. Lower line-scale RMS is useful evidence about subtraction, not proof of recovered LSF accuracy or a calibrated χ².

For 59750, an additional controlled comparison holds the Science-median LSF fixed and changes only the coefficient vector. High-pass Z RMS is 0.219 with the decomposition coefficients and 0.676 with ML coefficients; the saved OH correction leaves it at 0.676. Provenance records a complete OH rollback: `Median OH correction produces invalid line populations`. The unconstrained linear temperature-tilt correction would make some intrinsic populations negative; the current integration responds by restoring the entire original ML coefficient vector. This protection is physically justified, but the complete rollback also discards the fitted brightness corrections; their benefit requires a separate constrained test. The actual matched-fiber final median has RMS 0.648 versus DRP 0.358. Thus this regression is already present before OH correction, and changing individual-fiber LSF alone does not explain it. The fixed-median-LSF comparison is a diagnostic on the observed median, distinct from a median of individually subtracted fibers.

## Reproduce

Run `python analysis/sky_subtraction_validation/audit_median_fits.py`, then `python analysis/sky_subtraction_validation/report_fibers.py`. `median_fit_audit.json` stores measured median-fit RMS and local line probes. The shared compact renderer supplies model spectra; the audit does not alter FITS or apply local shifts to the final sky.
