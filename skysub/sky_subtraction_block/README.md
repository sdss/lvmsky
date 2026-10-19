# Sky-only block LSF subtraction

This method shares the median decomposition, ML predictor, blue continuum residual transfer and Science-median OH correction with `sky_subtraction`. Its LSF stage uses selected SkyE and SkyW fibers only, separately for each spectrograph in the R and Z channels, where sky lines constrain the profile. B is not fitted: every Science fiber keeps the survey LSF of the preparation there (it still enters the model, the flux and the residuals). Science spectra never enter this LSF fit.

1. Select hardware-good, Gaia-clean Sky fibers (Gaia / observed G < 0.01). Prepare masked telescope medians, decompose them, predict the Science component coefficients, and determine the B continuum and OH corrections.
2. Decompose each selected Sky fiber with its telescope-median LSF and O2 temperature fixed. Each fiber gets its own nonnegative component coefficients. Freeze these coefficients for the block fit.
3. Solve one native-pixel weighted least-squares QP for each of six spectrograph × (R, Z) blocks. The LSF density is `sum H[k,q,s] M[k](offset) B[q](wavelength) S[s](fiberid)`. The initial slit model uses four cubic B-spline controls over the full spectrograph fiber-ID range. The wavelength basis is the existing four-control cubic basis. SkyE and SkyW use their own component coefficients and telescope transmission arrays.
4. Enforce `H >= 0` and `sum_k H[k,q,s] = 1` for every wavelength/slit control pair. The centroid is unconstrained by default; `--centered` (`centered=True`) also enforces `sum_k mu[k] H[k,q,s] = 0`, with `mu[k]` the exact first moment of the normalized cubic M-spline. Partition of unity propagates the equalities to every wavelength and slit position. There is no wavelength-solution shift or added shape/width constraint. Native IVAR supplies weights; input/protected-line masks apply. Five Angstroms at channel edges are trimmed from the fit because the affine designs hold neighboring-channel leakage fixed.
5. Evaluate the fitted R and Z surfaces at each Science fiber, keep that fiber's survey B LSF, and render the common ML coefficients with this LSF. Apply the same B continuum transfer and Science-median OH correction to every fiber, through the shared renderer. Store the actual final SKY and FLUX = observed − SKY.

```sh
python -m skysub.sky_subtraction_block /path/lvmCFrame-00059750.fits \
  --ensemble /path/trusted_ensemble.pt --gaia-dir /path/gaia/fibers \
  --output /separate/products/lvmSFrame-skyml-block-00059750.fits --jobs 8
```

For a controlled LSF comparison, add `--seed /path/lvmSFrame-skyml-00059750.fits` to reuse that exact exposure's saved median/ML/correction preparation. The input SHA256 and row/wavelength identities must match. The seed's individual Science LSF fits are replaced completely. A run without a seed uses the existing survey-only preparation pipeline and stores its compact preparation beside the block product. Existing destinations are protected; working filenames have no version suffixes.

## Compact storage and diagnostics

The shared [compact product contract](../sky_subtraction/FITS_STRUCTURE.md) applies. `BLOCK_LSF_COEF` stores the six blocks × offset × wavelength × slit controls; `BLOCK_SLIT_KNOTS`, `BLOCK_WAVE_KNOTS` and `BLOCK_LSF_META` store exact domains, degrees, solver status, rank, input counts, chi-square, and moment checks. `SKY_FIBER_COEF` records every selected Sky fiber's identity, component vector and decomposition status. The usual `LSF_COEF/KNOTS/META` contain evaluations of the R/Z global solution at Science fibers together with their survey B LSF, so existing `SkyMLReconstructor.sky()` works directly. No intermediate model images are persisted. `SKYML_META` records preparation/input hashes, settings and implementation hashes. `SkyBlockReconstructor.global_state()` and `.sky_fiber()` reconstruct block profiles and selected Sky-fiber fits in memory.

Input IVAR describes observed noise only. Fixed decomposition coefficients and common calibration errors induce correlations not propagated by this QP. Strictly centered LSFs can broaden or distort when individual wavelength offsets exist; comparison with measured Science profiles therefore retains their original centroids. Slit intervals outside the selected Sky range are extrapolation relative to the data, even though the B-spline evaluation stays inside the full spectrograph domain. Rank and coverage are recorded; QP success establishes numerical feasibility, not recovery accuracy. Slit flexibility and additional regularization remain empirical choices for subsequent tests.

Validation: `python -m pytest skysub/sky_subtraction_block/tests -q`. The injected test recovers a varying global surface from multiple fibers and checks continuous positivity, normalization and centering between control nodes.
