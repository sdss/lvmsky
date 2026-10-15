# Median LSF metrics from saved coefficients

Example: exposure 59750, MJD 61192, SkyE and SkyW median decompositions. Source: the existing compact SkyML product recorded in `59750_median_lsf_coefficient_metrics.json`. B is fixed to the selected-fiber survey reference; R and Z are fitted median profiles. These are effective median LSFs, not individual-fiber measurements. No new spectral fit is performed.

## Direct coefficient contractions

Write `K(u, lambda) = sum(k,p) H[k,p] M[k](u) B[p](lambda)`, with unit-integral M-splines, normalized coefficients, and `u` in Angstrom. Precompute `v[r,k] = integral u**r M[k](u) du`. Then `m[r](lambda) = B(lambda) @ H.T @ v[r]`.

The coefficient contractions are evaluated before the wavelength basis. Cubic M-splines multiplied by `u**r`, for r <= 3, are integrated exactly up to floating-point precision by four-point Gauss quadrature on every knot interval. Tail vectors use the existing analytic M-spline bin integrals.

- Centroid: `mu = m1`, in Angstrom; linear in H.
- Raw second moment: `m2`, in Angstrom squared; linear in H.
- Central variance: `V = m2 - mu**2`; nonlinear unless mu is fixed, including mu = 0. Width `sigma = sqrt(V)` is nonlinear.
- Tail fraction: `T(a) = integral(abs(u) > a) K(u) du`; linear in H for a fixed origin and threshold. Plotted thresholds 1.5 and 2.0 Angstrom are illustrative, not acceptance criteria. Tails refer to the entire stored profile on support +/-2.5 Angstrom; they are not a separately fitted wing component and do not constrain unmodeled flux beyond this support.
- Raw third moment: `m3`, in Angstrom cubed; linear in H and sensitive to both offset and asymmetry.
- Central third moment: `c3 = m3 - 3*mu*m2 + 2*mu**3`; nonlinear for variable mu.
- Standardized central skewness: `gamma1 = c3 / V**1.5`; dimensionless and nonlinear. Positive values describe a positive third central moment toward larger wavelength offsets.

The plots preserve the original fitted profiles and their centroids. They compare raw and central moments rather than silently recentering the profiles. No coefficient covariance or uncertainty band is available in these products. Smooth curves describe the fitted spline parameterization; they do not independently establish LSF recovery accuracy.

## Relationship to a linear wavelength-correction formulation

With fixed unconvolved component spectra and fixed output sampling, a free effective kernel K is linear in its full coefficients. Its centroid is an exact linear functional of those coefficients, without a spectral Jacobian. Requiring the centroid to follow a correction spline is therefore a linear equality, `m1(lambda) = sum(j) d[j] G[j](lambda)`. A common wavelength basis makes this a coefficient-level equality.

A centered representation can subsequently be defined exactly as `L(v,lambda) = K(v + mu(lambda),lambda)`. The centroid of L is zero. This reparameterization does not imply that independently parameterized centered L and a finite shift have a jointly linear spectral model. In particular, enforcing smoothness of the centered shape across the slit while allowing independent fiber shifts requires a separate mathematical treatment. The effective-kernel moment formulation is the appropriate starting point for that discussion.

## Checks and reproduction

Run from the repository with the existing lvmdrp311 environment:

```bash
OPENBLAS_NUM_THREADS=1 /home/ik52/miniconda3/envs/lvmdrp311/bin/python analysis/sky_multifiber_method_discussion/plot_median_coefficient_metrics.py
```

The generator checks normalized nonnegative coefficients, positive central variance, physical tail fractions, four-point versus eight-point quadrature, and independent dense-profile integration at three wavelengths in every channel/telescope. CSV columns contain the physical values; tail fractions are multiplied by 100 only for plotting.
