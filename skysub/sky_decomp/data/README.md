# Sky decomposition data bundle

This directory is the portable data root for the production sky-decomposition
methods. It contains the physical Moon/Zodi model, the exact PALACE tables used
by the validated fits, the frozen full-grid residual PCA basis, and an unchanged
copy of the Meftah SOLAR-HRS text spectrum used by the historical decomposition.

The historical `SkyDecomp` and non-split `SkyDecompLSFSurfaceIterative`
defaults are not redirected here. Their existing `base_dir` and
`skysub/palace` contracts remain unchanged. The Moon/Zodi method and the
`lsf-surface-iterative-split-zodi` command-line mode use this directory by
default; any class can also use the copied historical layout explicitly with
`base_dir=.../data`.

## Directory layout

```text
data/
├── README.md
├── bundle_manifest.json
├── Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt
├── moon_zodi/
│   ├── README.md
│   ├── model_manifest.json
│   ├── jpl_de432s_short_planetary_ephemeris.bsp
│   ├── meftah_solar_hrs_disk_integrated_v1_1_vacuum_velocity_step_2kms.npz
│   ├── eso_skycalc_rolo_moon_albedo.dat
│   └── eso_skycalc_leinert_zodiacal_light.dat
├── residual_pca/
│   ├── palace_aijc_full_native_pca_v1.npz
│   ├── palace_aijc_line_amplitude_pca_v1.npz
│   ├── palace_aijc_vnf_coefficient_pca_prep_v1.npz
│   └── palace_aijc_vnf_coefficient_line_amplitude_pca_v1.npz
└── palace/PMD/
    ├── palace_line_telluric_r4m_v1.fits
    ├── pmd_popmodel_OH_telluric_upper_parity_lsf_adam_25000_v1.dat
    ├── pmd_popmodel_OH_h_family_default_ef_v1.dat
    ├── pmd_popmodel_OH_joint_v2_updated.dat
    ├── pmd_refcont_canonhyb_v1.dat
    ├── pmd_refcont_joint_native_adam_invsky_p2_10000iter.dat
    ├── pmd_intdata_atom.dat
    ├── pmd_intmodel_Orc.dat
    └── pmd_popmodel_O2.dat
```

`bundle_manifest.json` records the expected PALACE suffixes and SHA-256 digest
of every copied decomposition table. `moon_zodi/model_manifest.json` records
the four physical-asset digests, the frozen model parameters, training domain,
formula, grid, and checkpoint provenance. The new method validates both
contracts before fitting.

## Meftah SOLAR-HRS source

`Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt` is copied byte-for-byte from the
pre-existing `skysub/` source and retains its original filename and header. It
is the disk-integrated SOLAR-HRS v1.1 irradiance spectrum attributed in the file
to Meftah et al. The columns are wavelength in nm and disk-integrated solar
spectral irradiance in W m-2 nm-1. SHA-256:

```text
d141e880a3486c039610bcd8de2b6b08930796de6588857ba643c493fbd687bd
```

The authoritative catalogue record is CDS catalogue VI/159, associated with
Meftah et al. (2023), *Remote Sensing*, 15, 3560. The local source file itself
does not state a redistribution licence. Preserve the attribution and verify
redistribution terms before publishing the raw text file outside the project.

The compact Moon/Zodi NPZ is derived from this source. Its former suffix `dv2`
meant a constant logarithmic-grid velocity step of 2 km/s; it did not mean
version 2. The canonical filename now spells this out. The NPZ contains:

- `wave_vacuum_angstrom`: 157,365 float64 vacuum wavelengths;
- `flux_disk_integrated`: float64 disk-integrated irradiance samples;
- `grid_velocity_step_kms`: scalar value 2.0;
- source metadata and source SHA-256.

The transformation is documented in `moon_zodi/model_manifest.json`. The
observed 12,401-pixel LVM spectrum is never replaced, rebinned, cropped, or
interpolated by this asset preparation; the high-resolution carrier is
projected onto the unchanged native detector-pixel boundaries at prediction
time.

## JPL DE432s ephemeris

`jpl_de432s_short_planetary_ephemeris.bsp` is the unchanged JPL/NAIF
`de432s.bsp` binary SPK kernel. `DE432` identifies a JPL Development Ephemeris;
the trailing `s` identifies the shorter, compact distribution. It provides
Solar-System positions and velocities. This method selects it explicitly
through Astropy to compute exposure-midpoint Sun, Moon, and Earth geometry and
the reflected-light velocity correction. It is not a sky spectrum, trained
parameter file, or Moon/Zodi intensity table.

Original source:
`https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp`.

## ESO SkyCalc Moon and zodiacal-light tables

`eso_skycalc_rolo_moon_albedo.dat` is an unchanged ESO sky-model table of ROLO
wavelength-dependent lunar reflectance coefficients. The underlying lunar
model is attributed to Kieffer and Stone (2005). The predictor combines these
coefficients with the Meftah solar carrier and observing geometry to construct
the high-resolution lunar carrier.

`eso_skycalc_leinert_zodiacal_light.dat` is an unchanged ESO sky-model table of
zodiacal-light surface brightness as a function of heliocentric ecliptic
longitude and ecliptic latitude. The B500 table is attributed to Leinert et al.
The predictor interpolates this geometry table and applies the frozen colour,
extinction, velocity, and scale convention recorded in the model manifest.

The ESO Cerro Paranal Advanced Sky Model and SkyCalc documentation are the
software provenance for these local files. Preserve ESO, ROLO/Kieffer-Stone,
and Leinert attribution when redistributing them.

## PALACE tables

The canonical PALACE v1.0 data and code release is Noll et al. (2024), Zenodo
DOI `10.5281/zenodo.14064022`; the model description is Noll et al. (2025),
*Geoscientific Model Development*, 18, 4353-4398. PALACE data are published
under CC BY 4.0 and code under GPLv3.

The new method uses the line-specific binary coefficient asset below plus five
selected PALACE-compatible ASCII tables. This bundle also retains the preceding
OH tables for explicit backward-compatible runs:

- `palace_line_telluric_r4m_v1.fits`: the 11,552-row line catalogue used by
  the integrated decomposition, in exact model order. For every OH, atomic,
  oxygen-recombination, and O2 transition it stores the official PALACE v1.0
  line-specific `T_ref` and `fH2O` values and the equivalent dry and H2O
  reference optical depths. Noll et al. (2025) state that these line-specific
  values were derived from LBLRTM transmission spectra at maximum resolving
  power `4e6`; the source `palace_lines.fits` is byte-identical to Zenodo
  record `10.5281/zenodo.14064023`. Runtime line attenuation is evaluated
  directly from these rows and never interpolated from `palace_cont.fits`;

- `pmd_popmodel_OH_telluric_upper_parity_lsf_adam_25000_v1.dat`: current runtime
  OH table. It changes only `Aij` for 11,393 mapped PALACE rows using the frozen
  25,000-step telluric-aware Adam result. The optimization jointly varied
  upper-state/per-family parity weights and each spectrum's continuous 2-D LSF;
  it used 9 development spectra, excluded `expnum=45851`, and was checked with
  frozen OH weights on 3 separate validation spectra. Holdout was not evaluated;
- `pmd_popmodel_OH_h_family_default_ef_v1.dat`: previous runtime OH table with
  the latest export of the terminal 30,000-step refined upper-group/per-family
  weights and the PALACE e/f split retained. It reproduces the previously
  promoted corrected-decoder flat-family export to within text-rounding
  precision (maximum relative `Aij` difference `8.90e-13`). Its declared
  held-out [O I] 5577 local non-regression gate failed, and this status is
  retained in the manifest;
- `pmd_popmodel_OH_joint_v2_updated.dat`: preceding frozen repository OH
  population table, retained for backward-compatible explicit selection and
  recording the table used by the Moon/Zodi model training;
- `pmd_refcont_canonhyb_v1.dat`: **the default diffuse continuum table since
  2026-09-10.** Canonical PALACE v1.0 `fcHO2` and `fcFeO` interpolated onto the
  native LVM grid, with `fcO2Ac` taken verbatim from the native-LVM refit
  below. The canonical HO2/FeO vectors restore PALACE's species
  identification -- FeO peaks at 5966 A, matching the 595 nm FeO(VIS)
  component of Noll et al. (2024), and HO2 is correctly the blue tail of the
  1.51 um feature with only 4.8% of its emission below 9800 A. The refit
  O2Ac is kept because canonical O2Ac peaks at 3220 A, outside the LVM band,
  so only its tail is in range and it runs about twice too high through
  4200-5900 A. Measured on ten far-arm dark off-ecliptic rows, the fully
  canonical table costs a factor 1.43 in blue chi2 and biases the median
  residual to -0.26 sigma; this hybrid recovers that to -0.02 and has the
  best full-band chi2 of the three variants;
- `pmd_refcont_joint_native_adam_invsky_p2_10000iter.dat`: the previous
  default -- experimental native-LVM HO2, FeO, and O2Ac continuum export;
  its header records the unchanged grid, optimizer, source hash, and
  diagnostic status. Retained for backward-compatible explicit selection.
  Its HO2 vector was refit to a 595 nm-peaked shape, duplicating FeO rather
  than PALACE's near-IR species;
- `pmd_intdata_atom.dat`: canonical atomic-line/multiplet reference data;
- `pmd_intmodel_Orc.dat`: canonical oxygen-recombination line model;
- `pmd_popmodel_O2.dat`: canonical O2 population-model table used by the O2
  prefit and final `O2_b01` component.

The first three tables are repository-specific frozen derivatives, not
unmodified PALACE v1.0 products. Their scientific provenance is retained in
their headers and in `bundle_manifest.json`. The three remaining tables are
copied unchanged from the local PALACE v1.0 installation. The original files under
`skysub/palace/PMD/` remain in place for backward compatibility.

## Full-grid residual PCA

`residual_pca/palace_aijc_full_native_pca_v1.npz` contains the residual mean and
the first 30 orthonormal components from 1,000 PALACE-`Aijc` fits. The matrix was
built on all 12,401 native wavelength pixels in float64 with column-mean
subtraction only; `expnum=45851` was excluded. The 10- and 20-component models
are exact prefixes of the same basis, rather than duplicated assets.

Use the basis only with the PALACE `Aijc` OH strengths from which its residuals
were constructed:

```python
from skysub.sky_decomp.residual_pca import (
    SkyDecompTelluricCorrectedLinesResidualPCA,
)

decomposer = SkyDecompTelluricCorrectedLinesResidualPCA(
    wave,
    telluric_calculator=telluric_calculator,
    pwv_mm=pwv_mm,
    source_airmass=source_airmass,
    drp_transmission=drp_transmission,
    n_residual_pca_components=30,  # or 10/20
)
```

## Line-amplitude residual PCA

`residual_pca/palace_aijc_line_amplitude_pca_v1.npz` is the separate
line-space alternative. For every one of the same 1,000 PALACE-`Aijc` fits,
the full 12,401-pixel residual was fitted without sign constraints using the
saved final LSF and 11,552 separate PALACE catalogue transitions. The group ID
is metadata only: every OH parity/branch transition is an independent signed
coefficient. Each sparse column uses the exact telluric-corrected, integrated
continuous-2D-LSF profile and is normalized to unit native-grid integral, so
the coefficient is an observed integrated residual-line flux. A weak,
scale-invariant ridge (`lambda=1e-4`) stabilizes unresolved blends. Of the
11,552 catalogue transitions, 11,538 have line centres and nonzero support on
the native grid; the 14 edge transitions remain explicit zero columns. PCA
uses per-amplitude column-mean subtraction only.

At prediction time the baseline decomposition and its five LSF cycles run
first. The selected amplitude components are then projected through that
spectrum's final individual-line design and fitted to its residual with the
LSF fixed:

```python
from skysub.sky_decomp.residual_pca import (
    SkyDecompTelluricCorrectedLinesLineAmplitudePCA,
)

decomposer = SkyDecompTelluricCorrectedLinesLineAmplitudePCA(
    wave,
    telluric_calculator=telluric_calculator,
    pwv_mm=pwv_mm,
    source_airmass=source_airmass,
    drp_transmission=drp_transmission,
    n_line_amplitude_pca_components=30,  # or 10/20
)
```

Both PCA methods retain the compact fitted LSF in the result's `lsf_state`.
`results_to_fits` writes it as `LSF_COEF`, `LSF_KNOTS`, and `LSF_META`; use
`load_lsf_surface_state` to reconstruct it without refitting.

## VN-grouped OH and individual-line amplitude PCA

`SkyDecompTelluricCorrectedLinesPalaceAijcVN` ties OH transitions only by
upper `(v_upper, N_upper)`, reducing the OH fit from 357 to 188 groups. The
relative `F_upper`, branch, and e/f parity strengths inside each VN group stay
fixed at PALACE `Aijc * gi`; the 11,393 OH transitions themselves are unchanged.
With the five non-OH atomic groups, two oxygen-recombination groups, and one O2
group, the emission-line fit has 196 non-negative coefficients.

`residual_pca/palace_aijc_vn_line_amplitude_pca_v1.npz` contains the first 30
PCA components from signed amplitudes of all 11,552 individual catalogue lines
fitted to the full-grid VN residuals. The frozen corpus contains 1,000 spectra;
`expnum=45851` was excluded before selection and `expnum=45858` was removed by
the recorded 7-robust-sigma RMS rule, leaving 999 PCA-training spectra. The
cumulative amplitude-space variance fractions are 0.6581, 0.7234, and 0.7569
for 10, 20, and 30 components.

At prediction time, the VN baseline first determines the continuous 2-D LSF.
The selected signed amplitude components are projected through that exact
telluric-plus-LSF individual-line design. The physical VN model and projected
PCA columns are then solved together in one final linear fit while the LSF is
held fixed:

```python
from skysub.sky_decomp.residual_pca import (
    SkyDecompTelluricCorrectedLinesVNLineAmplitudePCA,
)

decomposer = SkyDecompTelluricCorrectedLinesVNLineAmplitudePCA(
    wave,
    telluric_calculator=telluric_calculator,
    pwv_mm=pwv_mm,
    source_airmass=source_airmass,
    drp_transmission=drp_transmission,
    n_line_amplitude_pca_components=30,  # or 10/20
)
```

## VNF coefficient PCA Prep and final line-amplitude PCA

`palace_aijc_vnf_coefficient_pca_prep_v1.npz` compresses the 357 non-negative
PALACE-`Aijc` VNF OH coefficients from 994 spectra retained after the recorded
7-robust-sigma cut. Eighteen signed PCs plus one non-negative mean scale explain
0.9999991015 of the centered training variance. The fit expands those latent
coefficients back to all 357 physical VNF amplitudes before each LSF update, so
the LSF is still inferred from the complete line catalogue.

`palace_aijc_vnf_coefficient_line_amplitude_pca_v1.npz` is trained from signed,
lightly ridge-regularized fits of all 11,552 individual line amplitudes to the
full-grid Prep residuals. Its first 200 PCs are stored. Fifty are the production
default selected from the held-out diminishing-returns checkpoint. The stricter
0.999999 amplitude-space variance target would require 993 PCs, and the asset
records that failure explicitly instead of presenting the 50-PC default as an
exact amplitude reconstruction.

```python
from skysub.sky_decomp.residual_pca import (
    SkyDecompTelluricCorrectedLinesPalaceAijcVNFPCAPrep,
    SkyDecompTelluricCorrectedLinesVNFPCALineAmplitudePCA,
)

prep = SkyDecompTelluricCorrectedLinesPalaceAijcVNFPCAPrep(
    wave,
    telluric_calculator=telluric_calculator,
    pwv_mm=pwv_mm,
    source_airmass=source_airmass,
    drp_transmission=drp_transmission,
)
final = SkyDecompTelluricCorrectedLinesVNFPCALineAmplitudePCA(
    wave,
    telluric_calculator=telluric_calculator,
    pwv_mm=pwv_mm,
    source_airmass=source_airmass,
    drp_transmission=drp_transmission,
    n_line_amplitude_pca_components=50,
)
```

Both results retain every fitted coefficient and the compact per-spectrum B/R/Z
LSF spline state. `results_to_fits` writes them through the existing `COEF`,
`LSF_COEF`, `LSF_KNOTS`, and `LSF_META` extensions without changing older output
schemas.

## Production split-zodi VNF PCA30

`palace_aijc_vnf_split_zodi_line_amplitude_pca30_v1.npz` is the production
PCA30 basis for the merged method. It was trained from 1,000 successful
PALACE-Aijc VNF fits using the production split Moon/Zodiacal and diffuse
continuum profile, each spectrum's telluric transmission and continuous LSF,
and the native 12,401-pixel grid. The recorded robust RMS cut retained 983
spectra. Thirty signed components explain 0.99972595 of the centered fitted
line-amplitude variance. The unchanged binary retains its original training
metadata for provenance; the bundle manifest records both its legacy and
production identifiers.

```python
from skysub.sky_decomp.residual_pca import (
    SkyDecompPalaceAijcVNFSplitZodiLineAmplitudePCA30,
)

decomposer = SkyDecompPalaceAijcVNFSplitZodiLineAmplitudePCA30(
    wave,
    telluric_calculator=telluric_calculator,
    pwv_mm=pwv_mm,
    source_airmass=source_airmass,
    drp_transmission=drp_transmission,
)
```

## Direct line-adjoint PCA

`palace_aijc_vnf_line_adjoint_pca_v1.npz` is the no-full-dictionary-fit
alternative trained on 994 retained residuals from the 1,000-spectrum
PALACE-`Aijc` 357-group VNF corpus. For spectrum `i`, its saved
telluric-plus-LSF line operator `A_i` maps the native residual
directly to diagonal matched-filter amplitudes,
`x_i = diag(A_i.T @ A_i)^-1 @ A_i.T @ residual_i`. No coupled 11,552-parameter
amplitude solve is performed. Consequently unresolved blends remain cross-talk
in `x_i`; these coordinates must not be described as deblended line fluxes.

The frozen float64 asset stores 500 PCs, which are also the held-out-selected
default. They explain 0.9975226 of the centered adjoint-amplitude variance. At
prediction time the vectors are projected through the fitted spectrum's exact
individual telluric-plus-LSF operator and joined to the physical VNF-PCA model
in one final linear solve with signed PCA coefficients:

```python
from skysub.sky_decomp.residual_pca import (
    SkyDecompTelluricCorrectedLinesPalaceAijcVNFLineAdjointPCA,
)

decomposer = SkyDecompTelluricCorrectedLinesPalaceAijcVNFLineAdjointPCA(
    wave,
    telluric_calculator=telluric_calculator,
    pwv_mm=pwv_mm,
    source_airmass=source_airmass,
    drp_transmission=drp_transmission,
)
```

The fitted result records these columns as `LineAdjointPCA_*`; FITS output uses
`DECOMPM=telluric-corrected-lines-palace-aijc-vnf-line-adjoint-pca` and
`LADPCAK=500`. The compact LSF output is unchanged.

Across the four frozen held-out spectra, the median full-grid RMS is 0.14796,
versus 0.32445 for the source VNF fit. The coupled full-dictionary amplitude
PCA remains more accurate (0.08344 at 50 PCs); the adjoint method exists to
avoid that training solve, not as a claim of superior residual reconstruction.

## Runtime selection

The new method uses this bundle without a caller-supplied `base_dir`:

```python
decomposer = SkyDecompMoonZodiLSFSurfaceIterative(
    wave,
    physical_to_fit_flux_scale=1e14,
)
```

The parallel split-zodi mode uses the same frozen OH and diffuse tables when no
legacy positional PALACE root is supplied:

```bash
python skysub/decompose_parallel.py input.fits \
    --fit-model lsf-surface-iterative-split-zodi
```

Both bundled command-line modes select
`pmd_popmodel_OH_telluric_upper_parity_lsf_adam_25000_v1.dat` and
`pmd_refcont_canonhyb_v1.dat` by default. Explicit `--palace-oh-suffix` and
`--palace-diffuse-suffix` values still override those defaults; pass
`--palace-oh-suffix _h_family_default_ef_v1` to restore the pre-telluric OH
table, or pass
`--palace-diffuse-suffix _joint_native_adam_invsky_p2_10000iter` to restore
the pre-2026-09-10 diffuse table.

The split-zodi mode also applies a **diffuse species-ratio bracket** by
default (`--diffuse-ratio-bound-dex 0.2`, `--diffuse-ratio-nominal
0.0396,0.7026,0.2578`): the three diffuse species are individually
unidentifiable in the LVM band, and the three arms of one exposure disagree
about `log10(FeO/HO2)` by 0.633 dex at the median when the ratios are free.
Pass `--diffuse-ratio-bound-dex 0` to disable it. The nominal is FLUX shares
measured on the corpus being fitted, not PALACE's own reference shares --
see `decompose_parallel.SPLIT_ZODI_DIFFUSE_RATIO_NOMINAL` for why, and
re-measure it if the basis or corpus changes.

For a remote clone, the packaged root can be selected explicitly from the Git
root. Omitting the suffix flags intentionally follows the versions declared by
the checked-out code and bundle manifest instead of pinning an older filename:

```bash
python skysub/decompose_parallel.py input.fits \
    --moon-zodi-data-root skysub/sky_decomp/data \
    --fit-model lsf-surface-iterative-split-zodi \
    --n-spline-knots 11 \
    --n-zodi-spline-knots 1 \
    --n-refinement-cycles 5 \
    --n-workers 32 \
    --output-dir moon_zodi_spline/
```

The optional positional root is also accepted when it points at this complete
bundle, but the named `--moon-zodi-data-root` form is unambiguous and preferred
for scripts.

An externally distributed copy of the complete `data/` directory can be used
without changing repository files:

```python
decomposer = SkyDecompMoonZodiLSFSurfaceIterative(
    wave,
    physical_to_fit_flux_scale=1e14,
    data_root="/shared/lvmsky/sky_decomp_data",
)
```

The supplied root must contain both `moon_zodi/` and `palace/PMD/`. A missing
or checksum-mismatched required file raises a specific exception; there is no
network download and no fallback to the historical solar-only Moon model.

## Packaging and redistribution

If binary data are committed, package configuration must include this complete
directory. If binaries are distributed separately, keep this README and the
two manifests with the code, publish one immutable bundle archive, and verify
all SHA-256 values after extraction. Do not distribute the 494 MB assessment
run directory as part of this production bundle.
