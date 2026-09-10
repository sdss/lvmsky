# Sky decomposition data bundle

This directory is the portable data root for
`SkyDecompMoonZodiLSFSurfaceIterative`. It contains the physical Moon/Zodi
model, the exact PALACE tables used by the validated fit, and an unchanged copy
of the Meftah SOLAR-HRS text spectrum used by the historical decomposition.

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
└── palace/PMD/
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

The new method needs five selected PALACE-compatible ASCII tables. This bundle
also retains the preceding OH table for explicit backward-compatible runs:

- `pmd_popmodel_OH_h_family_default_ef_v1.dat`: current runtime OH table with
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
`pmd_popmodel_OH_h_family_default_ef_v1.dat` and
`pmd_refcont_canonhyb_v1.dat` by default. Explicit `--palace-oh-suffix` and
`--palace-diffuse-suffix` values still override those defaults; pass
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
