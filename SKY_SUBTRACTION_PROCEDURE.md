# Sky subtraction procedure

Recorded on 2026-11-16. The user reported running the scan and median-build
commands below on Utah. This record does not assert their processing counts or
output quality. The decomposition command is the next step and has not been
run as part of this work.

## Runtime and working directory

Use the configured Python 3.11 scientific environment on Utah. The module
commands below run from the repository root:

```bash
cd /uufs/chpc.utah.edu/common/home/u6058164/lvmsky
```

Before decomposition, copy the shared-LSF code changes and the complete
`skysub/sky_decomp/data/lsf/` CSV calibration bundle to this checkout. The
decomposition command also runs from the repository root. Data, calibration,
work and output paths are absolute.

## 1. Scan the Knox SFrames (user-reported run)

```bash
python -m skysub.medians_computation.src.lvm_medians \
  --work-dir /uufs/chpc.utah.edu/common/home/sdss/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2_knox \
  scan \
  --sframes-root /uufs/chpc.utah.edu/common/home/sdss/sdsswork/lvm/spectro/redux/knox
```

This recursively finds and sorts `lvmSFrame-*.fits` and `.fits.gz` paths. It
creates the work directory and writes its `sframes.txt`; it does not compute
spectra or download Gaia data.

## 2. Build medians (user-reported run)

```bash
python -m skysub.medians_computation.src.lvm_medians \
  --work-dir /uufs/chpc.utah.edu/common/home/sdss/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2_knox \
  build-medians \
  --gaia-fibers-dir /uufs/chpc.utah.edu/common/home/sdss/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2/gaia/fibers \
  --workers 20 \
  --gaia-ratio-threshold 0.01 \
  --output /uufs/chpc.utah.edu/common/home/sdss/sdsswork/users/u6058164/lvmsframe_median_stack_1.3.2_knox_gaia1over100.fits
```

The command reads `sframes.txt` from the work directory and reuses individual
Gaia fiber caches. The `gaia1over100` selection is Gaia/LVM G-band flux ratio
at most 0.01; these ratios come from the existing cache rather than a new
measurement on Knox spectra. The median product includes `FIBER_ROLE`, which
records the physical fibers contributing to Science, Sky near and Sky far.

These commands do not apply the Knox header-quality criterion
`FMDIRECT == True`, `FMDRMODE == 'oh+bz'`, `FMQR1 < 1`, `FMQZ1 < 1`.
That selection must be established separately before using the resulting
decompositions as a training corpus.

## 3. Decompose Science and both Sky medians (next step)

```bash
python skysub/decompose_parallel.py \
  /uufs/chpc.utah.edu/common/home/sdss/sdsswork/users/u6058164/lvmsframe_median_stack_1.3.2_knox_gaia1over100.fits \
  --fit-model palacecorr-aijc-vnf-split-zodi-lsf-spline2d-solarblue \
  --blue-lsf reference \
  --lsf-calibration /uufs/chpc.utah.edu/common/home/u6058164/lvmsky/skysub/sky_decomp/data/lsf \
  --n-workers 20 \
  --pin-workers \
  --output-dir /uufs/chpc.utah.edu/common/home/sdss/sdsswork/users/u6058164/sky_decomposition_1.3.2_knox/
```

One invocation fits all three roles: `sci`, `sky1` (near), and `sky2` (far).
The shared Sky Subtraction code loads per-fiber CSV profiles and averages
the measured profiles selected by each row's `FIBER_ROLE`. Blue remains fixed;
R/Z start from their selected-fiber CSV means and refine on the target spectrum.
The default refinement limit is 30 cycles, with earlier convergence permitted.

`--pin-workers` binds each worker process to one allowed logical CPU on Linux
using `sched_setaffinity`. With 20 available CPUs, the 20 workers are assigned
different CPUs. The flag is optional; removing it lets the operating system
schedule workers across the allowed CPUs. It does not change the fit settings.

The output directory is created automatically. After a successful ordinary
run (without `--compact-only`), expect 11 FITS files. In the names below:

```text
S = lvmsframe_median_stack_1.3.2_knox_gaia1over100
M = palacecorr_aijc_vnf_split_zodi_lsf_spline2d_solarblue
```

| Filename pattern | Count | Contents |
| --- | ---: | --- |
| `S_decomp_{sci,sky1,sky2}_M.fits` | 3 | Full decomposition for each role |
| `S_{sci,sky1,sky2}_meta_coef_M.fits` | 3 | Compact coefficients, LSF states and fit metadata |
| `S_meta_only.fits` | 1 | Input metadata table |
| `S_every10_decomp_{sci,sky1,sky2}_M.fits` | 3 | Every tenth row of each full decomposition |
| `S_every10.fits` | 1 | Every tenth row of the input median stack |

The braces denote three separate filenames. The original full median stack
stays at its input path. Use a new output directory for a different
calibration or method. This is decomposition only; it neither trains the MLP
nor writes sky-subtracted SFrames.

## Sky Subtraction interface compatibility

The shared-LSF change preserves the public `run(...)` signatures of
`skysub.sky_subtraction.pipeline` and `skysub.sky_subtraction_block.pipeline`.
The Sky Subtraction command-line entry points and options are unchanged.
Internal LSF loading and averaging moved into `sky_decomp.lsf_reference`;
Sky Subtraction retains its existing CSV reference behavior.

See [CHANGELOG.md](CHANGELOG.md) for implementation and validation details.
