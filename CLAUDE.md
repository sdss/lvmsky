# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Research code for sky subtraction in SDSS-V LVM (Local Volume Mapper) spectra. Active development is almost entirely in `skysub/`. The other top-level directories are vendored or legacy:

- `skycorr/`, `skymodel/`: vendored ESO C code (Skycorr, Cerro Paranal Sky Model, autotools builds). Treat as third-party and leave alone unless asked.
- `telluric/`, `notebooks/`, `scratch/`: older analyses and one-off scripts.
- `model_tuning_analysis.md`, `moon_zodi_split_priors_2026-09-03.md`: dated design notes that record why priors and tuning are set the way they are. Read them before changing fit priors.
- `skysub/sky_subtraction_method.tex` / `.bib`: the paper-style method description. `skysub/publication_plots.ipynb` makes its figures.

There is no top-level package metadata, requirements file, or CI. The only installable sub-project is `skysub/medians_computation` (`lvm-medians`).

## Pipeline (big picture)

The pipeline has three stages, and each stage's output is the next stage's input:

1. **Median stacks** (`skysub/medians_computation`, CLI `lvm-medians`): scans LVM SFrames, optionally caches Gaia DR3 matches to reject star-contaminated fibers, and writes median spectra per exposure and arm (SCI, SkyE, SkyW) with a META table.
2. **Decomposition** (`skysub/sky_decomp` + driver `skysub/decompose_parallel.py`): fits each spectrum as a non-negative additive model (OH lines from PALACE tables, atomic airglow, O2, a Moon/Zodi continuum split into two B-spline families × solar spectrum, a diffuse HO2/FeO/O2Ac continuum, and ORC) with a curvature-penalised QP solved by Clarabel. It refines the LSF iteratively and applies per-row telluric transmission. The output is a coefficient FITS with META, coefficients, LSF HDUs, and reliability flags.
3. **Prediction** (`skysub/mlp_predictor`): a PyTorch ensemble (`DualEncoderGroupHeadMLPCompressed`) maps the two sky-arm coefficient vectors plus geometry/time context (39 features) to the science-arm coefficients. Then two post-prediction corrections run on the reconstructed spectrum:
   - `sky_arm_correction.py` transfers the sky arms' own decomposition residual (for example, solar absorption features) to the science prediction.
   - `sky_line_scaling.py` rescales OH bands and atomic-line templates on the science spectrum itself. It is designed to be robust to science continuum and science-line leakage.

`skysub/notebook_example_predict_sky.ipynb` is the end-to-end reference for production use. It decomposes one exposure via `decompose_parallel.decompose_in_process` so that the coefficients are produced *exactly* as the training corpus was. Keep that invariant: the predictor's inputs must come from the same fit model and defaults as the training corpus.

### Decomposition class hierarchy

`SkyDecomp` (`fit.py`) → `SkyDecompLSFSurfaceIterative` (`lsf_surface_iterative.py`) → `SkyDecompLSFSpline2D` (`lsf_spline2d.py`) → telluric line classes (`telluric_corrected_lines.py`) → PALACE-Aijc / VNF / PCA variants (`residual_pca.py`). Subclasses mostly change the design matrix or line tables. Shared fitting machinery lives in the base classes, so a change in `fit.py` or `lsf_surface_iterative.py` affects every model.

`decompose_parallel.py` exposes only the telluric fit models (`TELLURIC_FIT_MODELS`). The default and deployed model is `palacecorr-aijc-vnf-split-zodi-lsf-spline2d` (SkyFar ridge-corrected PALACE OH strengths, `_palacecorr_` output suffix). The non-telluric models still exist only because they are base classes. To reproduce pre-2026-09-17 corpora, check out an older commit.

`reliability.py` defines the per-row flags (moon/zodi reversal, diffuse collapse, etc.) once, for both the fitter that writes them and the analyses that read them. It is numpy-only on purpose. `result_io.py` defines the FITS HDU layout.

### Data bundle

`skysub/sky_decomp/data/` is the portable data root (PALACE tables under `palace/PMD/`, the Moon/Zodi physical assets, the residual-PCA bases, and the sensitivity curves). `bundle_manifest.json` and `moon_zodi/model_manifest.json` record SHA-256 digests, and the fitters validate them before running. If you replace a data file, update its manifest. `decompose_parallel` also fingerprints the input, data manifest, code, and parameters for its resumable per-row NPZ cache. An identical rerun resumes, and a code change invalidates the cache.

## Running things

There are two import conventions, and both are in use:

- `sky_decomp` and the driver use fully qualified `skysub.…` imports. Run them from the **repo root**, for example `python -m skysub.decompose_parallel <stack.fits> --n-workers 8 --output-dir out/`. Running `python skysub/decompose_parallel.py …` also works because it patches `sys.path`.
- `mlp_predictor` and the notebooks import `mlp_predictor`, `sky_decomp`, and `decompose_parallel` as top-level modules. Run notebooks with the working directory set to `skysub/`. Tests for `mlp_predictor` code insert `skysub/` into `sys.path` themselves.

Tests (pytest; from repo root):

```bash
python -m pytest skysub/sky_decomp/tests                       # decomposition + predictor-correction tests
python -m pytest skysub/sky_decomp/tests/test_sky_line_scaling.py::<test_name>
python -m pytest skysub/experiments_moon_zodi_baseline_decomp_jax_v1/tests
cd skysub/medians_computation && PYTHONPATH="$PWD/src" python -m pytest   # lvm-medians
```

Golden-file regression tests (`test_lsf_surface_iterative_golden.py`, data in `sky_decomp/tests/data/`) pin numerical fit output. If an intentional change breaks one, regenerate the golden file deliberately and say so. Don't loosen the tolerances.

Environment: none of the conda environments on this machine currently has the full stack (`torch` and `pytest` are missing). The upstream author runs `conda run -n lvmdrp_dev_311 …` (Python 3.11). The core dependencies are numpy, scipy, astropy, clarabel (the QP solver), torch (`mlp_predictor`), and tqdm. jax is needed only for the jax experiment. `lvm-medians` linting: `ruff` with line length 100.

`decompose_parallel.py` forces all BLAS/OpenMP/Rayon thread pools to one thread per worker before importing numpy. Keep that ordering if you edit the imports.

## Conventions

- Experiments live in `skysub/experiments_*_v1/` as self-contained modules (often with their own README, outputs, and tests). Production code lives in `sky_decomp/` and `mlp_predictor/`. Some experiment READMEs reference modules that no longer exist on this branch (for example, `sky_decomp/niv_continuum.py`).
- Code comments often carry dated rationale (for example, `# 2026-09-24: …`) that explains why a default was chosen, with the measured effect. Preserve these comments, and add one in the same style when you change a default.
- `mlp_predictor/diagnostics_cells/` holds one notebook-cell-sized diagnostic per file. `config.PipelineConfig()` defaults to an older corpus, so callers set the corpus explicitly.
