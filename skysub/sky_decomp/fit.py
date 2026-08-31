"""LVM sky-decomposition weighted-NNLS fit engine.

`SkyDecomp` builds a non-negative additive spectral model of one arm's flux
against a physically motivated design matrix and solves a curvature-penalised
quadratic program (QP) via Clarabel.  Public families in the default fit are:

    OH lines  +  Moon_bs (continuum spline * solar)  +  diffuse (HO2 + FeO + O2Ac)
    +  atomic airglow lines  +  ORC  +  O2

Two subclasses in this package add per-row iterative LSF refinement:

    `SkyDecompLSFSurfaceIterative` — smooth wavelength-dependent LSF surface,
        used for the p40_p70 corpus decompositions.
    `SkyDecompMoonZodiLSFSurfaceIterative` — same wrapper but with a fully
        physical Moon+Zodi predictor replacing the color-envelope spline.

`decompose_parallel.py` runs any of these at corpus scale via
`--fit-model {baseline, lsf-surface-iterative, lsf-surface-iterative-split-zodi,
moon-zodi-lsf-surface-iterative}`.

split_zodi mode
---------------

`SkyDecomp(split_zodi=True)` and the same flag on `SkyDecompLSFSurfaceIterative`
promote the single continuum spline `Moon_bs` into two color-tagged families:

    Moon_bs = solar_rb(λ) * moon_albedo(λ) * B_K_moon(λ) * c_moon[k]
    Zodi_bs = solar_rb(λ) * zodi_color(λ) * B_K_zodi(λ) * c_zodi[k]

where

    solar_rb(λ)     — LATMOS solar SED rebinned and LSF-convolved onto the LVM
                      wavelength grid (already used by the baseline `Moon_bs`).
    moon_albedo(λ)  — ROLO lunar disk-integrated albedo at a fiducial phase
                      angle (default 30°), median-normalised to unity.
                      Loaded from `moon_albedo_asset_path`, default
                      `<base_dir>/moon_zodi/eso_skycalc_rolo_moon_albedo.dat`.
    zodi_color(λ)   — Leinert-style reddening `(λ/5000Å)^zodi_color_exponent`,
                      default exponent 0.26 (matches `moon_zodi_model.py`),
                      median-normalised to unity.
    B_K(λ)          — cubic B-spline design on `K` uniformly-spaced interior
                      knots.  Default `n_zodi_spline_knots=3` (the zodi spline
                      is intentionally much smoother than the moon spline).

Both families remain non-negative in the QP.  The two color envelopes are
what makes them distinguishable: moon_albedo carries sharp ROLO mineral
bands while zodi_color is a smooth power-law reddening.

Curvature penalties
~~~~~~~~~~~~~~~~~~~
Each family carries a second-difference smoothness penalty:

    moon: ‖D² c_moon‖² * moon_smooth_lambda    (default 0.1)
    zodi: ‖D² c_zodi‖² * zodi_smooth_lambda    (default = 10 * moon)

Bump `zodi_smooth_lambda` (~1e-1) when using a small zodi knot count (K≤3);
a tighter zodi penalty makes the fitted zodi cleanly follow the red-slope
`zodi_color` envelope instead of picking up locally blue features.

Recommended settings
~~~~~~~~~~~~~~~~~~~~
Validated on the p40_p70 every10 corpus (100-row phase-stratified sample):

    n_spline_knots        = 29       (unchanged Moon_bs count)
    n_zodi_spline_knots   = 3
    moon_smooth_lambda    = 1e-3     (baseline value; kept for continuity)
    zodi_smooth_lambda    = 1e-1
    moon_albedo_fiducial_phase_deg = 30.0
    zodi_color_exponent   = 0.26

Backwards compatibility
~~~~~~~~~~~~~~~~~~~~~~~
`split_zodi=False` is the default and produces the exact same design matrix,
components, and solver output as prior versions of this module.  When
`split_zodi=True` the result gains a `"zodi"` component (and matching entry
in `_components_sigma_from_coef_err`); `SkyDecompResult` grows two kw-only
fields `zodi_names` and `zodi_knots`.  83 pre-existing tests pass unchanged.

Documentation and validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The design, identifiability analysis, and end-to-end validation on 100
phase-stratified rows live in
`notebooks/moon_zodi_split_identifiability.ipynb`.  Corpus-scale usage:

    python skysub/decompose_parallel.py <every10.fits> <palace_dir> \\
        --fit-model lsf-surface-iterative-split-zodi \\
        --n-zodi-spline-knots 3 --zodi-smooth-lambda 0.1 \\
        --n-refinement-cycles 5 --n-workers 8 --output-dir <out>

produces `..._decomp_{sci,sky1,sky2}_lsf_surface_iterative_split_zodi.fits`
each carrying a new `COMP_ZODI` HDU alongside the standard component HDUs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import pathlib
import re
import time
import tracemalloc

import clarabel
import numpy as np
import scipy.sparse as sp
from astropy.table import Table
from scipy.interpolate import BSpline
from scipy.optimize import least_squares

from lvmdrp.core.fluxcal import rebin_and_convolve

CAP_WAVE = 5.0
HC_OVER_KB_CMK = 1.4387769
OH_GROUP_KEYS = ("v_upper", "N_upper", "F_upper")
O2_MIN = 8600.0
O2_MAX = 8696.0
T_O2_REF = 191.5
T_O2_HALF_RANGE = 20.0
O2_MIN_VALID_FRAC = 0.5
SUFFIXED_PMD_TABLES = frozenset({"pmd_popmodel_OH.dat", "pmd_refcont.dat"})
LSF_KERNEL_SIZE = 11
LSF_MIN_VALID_FRAC = 0.5
LSF_CHANNELS = (
    ("B", None, 5787.0),
    ("R", 5787.0, 7454.0),
    ("Z", 7454.0, None),
)


def vac_to_air(lam_vac_a: np.ndarray) -> np.ndarray:
    lam = np.asarray(lam_vac_a, float)
    s2 = (1e4 / lam) ** 2
    n = 1.0 + 8.34254e-5 + 2.406147e-2 / (130.0 - s2) + 1.5998e-4 / (38.9 - s2)
    return lam / n


def decode_hitran_id(table: Table) -> Table:
    ids = np.asarray(table["ID"].astype(str), dtype=str)
    if np.min(np.char.str_len(ids)) < 13:
        raise ValueError("HITRAN ID strings are shorter than expected.")

    v_up = np.where(np.array([s[4:5] for s in ids]) == "X", "10", np.array([s[4:5] for s in ids])).astype(int)
    v_low = np.array([s[5:6] for s in ids], dtype=int)
    branch_n = np.array([s[6:7] for s in ids])
    branch_j = np.array([s[7:8] for s in ids])
    f_up = np.array([s[8:9] for s in ids], dtype=int)
    f_low = np.array([s[9:10] for s in ids], dtype=int)
    n_low = np.array([s[10:12] for s in ids], dtype=int)
    parity = np.array([s[12:13] for s in ids])

    delta_map = {"O": -2, "P": -1, "Q": 0, "R": 1, "S": 2}
    n_up = n_low + np.vectorize(delta_map.get)(branch_n)

    table["v_upper"] = v_up
    table["v_lower"] = v_low
    table["N_upper"] = n_up
    table["N_lower"] = n_low
    table["F_upper"] = f_up
    table["F_lower"] = f_low
    table["branch_N"] = branch_n
    table["branch_J"] = branch_j
    table["parity"] = parity
    table["is_main"] = (f_up == f_low) & np.isin(branch_n, ["P", "Q", "R"])
    return table


def grp2vector(
    line_wave: np.ndarray,
    line_amp: np.ndarray,
    wave: np.ndarray,
    lsf_sigma: np.ndarray | float,
) -> np.ndarray:
    cent = np.asarray(line_wave, float)
    amp = np.asarray(line_amp, float)
    sig = np.interp(cent, wave, lsf_sigma) if np.ndim(lsf_sigma) > 0 else float(lsf_sigma)
    yy = (wave[:, None] - cent) / sig
    return np.sum(amp[None, :] * np.exp(-0.5 * yy**2), axis=1)


def sticks2vector(
    line_wave: np.ndarray,
    line_amp: np.ndarray,
    wave: np.ndarray,
) -> np.ndarray:
    line_wave = np.asarray(line_wave, float)
    line_amp = np.asarray(line_amp, float)
    wave = np.asarray(wave, float)
    out = np.zeros_like(wave, dtype=float)
    if wave.size < 2 or line_wave.size == 0:
        return out

    inside = (line_wave >= wave[0]) & (line_wave <= wave[-1])
    if not np.any(inside):
        return out
    line_wave = line_wave[inside]
    line_amp = line_amp[inside]

    hi = np.searchsorted(wave, line_wave, side="left")
    inner = (hi > 0) & (hi < wave.size)
    if np.any(inner):
        i1 = hi[inner]
        i0 = i1 - 1
        span = wave[i1] - wave[i0]
        frac = np.divide(line_wave[inner] - wave[i0], span, out=np.zeros_like(i0, dtype=float), where=span > 0)
        frac = np.clip(frac, 0.0, 1.0)
        np.add.at(out, i0, line_amp[inner] * (1.0 - frac))
        np.add.at(out, i1, line_amp[inner] * frac)

    np.add.at(out, 0, np.sum(line_amp[hi == 0]))
    np.add.at(out, wave.size - 1, np.sum(line_amp[hi >= wave.size]))
    return out


def _build_d2_operator(n_par: int) -> np.ndarray:
    """Second-difference operator D such that (D c)_i = c_i - 2 c_{i+1} + c_{i+2}."""
    if n_par < 3:
        return np.zeros((0, n_par), dtype=float)
    d2 = np.zeros((n_par - 2, n_par), dtype=float)
    idx = np.arange(n_par - 2)
    d2[idx, idx] = 1.0
    d2[idx, idx + 1] = -2.0
    d2[idx, idx + 2] = 1.0
    return d2


@dataclass(slots=True)
class SkyDecompResult:
    coef: np.ndarray
    # Per-coefficient 1σ posterior uncertainty in the same units and order as
    # ``coef``. Computed from the active-set Fisher information of the final
    # weighted-NNLS solve (see ``SkyDecomp._coef_err_active_set`` for the full
    # derivation): the (regularised) Hessian sub-block of the strictly positive
    # coefficients is inverted, the internal column/data scaling is undone, and
    # the variance is inflated by ``max(reduced_chi², 1)`` before taking the
    # square root. Boundary coefficients pinned at c = 0 by the c ≥ 0 constraint
    # get ``NaN`` (no symmetric interval is defined there).
    coef_err: np.ndarray
    bestfit: np.ndarray
    resid: np.ndarray
    resid_level: float
    fit_status: str
    fit_summary: str
    reduced_chi2: float
    fit_elapsed_sec: float
    components: dict[str, np.ndarray]
    design_names: list[str]
    t_o2: float
    t_o2_err: float
    r2: float
    rms_resid: float
    peak_memory_mb: float
    o2_fit_status: str
    o2_fit_summary: str
    o2_fit_elapsed_sec: float
    o2_valid_frac: float
    lsf_kernels: dict[str, np.ndarray]
    lsf_metrics: dict[str, dict[str, object]]
    bestfit_lsf: np.ndarray
    moon_knots: np.ndarray
    moon_boosted_pixels: np.ndarray
    # Unit-integrated O2 template on `wave` and the prefit amplitude that
    # scales it; the fitted `coef['O2_b01']` absorbs the amplitude, so
    # `coef * vector_o2` reproduces the model O2 contribution.
    vector_o2: np.ndarray
    o2_prefit_amp: float
    # Per-pixel 1σ of the LSF-refined reconstruction, obtained by first-order
    # (Jacobian) propagation of ``coef_err`` through the LSF-convolved design
    # matrices used to build ``bestfit_lsf`` (see
    # ``SkyDecomp._components_sigma_from_coef_err``).  Same shape as
    # ``bestfit_lsf``.  Written to disk as the ``FLUX_SIGMA_TOTAL`` HDU by the
    # parallel decomposition writer.
    bestfit_lsf_sigma: np.ndarray
    # Zodi family (only populated when `SkyDecomp(split_zodi=True)`).  kw_only so
    # LSFSurfaceIterativeResult can still add positional fields after us.
    zodi_names: list[str] = field(default_factory=list, kw_only=True)
    zodi_knots: np.ndarray = field(default_factory=lambda: np.array([], float), kw_only=True)
    # Full posterior covariance sub-matrices for the moon and (optionally) zodi
    # B-spline blocks -- 2D arrays with shape (n_moon, n_moon) / (n_zodi, n_zodi)
    # in the same physical units as ``coef``. Entries corresponding to inactive
    # (boundary-zero) coefficients are ``NaN`` (no symmetric interval defined).
    # These carry the off-diagonal correlations that ``coef_err`` (= sqrt of the
    # diagonal) discards; downstream code that needs the joint uncertainty on
    # any linear combination of moon or zodi coefficients (e.g. total moon
    # amplitude, PCA scores) should use these matrices instead of assuming
    # independence.  ``None`` when the block was not present in the solve.
    coef_cov_moon: np.ndarray | None = field(default=None, kw_only=True)
    coef_cov_zodi: np.ndarray | None = field(default=None, kw_only=True)


class SkyDecomp:
    def __init__(
        self,
        wave: np.ndarray,
        *,
        lsf_sigma: np.ndarray | float = 0.5,
        n_spline_knots: int = 25,
        base_dir: str | Path | None = None,
        palace_suffix: str | None = None,
        palace_oh_suffix: str | None = None,
        palace_diffuse_suffix: str | None = None,
        o2_min_valid_frac: float = O2_MIN_VALID_FRAC,
        moon_smooth_lambda: float = 1e-3,
        moon_interline_boost: float = 0.0,
        moon_interline_red_min: float = 7454.0,
        moon_interline_exclusion_a: float = 3.0,
        moon_interline_line_flux_threshold: float = 0.0,
        split_zodi: bool = False,
        n_zodi_spline_knots: int = 3,
        zodi_smooth_lambda: float | None = None,
        moon_albedo_asset_path: str | pathlib.Path | None = None,
        moon_albedo_fiducial_phase_deg: float = 30.0,
        zodi_color_exponent: float = 0.26,
    ) -> None:
        self.wave = np.asarray(wave, float)
        self.lsf_sigma = np.asarray(lsf_sigma, float)
        self.n_spline_knots = int(n_spline_knots)
        self.o2_min_valid_frac = float(np.clip(o2_min_valid_frac, 0.0, 1.0))
        self.moon_smooth_lambda = max(float(moon_smooth_lambda), 0.0)
        self.moon_interline_boost = max(float(moon_interline_boost), 0.0)
        self.moon_interline_red_min = float(moon_interline_red_min)
        self.moon_interline_exclusion_a = max(float(moon_interline_exclusion_a), 0.0)
        self.moon_interline_line_flux_threshold = float(np.clip(moon_interline_line_flux_threshold, 0.0, 1.0))
        self.base_dir = (
            Path(base_dir).resolve() if base_dir is not None else Path(__file__).resolve().parent.parent
        )
        self.palace_suffix = self._validate_palace_suffix(palace_suffix)
        self.palace_oh_suffix = self._validate_palace_suffix(palace_oh_suffix)
        self.palace_diffuse_suffix = self._validate_palace_suffix(palace_diffuse_suffix)
        self.pmd_dir = self.base_dir / "palace" / "PMD"
        self.solar_path = self.base_dir / "Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt"

        self.t_o2 = T_O2_REF
        self.t_o2_err = np.nan
        self.o2_prefit_amp = np.nan
        self.vector_o2 = np.zeros_like(self.wave)
        self.matrix_o2 = self.vector_o2[None, :]
        self.vector_o2_stick = np.zeros_like(self.wave)
        self.matrix_o2_stick = self.vector_o2_stick[None, :]
        self.o2_prefit_bestfit = self.vector_o2.copy()
        self.bestfit = np.zeros_like(self.wave)
        self.coef = np.array([], float)
        self.coef_err = np.array([], float)
        self.fit_status = ""
        self.fit_summary = ""
        self.r2 = np.nan
        self.rms_resid = np.nan
        self.peak_memory_mb = np.nan
        self.o2_fit_status = ""
        self.o2_fit_summary = ""
        self.o2_fit_elapsed_sec = 0.0
        self.o2_valid_frac = np.nan
        self.lsf_kernels: dict[str, np.ndarray] = {}
        self.lsf_metrics: dict[str, dict[str, object]] = {}
        self.bestfit_lsf = np.zeros_like(self.wave)
        self.moon_knots_used = np.array([], float)
        self.moon_boosted_pixels_used = np.array([], float)

        # ---- Split-zodi configuration ----
        self.split_zodi = bool(split_zodi)
        self.n_zodi_spline_knots = int(n_zodi_spline_knots)
        if zodi_smooth_lambda is None:
            self.zodi_smooth_lambda = 10.0 * self.moon_smooth_lambda
        else:
            self.zodi_smooth_lambda = max(float(zodi_smooth_lambda), 0.0)
        self.moon_albedo_asset_path = (
            pathlib.Path(moon_albedo_asset_path)
            if moon_albedo_asset_path is not None
            else self.base_dir / 'moon_zodi'
              / 'eso_skycalc_rolo_moon_albedo.dat'
        )
        self.moon_albedo_fiducial_phase_deg = float(moon_albedo_fiducial_phase_deg)
        self.zodi_color_exponent = float(zodi_color_exponent)
        # Populated by _build_static_basis when split_zodi=True; otherwise empty.
        self.matrix_zodi = np.zeros((0, self.wave.size), dtype=float)
        self.matrix_zodi_hr = np.zeros((0, self.wave.size), dtype=float)
        self.zodi_names: list[str] = []
        self.zodi_knots_used = np.array([], float)
        self.moon_albedo_shape = np.ones_like(self.wave)
        self.zodi_color_shape = np.ones_like(self.wave)
        # Cached second-difference operators for the moon / zodi curvature penalties;
        # populated in _build_moon / _build_zodi.  Shape is (n_par - 2, n_par).
        self._d2_moon = np.zeros((0, 0), dtype=float)
        self._d2_zodi = np.zeros((0, 0), dtype=float)

        self._build_static_basis()

    def fit(
        self,
        flux: np.ndarray,
        ivar: np.ndarray,
        *,
        verbose: bool = False,
        n_lsf_refits: int = 1,
    ) -> SkyDecompResult:
        flux = np.asarray(flux, float)
        ivar = np.asarray(ivar, float)
        n_lsf_refits = max(int(n_lsf_refits), 0)
        trace_started = False
        if not tracemalloc.is_tracing():
            tracemalloc.start()
            trace_started = True
        tracemalloc.reset_peak()
        t0 = time.perf_counter()

        self._prefit_o2(flux, ivar)
        self.design_matrix = self._assemble_design_matrix()
        comp_slices = self._component_slices(
            self._matrix_bundle(
                self.matrix_oh,
                self.matrix_moon,
                self.matrix_diffuse,
                self.matrix_atom,
                self.matrix_orc,
                self.matrix_o2,
                matrix_zodi=self.matrix_zodi if self.split_zodi else None,
            )
        )
        first = self._fit_design(
            self.design_matrix,
            flux,
            ivar,
            moon_slice=comp_slices["moon"],
            diffuse_slice=comp_slices["diffuse"],
            zodi_slice=comp_slices.get("zodi"),
        )
        self.bestfit = first["bestfit"]
        refined = first
        final_matrices = self._matrix_bundle(
            self.matrix_oh,
            self.matrix_moon,
            self.matrix_diffuse,
            self.matrix_atom,
            self.matrix_orc,
            self.matrix_o2,
            matrix_zodi=self.matrix_zodi if self.split_zodi else None,
        )
        self.lsf_kernels = {}
        self.lsf_metrics = {}
        iter_logs: list[dict[str, object]] = []

        for i_refit in range(1, n_lsf_refits + 1):
            source_lsf, fixed_background = self._build_lsf_source(refined["coef"])
            self._fit_lsf_channels(flux, ivar, source_lsf, fixed_background)
            final_matrices = self._assemble_refined_matrices()
            _refined_names = ["oh", "moon"]
            if "zodi" in comp_slices:
                _refined_names.append("zodi")
            _refined_names.extend(["diffuse", "atom", "orc", "o2"])
            refined_design = np.vstack([final_matrices[name] for name in _refined_names])
            refined = self._fit_design(
                refined_design,
                flux,
                ivar,
                moon_slice=comp_slices["moon"],
                diffuse_slice=comp_slices["diffuse"],
                zodi_slice=comp_slices.get("zodi"),
            )
            iter_logs.append(
                {
                    "iter": i_refit,
                    "lsf": {
                        ch: {
                            "status": met["status"],
                            "chi2_red": met["chi2_red"],
                            "sigma_pix": met["sigma_pix"],
                            "center_pix": met["center_pix"],
                            "runtime_sec": met["runtime_sec"],
                        }
                        for ch, met in self.lsf_metrics.items()
                    },
                    "fit": {
                        "status": refined["status"],
                        "chi2_red": refined["reduced_chi2"],
                        "r2": refined["r2"],
                        "rms": refined["rms_resid"],
                        "qp_dt": refined["qp_elapsed_sec"],
                    },
                }
            )

        components = self._components_from_coef(refined["coef"], final_matrices)
        # `final_matrices['o2']` is a (1, n_wave) block that already incorporates
        # every LSF refit; keep `vector_o2` in step so what we persist matches
        # the O2 basis used to build `components['o2']`.
        if final_matrices["o2"].shape[0] == 1:
            self.vector_o2 = final_matrices["o2"][0].copy()

        # Per-pixel 1σ of the refined reconstruction (LSF-aware Jacobian
        # propagation).  Uses the SAME matrix bundle as `_components_from_coef`
        # above so the mean and σ come from an identical LSF/O2 basis.  Needs
        # the finalised coef_err from `refined`, not `self.coef_err` (which is
        # not filled until a few lines below).
        _sigma_coef_err = np.asarray(
            refined.get("coef_err", np.full(refined["coef"].shape, np.nan)),
            dtype=float,
        )
        _sigma_comps = self._components_sigma_from_coef_err(_sigma_coef_err, final_matrices)
        bestfit_lsf_sigma = np.sqrt(
            _sigma_comps["oh"] ** 2
            + _sigma_comps["moon"] ** 2
            + _sigma_comps["diffuse"] ** 2
            + _sigma_comps["atom"] ** 2
            + _sigma_comps["orc"] ** 2
            + _sigma_comps["o2"] ** 2
        )

        fit_elapsed_sec = time.perf_counter() - t0
        peak_memory_mb = tracemalloc.get_traced_memory()[1] / 1024**2
        if trace_started:
            tracemalloc.stop()
        self.peak_memory_mb = peak_memory_mb
        self.bestfit_lsf = refined["bestfit"]
        self.coef = refined["coef"]
        self.coef_err = np.asarray(
            refined.get("coef_err", np.full_like(self.coef, np.nan)),
            dtype=float,
        )
        self.fit_status = refined["status"]
        self.fit_summary = (
            f"status={refined['status']} | npar={refined['n_par']} | ngood={refined['n_good']} | "
            f"chi2_red={refined['reduced_chi2']:.4g} | R2={refined['r2']:.5f} | "
            f"qp_dt={refined['qp_elapsed_sec']:.2f}s | moon_smooth_lambda={self.moon_smooth_lambda:.3g} | "
            f"moon_mask_downweight_boost={self.moon_interline_boost:.3g} | "
            f"moon_interline_line_flux_threshold={self.moon_interline_line_flux_threshold:.3g} | "
            f"n_lsf_refits={n_lsf_refits} | dt={fit_elapsed_sec:.2f}s"
        )
        self.r2 = refined["r2"]
        self.rms_resid = refined["rms_resid"]

        if verbose:
            print("O2")
            print(f"  T         {self.t_o2:.2f} +/- {self.t_o2_err:.2f} K")
            print(f"  chi2_red  {self._fmt_num(self._extract_o2_chi2_red())}")
            print(f"  dt        {self.o2_fit_elapsed_sec:.3f} s")
            print()
            print("decomp")
            print(f"  init      {self._fmt_num(first['reduced_chi2'])}")
            print()
            print("iterations")
            for log in iter_logs:
                lsf_dt = 0.0
                vals = {}
                for ch in ("B", "R", "Z"):
                    met = log["lsf"].get(ch)
                    vals[ch] = self._fmt_num(np.nan if met is None else met["chi2_red"])
                    if met is not None:
                        lsf_dt += float(met["runtime_sec"])
                fit_log = log["fit"]
                print(
                    f"  [{log['iter']}] LSF     "
                    f"B={vals['B']}   R={vals['R']}   Z={vals['Z']}   dt={lsf_dt:.3f} s"
                )
                print(
                    f"      decomp  chi2_red={self._fmt_num(fit_log['chi2_red'])}   "
                    f"qp={fit_log['qp_dt']:.3f} s"
                )
            if not iter_logs:
                print("  none")
            print()
            print("final")
            print(
                f"  decomp    {self._fmt_num(refined['reduced_chi2'])}\n"
                f"  refits    {n_lsf_refits}\n"
                f"  total_dt  {fit_elapsed_sec:.3f} s\n"
                f"  peak_mem  {peak_memory_mb:.2f} MB"
            )

        return SkyDecompResult(
            coef=self.coef,
            coef_err=self.coef_err,
            bestfit=first["bestfit"],
            resid=first["resid"],
            resid_level=first["resid_level"],
            fit_status=self.fit_status,
            fit_summary=self.fit_summary,
            reduced_chi2=refined["reduced_chi2"],
            fit_elapsed_sec=fit_elapsed_sec,
            components=components,
            design_names=self.design_names,
            t_o2=self.t_o2,
            t_o2_err=self.t_o2_err,
            r2=self.r2,
            rms_resid=self.rms_resid,
            peak_memory_mb=peak_memory_mb,
            o2_fit_status=self.o2_fit_status,
            o2_fit_summary=self.o2_fit_summary,
            o2_fit_elapsed_sec=self.o2_fit_elapsed_sec,
            o2_valid_frac=self.o2_valid_frac,
            lsf_kernels=self.lsf_kernels,
            lsf_metrics=self.lsf_metrics,
            bestfit_lsf=self.bestfit_lsf,
            moon_knots=self.moon_knots_used.copy(),
            moon_boosted_pixels=self.moon_boosted_pixels_used.copy(),
            vector_o2=self.vector_o2.copy(),
            o2_prefit_amp=float(self.o2_prefit_amp),
            bestfit_lsf_sigma=bestfit_lsf_sigma,
            zodi_names=list(self.zodi_names),
            zodi_knots=self.zodi_knots_used.copy(),
            coef_cov_moon=refined.get("coef_cov_moon"),
            coef_cov_zodi=refined.get("coef_cov_zodi"),
        )

    @staticmethod
    def _matrix_bundle(
        matrix_oh: np.ndarray,
        matrix_moon: np.ndarray,
        matrix_diffuse: np.ndarray,
        matrix_atom: np.ndarray,
        matrix_orc: np.ndarray,
        matrix_o2: np.ndarray,
        matrix_zodi: np.ndarray | None = None,
    ) -> dict[str, np.ndarray]:
        out: dict[str, np.ndarray] = {
            "oh": matrix_oh,
            "moon": matrix_moon,
            "diffuse": matrix_diffuse,
            "atom": matrix_atom,
            "orc": matrix_orc,
            "o2": matrix_o2,
        }
        if matrix_zodi is not None and matrix_zodi.shape[0] > 0:
            out["zodi"] = matrix_zodi
        return out

    @staticmethod
    def _per_column_chi2_from_residuals(
        design_matrix: np.ndarray,
        resid: np.ndarray,
        ivar: np.ndarray,
        good: np.ndarray,
        *,
        fallback: float,
    ) -> np.ndarray:
        """Per-column pseudo-chi² from residual power weighted by column support.

        For each column ``j`` of the physical design matrix, compute a
        heteroscedasticity-consistent (HC1-style) chi² by weighting the
        pixel-wise weighted residual power by that column's spectral support:

            w_j[i]        = |A[i, j]| / max_i |A[i, j]|
            chi²_j        = Σ_i w_j[i] · ivar[i] · resid[i]²  /  Σ_i w_j[i]

        The ``max_i |A[i, j]|`` normalisation cancels in the ratio, so the
        implementation collapses to two BLAS gemv calls over ``|A[good, :]|``:

            chi²_j        = (Σ_i |A[i, j]| · ivar[i] · resid[i]²) / Σ_i |A[i, j]|

        Columns with vanishing support (all-zero design over the good pixels)
        fall back to ``fallback`` (typically the aggregate ``reduced_chi²``).
        This gives ``_coef_err_active_set`` a per-column inflation factor
        instead of the scalar ``max(reduced_chi², 1)`` used historically,
        matching the MLP-side observation (top methods §7) that the
        miscalibration is column-dependent.
        """
        # design_matrix here is (n_col, n_pix) -- see _fit_design.
        n_par = int(design_matrix.shape[0])
        per_col = np.full(n_par, float(fallback), dtype=float)
        if not good.any():
            return per_col
        # Take absolute value once on the good-pixel slice (n_col, n_good);
        # the transpose in the historical implementation is skipped because
        # gemv over rows is just as fast as gemv over columns and saves the
        # extra (n_good, n_col) copy.
        abs_slice = np.abs(design_matrix[:, good])  # (n_col, n_good)
        r2w = (resid[good] * resid[good]) * ivar[good]  # (n_good,)
        num = abs_slice @ r2w  # (n_col,) -- single gemv
        den = abs_slice.sum(axis=1)  # (n_col,)
        support_mask = den > 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            chi2_j = num / np.where(support_mask, den, 1.0)
        chi2_j = np.where(np.isfinite(chi2_j), chi2_j, float(fallback))
        per_col[support_mask] = chi2_j[support_mask]
        return per_col

    @staticmethod
    def _coef_err_active_set(
        coef: np.ndarray,
        source_p: list[np.ndarray],
        source_col_scale: list[np.ndarray],
        source_data_scale: list[float],
        source_local_index: list[int],
        reduced_chi2: float,
        per_column_chi2: np.ndarray | None = None,
        cov_block_slices: dict[str, slice] | None = None,
    ) -> np.ndarray | tuple[np.ndarray, dict[str, np.ndarray]]:
        """Active-set posterior 1σ uncertainty for a nonnegative WLS solve.

        Notes
        -----
        The QP minimises

            J(c) = ½‖W(A c − y)‖² + ½ λ ‖L c‖²    subject to    c ≥ 0

        where ``W = diag(ivar_good)`` is the weighting used by ``_fit_design``,
        ``A`` is the physical design matrix and ``L`` is the moon-spline
        second-difference operator (with strength ``λ = 2·moon_smooth_lambda``).
        For readability the internal solve additionally divides the data by a
        scalar ``data_scale`` and rescales columns by ``col_scale``; those
        transformations cancel out of the KKT solution but must be undone here
        before the covariance is reported in native coefficient units.

        At the KKT solution, split coefficients into

            active   := {j : c_j > tol}      inactive := {j : c_j = 0}

        Holding the inactive set at zero collapses the constrained problem to
        an unconstrained WLS on the active columns

            c_a ← argmin ‖W (A_a c_a − y)‖² + λ ‖L_a c_a‖²

        whose Hessian is the active sub-block of the (regularised) Fisher
        information ``H = AᵀWA + λ LᵀL``. The corresponding maximum a
        posteriori covariance is

            Cov(c_a) = H_{a,a}⁻¹                                    (1)

        and it is what the routine reports.  Equivalently, this is the
        classical Cramér–Rao unconstrained covariance of the sub-model that
        contains only the active coefficients.

        Smoothness prior. ``H`` above is the *regularised* Fisher information,
        i.e. it already includes the ``λ LᵀL`` term added at solve time
        (`` _fit_design`` line "p_dense_local[..] += 2·smooth_λ·D²ᵀD² ").  So
        the σ reported here reflects both the data likelihood and the
        smoothness prior on adjacent B-spline coefficients; no separate
        Bayesian tightening step is needed.

        Undoing the internal scaling. What is stored in the QP is the dense
        matrix

            P_scaled = diag(1/col_scale) · (H / data_scale²) · diag(1/col_scale)
                       + regulariser

        so ``H⁻¹ = diag(1/col_scale) · P_scaled⁻¹ · diag(1/col_scale) · data_scale²``.
        The active sub-block of ``P_scaled`` is inverted once per solve group
        (Hessian identity is the group key so a two-stage moon/diffuse refit
        is handled cleanly), the column-scale factors are folded in via
        ``inv_col = 1/col_scale[a]``, and the diagonal gives ``Var(c_j)`` in
        physical units. The 1σ error is ``sqrt(Var)``.

        Noise inflation. Empirical ``reduced_chi²`` above unity indicates that
        the input ``ivar`` under-estimates the noise (or the model is
        mildly under-specified); the variance is therefore inflated by
        ``max(reduced_chi², 1)`` before the square root. This is the standard
        "sigma-hat" adjustment applied to weighted least squares.

        Per-column inflation (2026-08-16). The MLP sigma-calibration diagnostic
        (top methods §7) showed the aggregate ``reduced_chi²`` scalar
        systematically over-inflates well-fit columns and under-inflates poorly
        fit ones. When ``per_column_chi2`` is provided (see ``_fit_design``),
        the inflation factor becomes column-specific ``max(per_column_chi2[j], 1)``
        computed from residual power weighted by each column's spectral support.
        This is a lightweight heteroscedasticity-consistent (HC1-style) fix that
        preserves the active-set covariance structure and matches the physical
        intuition that emission-line columns and continuum-spline columns have
        very different residual behaviour. If ``per_column_chi2`` is ``None`` the
        routine falls back to the scalar ``reduced_chi²`` for backward
        compatibility.

        Boundary treatment. Inactive coefficients are pinned at c_j = 0 by
        the constraint and have no symmetric 1σ interval defined; they are
        reported as ``NaN``. A one-sided upper limit is easy to derive from
        ``|A_jᵀ W (y − A c)| / √(A_jᵀ W A_j)`` but is not written here to
        keep the array shape aligned with ``coef`` and avoid mixing
        interval-type semantics on a per-row basis.

        Caveats. Active-set covariance ignores the constraint's effect on
        coefficients that are close to zero: a small perturbation of the data
        can flip a boundary coefficient into the active set. A bootstrap or
        profile-likelihood analysis is more faithful in that regime but not
        needed for the calibration MLP downstream.

        Joint covariance blocks (2026-08-25). When ``cov_block_slices`` is
        supplied -- a mapping ``{name: slice}`` selecting one or more
        contiguous coefficient ranges (typically ``{'moon': moon_slice,
        'zodi': zodi_slice}``) -- the routine also assembles the full active-
        set covariance restricted to each named block and returns it as a 2D
        matrix (n_block × n_block) alongside the usual 1D ``coef_err``.  Off-
        diagonal correlations are computed with the same column-scale-undo
        and per-column-χ² inflation logic as the diagonal (off-diagonal
        inflation uses ``sqrt(inflate[i] * inflate[j])`` so the correlation
        matrix is preserved and the diagonal reproduces the vector ``coef_err``).
        Inactive-coefficient rows and columns are ``NaN``.  Return type is
        ``(coef_err, {name: cov_2d})`` in this case; when
        ``cov_block_slices`` is None (default) the routine returns just
        ``coef_err`` for backwards compatibility.
        """
        n_par = int(coef.size)
        coef_err = np.full(n_par, np.nan, dtype=float)
        want_blocks = cov_block_slices is not None
        if want_blocks:
            cov_blocks = {
                name: np.full(
                    (sl.stop - sl.start, sl.stop - sl.start), np.nan, dtype=float)
                for name, sl in cov_block_slices.items()
            }
        if n_par == 0:
            if want_blocks:
                return coef_err, cov_blocks
            return coef_err

        active_tol = 1e-8 * float(max(np.max(coef), 1.0))
        active_mask = coef > active_tol
        coef_err[~active_mask] = np.nan

        scalar_inflate = float(max(reduced_chi2, 1.0))
        if per_column_chi2 is not None:
            per_col = np.asarray(per_column_chi2, dtype=float)
            if per_col.shape != coef.shape:
                raise ValueError(
                    f"per_column_chi2 shape {per_col.shape} != coef shape {coef.shape}"
                )
            per_col = np.where(np.isfinite(per_col) & (per_col > 0.0), per_col,
                                scalar_inflate)
            inflate_arr = np.maximum(per_col, 1.0)
        else:
            inflate_arr = np.full(coef.shape, scalar_inflate, dtype=float)

        # Group active coefficients by the solve they came from (identity of the
        # Hessian array); within each group we invert one small sub-block once.
        groups: dict[int, list[int]] = {}
        for j in np.flatnonzero(active_mask):
            groups.setdefault(id(source_p[int(j)]), []).append(int(j))

        for solve_indices in groups.values():
            first = solve_indices[0]
            p_solve = source_p[first]
            cs_solve = source_col_scale[first]
            ds_solve = source_data_scale[first]
            local_idx = np.asarray(
                [source_local_index[j] for j in solve_indices],
                dtype=int,
            )
            H_a = p_solve[np.ix_(local_idx, local_idx)]
            H_a = 0.5 * (H_a + H_a.T) + 1e-12 * np.eye(H_a.shape[0])
            try:
                cov_a = np.linalg.inv(H_a)
            except np.linalg.LinAlgError:
                cov_a = np.linalg.pinv(H_a)
            inv_col = 1.0 / cs_solve[local_idx]
            # Undo column and data scaling; inflate by per-column χ² (falling
            # back to the scalar reduced_chi² when per_column_chi2 is None); clip
            # to nonnegative before the square root to guard against tiny
            # negative eigenvalues from finite-precision arithmetic.
            #
            # diag(D · cov_a · D) = d² · diag(cov_a) when D is diagonal, so the
            # historical outer-product materialisation is skipped.
            solve_indices_arr = np.asarray(solve_indices, dtype=int)
            var_a = np.clip(
                (inv_col * inv_col)
                * np.diagonal(cov_a)
                * (ds_solve ** 2)
                * inflate_arr[solve_indices_arr],
                0.0,
                np.inf,
            )
            coef_err[solve_indices_arr] = np.sqrt(var_a)

            if want_blocks:
                # Materialise the physical-unit covariance for this solve group
                # (small square matrix; only done here in the block path).  This
                # applies the same column-scale-undo and χ² inflation as the
                # diagonal, so the sqrt of its diagonal reproduces coef_err.
                inflate_row = inflate_arr[solve_indices_arr]
                scale_row = inv_col * ds_solve * np.sqrt(inflate_row)
                cov_a_phys = (scale_row[:, None] * scale_row[None, :]) * cov_a
                for name, sl in cov_block_slices.items():
                    block_mask = (
                        (solve_indices_arr >= sl.start)
                        & (solve_indices_arr < sl.stop)
                    )
                    if not block_mask.any():
                        continue
                    block_local = solve_indices_arr[block_mask] - sl.start
                    block_source = np.flatnonzero(block_mask)
                    sub = cov_a_phys[np.ix_(block_source, block_source)]
                    cov_blocks[name][np.ix_(block_local, block_local)] = sub

        if want_blocks:
            return coef_err, cov_blocks
        return coef_err

    def _fit_design(
        self,
        design_matrix: np.ndarray,
        flux: np.ndarray,
        ivar: np.ndarray,
        *,
        moon_slice: slice | None = None,
        diffuse_slice: slice | None = None,
        zodi_slice: slice | None = None,
    ) -> dict[str, object]:
        """Weighted nonnegative least-squares fit of one design matrix.

        Solves the quadratic program

            c* = argmin_{c ≥ 0} ½‖diag(√ivar) (A c − y)‖²
                                + λ ‖L c‖²    (moon-spline curvature)

        on the finite/positive-ivar subset of the data. ``A = design_matrix``
        holds the physical basis; ``λ = moon_smooth_lambda`` (from the model
        config) is nonzero only for the moon block, whose curvature operator
        ``L`` penalises second differences of adjacent spline coefficients.

        For numerical stability the solve is performed on a doubly-rescaled
        version of the problem: ``y'' = y / data_scale`` and
        ``A'' = A / (data_scale · col_scale)`` where ``data_scale`` normalises
        the target RMS to unity and ``col_scale`` normalises each column to
        unit ℓ²-norm. The nonnegative QP is delegated to Clarabel. The
        returned coefficient is the fitted vector rescaled back into physical
        units, ``c = x* / col_scale``.

        An optional second solve refits the moon+diffuse (``target_cols``)
        block with an interline-boosted weight vector while the remaining
        coefficients are held at their first-solve values (subtracted from
        the data). When this happens each moon/diffuse coefficient's covariance
        is taken from the target-refit Hessian, not the first-solve Hessian —
        see the per-column provenance arrays ``source_p``, ``source_col_scale``,
        ``source_data_scale`` and ``source_local_index`` prepared here.

        Uncertainties
        -------------
        ``coef_err`` is filled by ``_coef_err_active_set`` right before the
        result dict is built. See that method for the full derivation; briefly,
        it inverts the active sub-block of the regularised Fisher information
        stored on ``source_p`` (one small linear-algebra call per solve group),
        undoes the ``col_scale`` / ``data_scale`` transformations, inflates
        variance by ``max(reduced_chi², 1)`` for a hat-sigma adjustment and
        reports the square root. Coefficients pinned at c = 0 by the c ≥ 0
        constraint receive ``NaN``.

        Returns
        -------
        dict with ``coef``, ``coef_err``, ``status``, ``bestfit``, ``resid``,
        ``resid_level``, ``n_good``, ``n_par``, ``chi2``, ``reduced_chi2``,
        ``r2``, ``rms_resid``, ``qp_elapsed_sec``.
        """
        good = np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0)
        y = flux[good]
        base_w = np.sqrt(ivar[good])
        self.moon_boosted_pixels_used = np.array([], float)
        boosted_w = None
        if self.moon_interline_boost > 0.0:
            interline_weight = self._moon_interline_weights()
            masked_mask = good & (interline_weight < (1.0 - 1e-12))
            self.moon_boosted_pixels_used = self.wave[masked_mask]
            boosted_w = base_w * np.sqrt(interline_weight[good])
        n_good = int(np.sum(good))

        a = design_matrix[:, good].T

        def _solve_nonnegative_weighted(
            a_mat: np.ndarray,
            y_vec: np.ndarray,
            w_vec: np.ndarray,
            moon_slice_local: slice | None,
            zodi_slice_local: slice | None = None,
        ) -> tuple[np.ndarray, str, float, np.ndarray, np.ndarray, float]:
            aw = a_mat * w_vec[:, None]
            yw = y_vec * w_vec
            data_scale = max(float(np.sqrt(np.nanmean(yw**2))) if yw.size else 1.0, 1.0)
            aw = aw / data_scale
            yw = yw / data_scale
            col_scale = np.sqrt(np.sum(aw**2, axis=0))
            col_scale = np.where(np.isfinite(col_scale) & (col_scale > 0), col_scale, 1.0)
            aw = aw / col_scale[None, :]

            n_par_local = int(aw.shape[1])
            p_dense_local = aw.T @ aw

            # Penalize curvature of moon-spline coefficients to suppress oscillations.
            if self.moon_smooth_lambda > 0.0 and moon_slice_local is not None:
                i0 = moon_slice_local.start or 0
                i1 = moon_slice_local.stop or i0
                n_moon = i1 - i0
                # Use cached D2 when it matches; else build on the fly (interline refit
                # solves a submatrix whose moon sub-block may not equal the cached size).
                d2 = self._d2_moon if self._d2_moon.shape == (max(n_moon - 2, 0), n_moon) else _build_d2_operator(n_moon)
                if d2.shape[0] > 0:
                    d2_scaled = d2 / col_scale[i0:i1][None, :]
                    p_dense_local[i0:i1, i0:i1] += 2.0 * self.moon_smooth_lambda * (d2_scaled.T @ d2_scaled)

            # Zodi-spline curvature penalty (same construction; usually a heavier lambda).
            if getattr(self, 'zodi_smooth_lambda', 0.0) > 0.0 and zodi_slice_local is not None:
                i0z = zodi_slice_local.start or 0
                i1z = zodi_slice_local.stop or i0z
                n_zodi = i1z - i0z
                d2z = self._d2_zodi if self._d2_zodi.shape == (max(n_zodi - 2, 0), n_zodi) else _build_d2_operator(n_zodi)
                if d2z.shape[0] > 0:
                    d2z_scaled = d2z / col_scale[i0z:i1z][None, :]
                    p_dense_local[i0z:i1z, i0z:i1z] += 2.0 * self.zodi_smooth_lambda * (d2z_scaled.T @ d2z_scaled)

            # Amplitude priors were removed 2026-08-26; ecliptic split alone carries identifiability.
            q_local = -(aw.T @ yw)
            p_local = sp.csc_matrix((p_dense_local + p_dense_local.T) / 2.0)
            p_local = sp.triu(p_local).tocsc()
            a_con = -sp.eye(n_par_local, format="csc")
            b_con = np.zeros(n_par_local, dtype=np.float64)
            cones = [clarabel.NonnegativeConeT(n_par_local)]
            settings = clarabel.DefaultSettings()
            settings.verbose = False

            t_qp_local = time.perf_counter()
            solver = clarabel.DefaultSolver(p_local, np.asarray(q_local, dtype=np.float64), a_con, b_con, cones, settings)
            qp_result_local = solver.solve()
            qp_dt_local = time.perf_counter() - t_qp_local
            coef_local = np.asarray(qp_result_local.x, float) / col_scale
            # `p_dense_local` is the column- and data-scaled Fisher information
            # (plus any regulariser); the caller needs it, `col_scale`, and
            # `data_scale` to back out the native-coordinate covariance.
            return (
                coef_local,
                str(qp_result_local.status),
                qp_dt_local,
                p_dense_local,
                col_scale,
                float(data_scale),
            )

        n_par = int(a.shape[1])
        (
            coef,
            status,
            qp_elapsed_sec,
            p_main,
            col_scale_main,
            data_scale_main,
        ) = _solve_nonnegative_weighted(
            a, y, base_w, moon_slice,
            zodi_slice_local=zodi_slice,
        )

        # Track whichever solve produced the FINAL value of each coefficient,
        # so uncertainties come from that solve's Hessian.
        source_p: list[np.ndarray] = [p_main] * n_par
        source_col_scale: list[np.ndarray] = [col_scale_main] * n_par
        source_data_scale: list[float] = [data_scale_main] * n_par
        source_local_index: list[int] = list(range(n_par))

        target_cols = []
        for comp_slice in (moon_slice, diffuse_slice):
            if comp_slice is None:
                continue
            i0 = comp_slice.start or 0
            i1 = comp_slice.stop or i0
            if i1 > i0:
                target_cols.extend(range(i0, i1))
        target_cols = sorted(set(target_cols))

        if boosted_w is not None and target_cols:
            target_cols_arr = np.asarray(target_cols, dtype=int)
            fixed_cols_arr = np.setdiff1d(np.arange(a.shape[1], dtype=int), target_cols_arr, assume_unique=True)
            y_target = y.copy()
            if fixed_cols_arr.size > 0:
                y_target = y_target - a[:, fixed_cols_arr] @ coef[fixed_cols_arr]

            a_target = a[:, target_cols_arr]
            moon_slice_local = None
            if moon_slice is not None:
                m0 = moon_slice.start or 0
                m1 = moon_slice.stop or m0
                if m1 > m0:
                    moon_global = np.arange(m0, m1, dtype=int)
                    local_pos = np.searchsorted(target_cols_arr, moon_global)
                    in_bounds = (local_pos >= 0) & (local_pos < target_cols_arr.size)
                    local_pos = local_pos[in_bounds]
                    local_pos = local_pos[target_cols_arr[local_pos] == moon_global[in_bounds]]
                    if local_pos.size > 0:
                        moon_slice_local = slice(int(local_pos.min()), int(local_pos.max()) + 1)
            coef_target, status_target, qp_dt_target, p_target, col_scale_target, data_scale_target = _solve_nonnegative_weighted(
                a_target,
                y_target,
                boosted_w,
                moon_slice_local,
                zodi_slice_local=None,
            )
            coef[target_cols_arr] = coef_target
            for local_pos_i, global_col in enumerate(target_cols_arr):
                source_p[int(global_col)] = p_target
                source_col_scale[int(global_col)] = col_scale_target
                source_data_scale[int(global_col)] = data_scale_target
                source_local_index[int(global_col)] = int(local_pos_i)
            status = f"{status} | md_refit={status_target}"
            qp_elapsed_sec += qp_dt_target

        bestfit = design_matrix.T @ coef
        resid = flux - bestfit
        resid_level = -3.0 * np.nanstd(resid)
        chi2 = float(np.sum(resid[good] ** 2 * ivar[good]))
        dof = max(n_good - n_par, 1)
        reduced_chi2 = chi2 / dof
        rms_resid = float(np.sqrt(np.nanmean(resid[good] ** 2))) if n_good else np.nan
        y_mean = float(np.average(y, weights=ivar[good])) if n_good else np.nan
        sst = float(np.sum((y - y_mean) ** 2 * ivar[good])) if n_good else np.nan
        r2 = 1.0 - chi2 / sst if np.isfinite(sst) and sst > 0 else np.nan
        per_column_chi2 = self._per_column_chi2_from_residuals(
            design_matrix, resid, ivar, good, fallback=reduced_chi2,
        )
        cov_block_slices: dict[str, slice] = {}
        if moon_slice is not None and (moon_slice.stop - moon_slice.start) > 0:
            cov_block_slices["moon"] = moon_slice
        if zodi_slice is not None and (zodi_slice.stop - zodi_slice.start) > 0:
            cov_block_slices["zodi"] = zodi_slice
        _err_ret = self._coef_err_active_set(
            coef,
            source_p,
            source_col_scale,
            source_data_scale,
            source_local_index,
            reduced_chi2,
            per_column_chi2=per_column_chi2,
            cov_block_slices=cov_block_slices or None,
        )
        if cov_block_slices:
            coef_err, cov_blocks = _err_ret
        else:
            coef_err = _err_ret
            cov_blocks = {}
        return {
            "coef": coef,
            "coef_err": coef_err,
            "coef_cov_moon": cov_blocks.get("moon"),
            "coef_cov_zodi": cov_blocks.get("zodi"),
            "status": status,
            "bestfit": bestfit,
            "resid": resid,
            "resid_level": resid_level,
            "n_good": n_good,
            "n_par": n_par,
            "chi2": chi2,
            "reduced_chi2": reduced_chi2,
            "r2": r2,
            "rms_resid": rms_resid,
            "qp_elapsed_sec": qp_elapsed_sec,
        }

    @staticmethod
    def _component_slices(mats: dict[str, np.ndarray]) -> dict[str, slice]:
        i0 = 0
        i1 = i0 + mats["oh"].shape[0]
        i2 = i1 + mats["moon"].shape[0]
        has_zodi = ("zodi" in mats) and (mats["zodi"].shape[0] > 0)
        i2z = i2 + (mats["zodi"].shape[0] if has_zodi else 0)
        i3 = i2z + mats["diffuse"].shape[0]
        i4 = i3 + mats["atom"].shape[0]
        i5 = i4 + mats["orc"].shape[0]
        i6 = i5 + mats["o2"].shape[0]
        out = {
            "oh": slice(i0, i1),
            "moon": slice(i1, i2),
            "diffuse": slice(i2z, i3),
            "atom": slice(i3, i4),
            "orc": slice(i4, i5),
            "o2": slice(i5, i6),
        }
        if has_zodi:
            out["zodi"] = slice(i2, i2z)
        return out

    def _components_from_coef(self, coef: np.ndarray, mats: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        sl = self._component_slices(mats)
        diffuse_coef = coef[sl["diffuse"]]
        components = {
            "oh": mats["oh"].T @ coef[sl["oh"]],
            "moon": mats["moon"].T @ coef[sl["moon"]],
            "ho2": diffuse_coef[0] * self.vector_ho2,
            "feo": diffuse_coef[1] * self.vector_feo,
            "o2ac": diffuse_coef[2] * self.vector_o2ac,
            "atom": mats["atom"].T @ coef[sl["atom"]],
            "orc": mats["orc"].T @ coef[sl["orc"]],
            "o2": mats["o2"].T @ coef[sl["o2"]],
        }
        components["diffuse"] = components["ho2"] + components["feo"] + components["o2ac"]
        if "zodi" in sl:
            components["zodi"] = mats["zodi"].T @ coef[sl["zodi"]]
        else:
            components["zodi"] = np.zeros_like(components["moon"])
        return components

    def _components_sigma_from_coef_err(
        self,
        coef_err: np.ndarray,
        mats: dict[str, np.ndarray],
    ) -> dict[str, np.ndarray]:
        """Propagate per-coefficient 1σ to per-component flux 1σ.

        For every matrix component (``oh``, ``moon``, ``atom``, ``orc``,
        ``o2``) the flux variance at wavelength λ is

            σ²_comp(λ) = Σ_j M[j, λ]² · σ_c[j]²                    (1)

        where ``M`` is the component's design matrix (already
        LSF-convolved by ``_assemble_refined_matrices``) and ``σ_c`` is the
        native-space per-coefficient uncertainty. This is the exact
        first-order (Jacobian) propagation because the reconstruction
        ``flux = M.T @ coef`` is linear in the coefficients.

        Implementation note. The sum is evaluated as ``(M*M).T @ σ²_c`` which
        dispatches to BLAS gemv on the squared matrix. Micro-benchmarks show
        this is faster than ``np.einsum('ij,ij,i->j', M, M, σ²_c)`` for
        typical shapes (n_coef ~ 400, n_pix ~ 3500) because einsum falls back
        to a general-purpose Python-side path for 3-tensor contractions
        instead of the tuned BLAS kernel.

        The diffuse sub-component uncertainties use the rows of the supplied
        diffuse design matrix. This keeps uncertainty propagation tied to the
        same final LSF-aware matrix bundle as the reconstructed mean and also
        supports reduced synthetic matrices used by workflow tests.

        NaN or non-finite ``coef_err`` values are treated as zero
        contribution, matching the "no measured uncertainty available"
        semantics used elsewhere in the pipeline.
        """
        sl = self._component_slices(mats)
        err = np.asarray(coef_err, dtype=np.float64).ravel()
        if err.size != sum(m.shape[0] for m in mats.values()):
            raise ValueError(
                f"coef_err length mismatch: expected "
                f"{sum(m.shape[0] for m in mats.values())}, got {err.size}"
            )
        err2 = np.where(np.isfinite(err), err ** 2, 0.0)

        def _matrix_sigma(name: str) -> np.ndarray:
            m = mats[name]
            # (m*m).T @ vec dispatches to BLAS gemv on the squared matrix.
            var = (m * m).T @ err2[sl[name]]
            return np.sqrt(np.maximum(var, 0.0))

        diffuse_err2 = err2[sl["diffuse"]]
        diffuse_matrix = mats["diffuse"]

        def _diffuse_row_sigma(index: int) -> np.ndarray:
            if index >= diffuse_matrix.shape[0]:
                return np.zeros(diffuse_matrix.shape[1], dtype=np.float64)
            return np.sqrt(diffuse_err2[index]) * np.abs(diffuse_matrix[index])

        sigmas = {
            "oh": _matrix_sigma("oh"),
            "moon": _matrix_sigma("moon"),
            "ho2": _diffuse_row_sigma(0),
            "feo": _diffuse_row_sigma(1),
            "o2ac": _diffuse_row_sigma(2),
            "diffuse": _matrix_sigma("diffuse"),
            "atom": _matrix_sigma("atom"),
            "orc": _matrix_sigma("orc"),
            "o2": _matrix_sigma("o2"),
        }
        sigmas["zodi"] = (
            _matrix_sigma("zodi") if "zodi" in sl
            else np.zeros_like(sigmas["moon"])
        )
        return sigmas

    def _build_lsf_source(self, coef: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        zodi_bundle = self.matrix_zodi_hr if (self.split_zodi and self.matrix_zodi.shape[0] > 0) else None
        mats = self._matrix_bundle(
            self.matrix_oh_stick,
            self.matrix_moon_hr,
            self.matrix_diffuse,
            self.matrix_atom_stick,
            self.matrix_orc_stick,
            self.matrix_o2_stick,
            matrix_zodi=zodi_bundle,
        )
        sl = self._component_slices(mats)
        emission = (
            mats["oh"].T @ coef[sl["oh"]]
            + mats["atom"].T @ coef[sl["atom"]]
            + mats["orc"].T @ coef[sl["orc"]]
            + mats["o2"].T @ coef[sl["o2"]]
        )
        moon = mats["moon"].T @ coef[sl["moon"]]
        if "zodi" in sl:
            moon = moon + mats["zodi"].T @ coef[sl["zodi"]]
        fixed_background = mats["diffuse"].T @ coef[sl["diffuse"]]
        return emission + moon, fixed_background

    def _assemble_refined_matrices(self) -> dict[str, np.ndarray]:
        zodi_hr = None
        if self.split_zodi and self.matrix_zodi.shape[0] > 0:
            zodi_hr = self._convolve_matrix_channelwise(self.matrix_zodi_hr)
        return self._matrix_bundle(
            self._convolve_matrix_channelwise(self.matrix_oh_stick),
            self._convolve_matrix_channelwise(self.matrix_moon_hr),
            self.matrix_diffuse,
            self._convolve_matrix_channelwise(self.matrix_atom_stick),
            self._convolve_matrix_channelwise(self.matrix_orc_stick),
            self._convolve_matrix_channelwise(self.matrix_o2_stick),
            matrix_zodi=zodi_hr,
        )

    def _build_static_basis(self) -> None:
        self.matrix_oh = self._build_oh()
        self.matrix_atom, self.atom_names = self._build_atom()
        self.matrix_orc, self.orc_names = self._build_orc()
        self.o2_names = ["O2_b01"]
        self._load_o2_model()
        self.matrix_moon, self.moon_names = self._build_moon()
        self.matrix_zodi, self.zodi_names = self._build_zodi()
        self.matrix_diffuse, self.diffuse_names = self._build_diffuse()
        self.design_names = (
            [f"OH_{i:03d}" for i in range(self.matrix_oh.shape[0])]
            + self.moon_names
            + self.zodi_names
            + self.diffuse_names
            + self.atom_names
            + self.orc_names
            + self.o2_names
        )
        self.design_matrix = self._assemble_design_matrix()

    def _assemble_design_matrix(self) -> np.ndarray:
        parts = [self.matrix_oh, self.matrix_moon]
        if self.split_zodi and self.matrix_zodi.shape[0] > 0:
            parts.append(self.matrix_zodi)
        parts.extend([self.matrix_diffuse, self.matrix_atom, self.matrix_orc, self.matrix_o2])
        return np.vstack(parts)

    def _build_oh(self) -> np.ndarray:
        oh = Table.read(self._pmd_path("pmd_popmodel_OH.dat"), format="ascii.basic", guess=False, comment="#", fast_reader=False)
        oh["wave"] = vac_to_air(np.asarray(oh["lam"], float) * 1e4)
        mask = (oh["wave"] >= self.wave.min() - CAP_WAVE) & (oh["wave"] <= self.wave.max() + CAP_WAVE)
        oh = decode_hitran_id(oh[mask])
        groups = oh.group_by(OH_GROUP_KEYS).groups

        matrix = np.zeros((len(groups), self.wave.size))
        self.matrix_oh_stick = np.zeros_like(matrix)
        for idx, grp in enumerate(groups):
            amp = np.asarray(grp["Aij"] * grp["gi"], float)
            matrix[idx] = grp2vector(grp["wave"], amp, self.wave, self.lsf_sigma)
            self.matrix_oh_stick[idx] = sticks2vector(grp["wave"], amp, self.wave)
        return matrix

    def _build_moon(self) -> tuple[np.ndarray, list[str]]:
        sol = np.loadtxt(self._require_path(self.solar_path), comments=";")
        solar_wave = vac_to_air(sol[:, 0] * 10.0)
        solar_flux = np.asarray(sol[:, 1], float)
        solar_flux /= np.nanmedian(solar_flux)
        solar_hr = np.nan_to_num(np.interp(self.wave, solar_wave, solar_flux, left=0.0, right=0.0))

        solar_rb = rebin_and_convolve(
            self.wave,
            solar_wave,
            solar_flux,
            self.lsf_sigma * 2.355,
            lsf_in_wavelength=True,
        )
        solar_rb /= np.nanmedian(solar_rb)
        self.vector_moon = solar_rb

        w0, w1 = self.wave[0], self.wave[-1]
        interior = self._uniform_moon_knots(w0, w1, self.n_spline_knots)

        self.moon_knots_used = interior.copy()

        t_knots = np.r_[(w0,) * 4, interior, (w1,) * 4]
        matrix_bspl = BSpline.design_matrix(self.wave, t_knots, 3).toarray()
        # Colour envelope: solar × moon_albedo when split_zodi, else solar alone.
        if getattr(self, 'split_zodi', False):
            self.moon_albedo_shape = self._load_moon_albedo_shape()
            envelope = solar_rb * self.moon_albedo_shape
            envelope_hr = solar_hr * self.moon_albedo_shape
        else:
            envelope = solar_rb
            envelope_hr = solar_hr
        matrix_moon = (envelope[:, None] * matrix_bspl).T
        self.matrix_moon_hr = (envelope_hr[:, None] * matrix_bspl).T
        moon_names = [f"Moon_bs{i:02d}" for i in range(matrix_moon.shape[0])]
        # Precompute the moon curvature penalty operator once (depends only on n_par).
        self._d2_moon = _build_d2_operator(matrix_moon.shape[0])
        return matrix_moon, moon_names

    def _load_moon_albedo_shape(self) -> np.ndarray:
        """Interpolate ROLO lunar albedo at the fiducial phase onto the fit grid."""
        path = self.moon_albedo_asset_path
        with pathlib.Path(path).open('r', encoding='utf-8') as fh:
            lines = [ln.rstrip() for ln in fh if ln.strip() and not ln.startswith('#')]
        constants = np.fromstring(lines[0], sep=' ')
        n_lam = int(lines[1])
        coefs = np.array([np.fromstring(ln, sep=' ') for ln in lines[2:2 + n_lam]])
        phase_deg_abs = abs(self.moon_albedo_fiducial_phase_deg)
        phase_rad = np.deg2rad(phase_deg_abs)
        signed_rad = np.deg2rad(self.moon_albedo_fiducial_phase_deg)
        signed_lim = (signed_rad if phase_deg_abs < 97.0
                      else 97.0 * signed_rad / phase_deg_abs)
        v = coefs[:, 1:]
        poly = (v[:, 0] + v[:, 1] * phase_rad + v[:, 2] * phase_rad**2 + v[:, 3] * phase_rad**3
                + v[:, 4] * signed_lim + v[:, 5] * signed_lim**3 + v[:, 6] * signed_lim**5)
        opp = (v[:, 7] * np.exp(-phase_deg_abs / constants[0])
               + v[:, 8] * np.exp(-phase_deg_abs / constants[1])
               + v[:, 9] * np.cos((phase_deg_abs - constants[2]) / constants[3]))
        tabulated = np.exp(poly + opp) / 0.87
        wave_nm = self.wave / 10.0
        albedo = np.interp(wave_nm, coefs[:, 0], tabulated)
        # Normalise to unit median so Moon_bs coefficient magnitudes
        # stay comparable to the current (no-albedo) fit-time scale.
        med = float(np.nanmedian(albedo))
        return (albedo / med).astype(float) if np.isfinite(med) and med > 0 else albedo

    def _zodi_color_shape_vec(self) -> np.ndarray:
        """Parametric zodi reddening: (λ/5000Å)^exponent, normalised to unit median."""
        shape = (self.wave / 5000.0) ** self.zodi_color_exponent
        med = float(np.nanmedian(shape))
        return (shape / med).astype(float) if np.isfinite(med) and med > 0 else shape

    def _build_zodi(self) -> tuple[np.ndarray, list[str]]:
        """Build Zodi_bs family: solar_rb × zodi_color × B-spline(K_zodi).

        Called only when `split_zodi=True`.  Uses the same solar template as the moon
        family (already computed in `_build_moon`), so must be invoked AFTER `_build_moon`.
        """
        if not self.split_zodi:
            return np.zeros((0, self.wave.size), dtype=float), []
        # Solar template is stored on self.vector_moon (populated in _build_moon).
        solar_rb = self.vector_moon
        self.zodi_color_shape = self._zodi_color_shape_vec()
        w0, w1 = self.wave[0], self.wave[-1]
        interior = self._uniform_moon_knots(w0, w1, self.n_zodi_spline_knots)
        self.zodi_knots_used = interior.copy()
        t_knots = np.r_[(w0,) * 4, interior, (w1,) * 4]
        matrix_bspl = BSpline.design_matrix(self.wave, t_knots, 3).toarray()
        envelope = solar_rb * self.zodi_color_shape
        matrix_zodi = (envelope[:, None] * matrix_bspl).T
        # HR twin: identical for now; LSF refit convolves it channelwise if used.
        self.matrix_zodi_hr = matrix_zodi.copy()
        zodi_names = [f'Zodi_bs{i:02d}' for i in range(matrix_zodi.shape[0])]
        # Precompute the zodi curvature penalty operator (depends only on n_par).
        self._d2_zodi = _build_d2_operator(matrix_zodi.shape[0])
        return matrix_zodi, zodi_names

    @staticmethod
    def _uniform_moon_knots(w0: float, w1: float, n_knots: int) -> np.ndarray:
        n_knots = max(int(n_knots), 0)
        if n_knots == 0:
            return np.array([], float)
        return np.linspace(float(w0), float(w1), n_knots + 2)[1:-1]

    def _moon_interline_weights(self) -> np.ndarray:
        line_signal = self._line_density_vector()
        finite = np.isfinite(line_signal)
        positive = finite & (line_signal > 0)
        line_signal_norm = np.zeros_like(line_signal, dtype=float)
        if np.any(positive):
            sig_max = float(np.nanmax(line_signal[positive]))
            if np.isfinite(sig_max) and sig_max > 0:
                line_signal_norm[positive] = line_signal[positive] / sig_max
        # Use a strict binary mask around detected lines; masked pixels are
        # down-weighted while unmasked pixels keep nominal weight.
        raw_line_mask = finite & (line_signal_norm > self.moon_interline_line_flux_threshold)
        line_mask = raw_line_mask.copy()
        if self.moon_interline_exclusion_a > 0 and self.wave.size > 2:
            dw = np.gradient(self.wave)
            dw_med = max(float(np.nanmedian(np.abs(dw))), 1e-6)
            n_pix = int(np.ceil(self.moon_interline_exclusion_a / dw_med))
            if n_pix > 0:
                ker = np.ones(2 * n_pix + 1, dtype=int)
                line_mask = np.convolve(raw_line_mask.astype(int), ker, mode="same") > 0

        red_mask = self.wave >= self.moon_interline_red_min
        masked = line_mask & red_mask

        # Reuse the existing knob: larger values mean stronger suppression of
        # masked pixels, with floor weight 1 / (1 + boost).
        w = np.ones_like(self.wave, dtype=float)
        if self.moon_interline_boost > 0.0:
            w[masked] = 1.0 / (1.0 + self.moon_interline_boost)
        return w

    def _line_density_vector(self) -> np.ndarray:
        out = np.zeros_like(self.wave, dtype=float)
        for name in ("matrix_oh_stick", "matrix_atom_stick", "matrix_orc_stick"):
            mat = getattr(self, name, None)
            if mat is not None and np.size(mat) > 0:
                out += np.nansum(np.abs(np.asarray(mat, float)), axis=0)

        o2_stick = getattr(self, "vector_o2_stick", None)
        if o2_stick is not None and np.size(o2_stick) > 0:
            out += np.abs(np.asarray(o2_stick, float))
        return out

    def _build_diffuse(self) -> tuple[np.ndarray, list[str]]:
        ref = Table.read(self._pmd_path("pmd_refcont.dat"), format="ascii")
        lam_ref = np.asarray(ref["lam"], float) * 1e4
        exact_native_grid = (
            lam_ref.shape == self.wave.shape
            and np.max(np.abs(lam_ref - self.wave)) <= 1.0e-10
        )
        if exact_native_grid:
            vectors = [np.asarray(ref[name], float) for name in ("fcHO2", "fcFeO", "fcO2Ac")]
        else:
            vectors = [
                np.interp(self.wave, lam_ref, np.asarray(ref[name], float))
                for name in ("fcHO2", "fcFeO", "fcO2Ac")
            ]
        self.vector_ho2, self.vector_feo, self.vector_o2ac = [
            np.nan_to_num(vector, nan=0.0, posinf=0.0, neginf=0.0)
            for vector in vectors
        ]

        matrix = np.vstack([self.vector_ho2, self.vector_feo, self.vector_o2ac])
        return matrix, ["HO2", "FeO", "O2Ac"]

    def _build_atom(self) -> tuple[np.ndarray, list[str]]:
        atom = Table.read(self._pmd_path("pmd_intdata_atom.dat"), format="ascii.basic", guess=False, comment="#", fast_reader=False)
        atom["wave"] = vac_to_air(np.asarray(atom["lam"], float) * 1e4)
        mask = (atom["wave"] >= self.wave.min() - CAP_WAVE) & (atom["wave"] <= self.wave.max() + CAP_WAVE)
        atom = atom[mask]
        atom = atom[~np.isin(np.asarray(atom["class"], str), ["H", "Orc"])]
        groups = atom.group_by("class").groups

        matrix = np.zeros((len(groups), self.wave.size))
        self.matrix_atom_stick = np.zeros_like(matrix)
        names = []
        for idx, grp in enumerate(groups):
            amp = np.asarray(grp["I"], float)
            amp /= np.nansum(amp)
            matrix[idx] = grp2vector(grp["wave"], amp, self.wave, self.lsf_sigma)
            self.matrix_atom_stick[idx] = sticks2vector(grp["wave"], amp, self.wave)
            names.append(f"ATOM_{grp['class'][0]}")
        return matrix, names

    def _build_orc(self) -> tuple[np.ndarray, list[str]]:
        orc = Table.read(self._pmd_path("pmd_intmodel_Orc.dat"), format="ascii.basic", guess=False, comment="#", fast_reader=False)
        orc["wave"] = vac_to_air(np.asarray(orc["lam"], float) * 1e4)
        mask = (orc["wave"] >= self.wave.min() - CAP_WAVE) & (orc["wave"] <= self.wave.max() + CAP_WAVE)
        groups = orc[mask].group_by("reffeat").groups

        matrix = np.zeros((len(groups), self.wave.size))
        self.matrix_orc_stick = np.zeros_like(matrix)
        names = []
        for idx, grp in enumerate(groups):
            amp = np.asarray(grp["I"], float)
            matrix[idx] = grp2vector(grp["wave"], amp, self.wave, self.lsf_sigma)
            self.matrix_orc_stick[idx] = sticks2vector(grp["wave"], amp, self.wave)
            names.append(f"ATOM_Orc_{grp['reffeat'][0]}")
        return matrix, names

    def _load_o2_model(self) -> None:
        pop_o2 = Table.read(self._pmd_path("pmd_popmodel_O2.dat"), format="ascii.basic", guess=False, comment="#", fast_reader=False)
        pop_o2["wave"] = vac_to_air(np.asarray(pop_o2["lam"], float) * 1e4)
        o2 = pop_o2[
            (pop_o2["wave"] >= O2_MIN)
            & (pop_o2["wave"] <= O2_MAX)
            & (np.asarray(pop_o2["vi"], int) == 0)
        ]
        self.lam_o2 = np.asarray(o2["wave"], float)
        self.ei_o2 = np.asarray(o2["Ei"], float)
        self.aij_o2 = np.asarray(o2["Aij"], float)
        self.gi_o2 = np.asarray(o2["gi"], float)
        self.e0_o2 = float(np.nanmin(self.ei_o2)) if self.ei_o2.size else 0.0
        self.o2_band = (self.wave >= O2_MIN) & (self.wave <= O2_MAX)

    def _prefit_o2(self, flux: np.ndarray, ivar: np.ndarray) -> None:
        if not np.any(self.o2_band) or self.lam_o2.size == 0:
            self.t_o2 = T_O2_REF
            self.t_o2_err = np.nan
            self.o2_prefit_amp = np.nan
            self.vector_o2 = np.zeros_like(self.wave)
            self.matrix_o2 = self.vector_o2[None, :]
            self.vector_o2_stick = np.zeros_like(self.wave)
            self.matrix_o2_stick = self.vector_o2_stick[None, :]
            self.o2_prefit_bestfit = self.vector_o2.copy()
            self.o2_fit_status = "fallback"
            self.o2_fit_summary = "O2 prefit | status=fallback | reason=no_o2_band"
            self.o2_fit_elapsed_sec = 0.0
            self.o2_valid_frac = 0.0
            return

        x = self.wave[self.o2_band]
        y = flux[self.o2_band]
        iv = ivar[self.o2_band]
        moon = self.vector_moon[self.o2_band]
        base_valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(iv) & (iv > 0) & np.isfinite(moon)
        min_valid = int(np.ceil(self.o2_min_valid_frac * x.size))
        self.o2_valid_frac = float(np.count_nonzero(base_valid) / x.size) if x.size else 0.0
        xv = x[base_valid]
        yv = y[base_valid]
        ivv = iv[base_valid]
        mv = moon[base_valid]
        x0 = float(np.nanmean(xv)) if xv.size else 0.0
        xs = max(float(np.nanmax(xv) - np.nanmin(xv)), 1.0) if xv.size else 1.0
        sqrt_iv = np.sqrt(ivv) if ivv.size else np.array([], float)
        xlin = (xv - x0) / xs if xv.size else np.array([], float)

        def o2_shape(temp: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            rel = self.aij_o2 * self.gi_o2 * np.exp(-HC_OVER_KB_CMK * (self.ei_o2 - self.e0_o2) / temp)
            rel_sum = np.nansum(rel)
            if not np.isfinite(rel_sum) or rel_sum <= 0:
                return np.zeros_like(self.wave), np.zeros_like(self.wave), np.zeros_like(xv)
            rel /= rel_sum
            stick_shape = sticks2vector(self.lam_o2, rel, self.wave)
            full_shape = grp2vector(self.lam_o2, rel, self.wave, self.lsf_sigma)
            return full_shape, stick_shape, full_shape[self.o2_band][base_valid]

        def linear_init(temp: float) -> np.ndarray:
            _, _, sv = o2_shape(temp)
            if sv.size < 3:
                return np.array([0.0, 0.0, 0.0], float)
            design = np.c_[sv, mv, mv * xlin]
            aw = design * sqrt_iv[:, None]
            coef, _, _, _ = np.linalg.lstsq(aw, yv * sqrt_iv, rcond=None)
            return coef

        def set_o2_summary(status: str, dt: float, chi2: float, par: np.ndarray, err: np.ndarray, nvalid: int) -> None:
            dof = max(nvalid - 4, 1)
            chi2_red = chi2 / dof if np.isfinite(chi2) else np.nan
            y_mean = float(np.average(yv, weights=ivv)) if nvalid else np.nan
            sst = float(np.sum((yv - y_mean) ** 2 * ivv)) if nvalid else np.nan
            r2 = 1.0 - chi2 / sst if np.isfinite(sst) and sst > 0 else np.nan
            self.o2_fit_summary = (
                "O2 prefit | "
                f"status={status} | nvalid={nvalid}/{x.size} ({self.o2_valid_frac:.1%}) | "
                f"T={self.t_o2:.2f}+/-{self.t_o2_err:.2f} K | "
                f"chi2_red={chi2_red:.4g} | R2={r2:.5f} | dt={dt:.3f}s"
            )

        init_coef = linear_init(T_O2_REF)

        if np.sum(base_valid) < max(min_valid, 3):
            full_shape, stick_shape, _ = o2_shape(T_O2_REF)
            amp0 = max(float(init_coef[0]), 0.0)
            self.t_o2 = T_O2_REF
            self.t_o2_err = np.nan
            self.o2_prefit_amp = amp0
            self.vector_o2 = np.zeros_like(self.wave)
            self.vector_o2[self.o2_band] = full_shape[self.o2_band]
            self.matrix_o2 = self.vector_o2[None, :]
            self.vector_o2_stick = stick_shape
            self.matrix_o2_stick = self.vector_o2_stick[None, :]
            self.o2_prefit_bestfit = amp0 * self.vector_o2
            self.o2_fit_status = "fallback"
            self.o2_fit_elapsed_sec = 0.0
            set_o2_summary("fallback", 0.0, np.nan, np.r_[self.t_o2, init_coef], np.full(4, np.nan), int(np.sum(base_valid)))
            return

        bounds = ([T_O2_REF - T_O2_HALF_RANGE], [T_O2_REF + T_O2_HALF_RANGE])
        p0 = np.array([T_O2_REF, max(init_coef[0], 0.0), max(init_coef[1], 0.0), init_coef[2]], float)

        def residual(par: np.ndarray) -> np.ndarray:
            _, _, sv = o2_shape(float(par[0]))
            model = par[1] * sv + mv * (par[2] + par[3] * xlin)
            return (yv - model) * sqrt_iv

        t0 = time.perf_counter()
        res = least_squares(
            residual,
            x0=p0,
            bounds=([bounds[0][0], 0.0, 0.0, -np.inf], [bounds[1][0], np.inf, np.inf, np.inf]),
            method="trf",
        )
        dt = time.perf_counter() - t0
        if not res.success or not np.isfinite(res.x[0]):
            full_shape, stick_shape, _ = o2_shape(T_O2_REF)
            amp0 = max(float(init_coef[0]), 0.0)
            self.t_o2 = T_O2_REF
            self.t_o2_err = np.nan
            self.o2_prefit_amp = amp0
            self.vector_o2 = np.zeros_like(self.wave)
            self.vector_o2[self.o2_band] = full_shape[self.o2_band]
            self.matrix_o2 = self.vector_o2[None, :]
            self.vector_o2_stick = stick_shape
            self.matrix_o2_stick = self.vector_o2_stick[None, :]
            self.o2_prefit_bestfit = amp0 * self.vector_o2
            self.o2_fit_status = "fallback"
            self.o2_fit_elapsed_sec = dt
            set_o2_summary("fallback", dt, np.nan, np.r_[self.t_o2, init_coef], np.full(4, np.nan), int(np.sum(base_valid)))
            return

        full_shape, stick_shape, _ = o2_shape(float(res.x[0]))
        self.t_o2 = float(res.x[0])
        jtj = res.jac.T @ res.jac
        dof = max(res.fun.size - res.x.size, 1)
        cov = (2.0 * res.cost / dof) * np.linalg.pinv(jtj)
        err = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))
        self.t_o2_err = float(err[0]) if err.size else np.nan
        self.o2_prefit_amp = float(res.x[1])
        self.vector_o2 = np.zeros_like(self.wave)
        self.vector_o2[self.o2_band] = full_shape[self.o2_band]
        self.matrix_o2 = self.vector_o2[None, :]
        self.vector_o2_stick = stick_shape
        self.matrix_o2_stick = self.vector_o2_stick[None, :]
        self.o2_prefit_bestfit = float(res.x[1]) * self.vector_o2
        self.o2_fit_status = str(res.status)
        self.o2_fit_elapsed_sec = dt
        chi2 = float(np.sum(res.fun**2))
        set_o2_summary(self.o2_fit_status, dt, chi2, res.x, err, int(np.sum(base_valid)))

    def _fit_lsf_channels(self, flux: np.ndarray, ivar: np.ndarray, source: np.ndarray, fixed_background: np.ndarray) -> None:
        offsets = np.arange(-(LSF_KERNEL_SIZE // 2), LSF_KERNEL_SIZE // 2 + 1, dtype=int)

        self.lsf_kernels = {}
        self.lsf_metrics = {}

        for channel, lo, hi in LSF_CHANNELS:
            mask = np.ones_like(self.wave, dtype=bool)
            if lo is not None:
                mask &= self.wave >= lo
            if hi is not None:
                mask &= self.wave < hi

            n_pix = int(np.sum(mask))
            if n_pix == 0:
                continue

            source_ch = np.asarray(source[mask], float)
            obs = np.asarray(flux[mask], float)
            obs_ivar = np.asarray(ivar[mask], float)
            bg = np.asarray(fixed_background[mask], float)
            target = obs - bg
            valid = np.isfinite(source_ch) & np.isfinite(target) & np.isfinite(obs_ivar) & (obs_ivar > 0)
            n_valid = int(np.sum(valid))
            valid_frac = n_valid / n_pix
            min_valid = int(np.ceil(LSF_MIN_VALID_FRAC * n_pix))
            fallback_kernel = self._default_channel_kernel(mask)

            metric: dict[str, object] = {
                "status": "ok",
                "reason": "",
                "n_pixels": n_pix,
                "n_valid_pixels": n_valid,
                "valid_frac": valid_frac,
                "chi2": np.nan,
                "chi2_red": np.nan,
                "r2": np.nan,
                "rms_resid": np.nan,
                "runtime_sec": 0.0,
                "sum_kernel": 1.0,
                "center_pix": 0.0,
                "sigma_pix": 0.0,
            }

            def finalize(kernel: np.ndarray, reason: str = "") -> None:
                kernel = np.asarray(kernel, float)
                model_lsf = self._convolve_with_kernel(source_ch, kernel, offsets)
                resid = target[valid] - model_lsf[valid]
                chi2 = float(np.sum(resid**2 * obs_ivar[valid])) if n_valid else np.nan
                dof = max(n_valid - LSF_KERNEL_SIZE, 1)
                y_mean = float(np.average(target[valid], weights=obs_ivar[valid])) if n_valid else np.nan
                sst = float(np.sum((target[valid] - y_mean) ** 2 * obs_ivar[valid])) if n_valid else np.nan
                metric["reason"] = reason
                metric["chi2"] = chi2
                metric["chi2_red"] = chi2 / dof if np.isfinite(chi2) else np.nan
                metric["r2"] = 1.0 - chi2 / sst if np.isfinite(sst) and sst > 0 else np.nan
                metric["rms_resid"] = float(np.sqrt(np.nanmean(resid**2))) if n_valid else np.nan
                metric["sum_kernel"] = float(np.sum(kernel))
                metric["center_pix"] = float(np.sum(offsets * kernel))
                metric["sigma_pix"] = float(np.sqrt(np.sum((offsets - metric["center_pix"]) ** 2 * kernel)))
                self.lsf_kernels[channel] = kernel
                self.lsf_metrics[channel] = metric

            if n_valid < max(min_valid, 1):
                metric["status"] = "fallback"
                finalize(fallback_kernel, "not_enough_valid_pixels")
                continue

            if not np.isfinite(np.sum(source_ch[valid] ** 2)) or np.sum(source_ch[valid] ** 2) <= 0:
                metric["status"] = "fallback"
                finalize(fallback_kernel, "degenerate_source_model")
                continue

            x = np.column_stack([self._shift_with_zeros(source_ch, int(off)) for off in offsets])
            x_valid = x[valid]
            w = np.sqrt(obs_ivar[valid])
            xw = x_valid * w[:, None]
            yw = target[valid] * w
            data_scale = max(float(np.sqrt(np.nanmean(yw**2))) if yw.size else 1.0, 1.0)
            xw = xw / data_scale
            yw = yw / data_scale

            p_dense = xw.T @ xw
            q = -(xw.T @ yw)
            p = sp.csc_matrix((p_dense + p_dense.T) / 2.0)
            p = sp.triu(p).tocsc()
            aeq = sp.csc_matrix(np.ones((1, LSF_KERNEL_SIZE), dtype=float))
            anonneg = -sp.eye(LSF_KERNEL_SIZE, format="csc")
            acon = sp.vstack([aeq, anonneg], format="csc")
            bcon = np.r_[1.0, np.zeros(LSF_KERNEL_SIZE, dtype=float)]
            cones = [clarabel.ZeroConeT(1), clarabel.NonnegativeConeT(LSF_KERNEL_SIZE)]

            settings = clarabel.DefaultSettings()
            settings.verbose = False

            t0 = time.perf_counter()
            solver = clarabel.DefaultSolver(p, np.asarray(q, dtype=np.float64), acon, bcon, cones, settings)
            qp_result = solver.solve()
            metric["runtime_sec"] = time.perf_counter() - t0

            kernel = np.asarray(qp_result.x, float)
            if str(qp_result.status) != "Solved" or (not np.isfinite(kernel).all()) or float(np.sum(kernel)) <= 0:
                metric["status"] = "fallback"
                finalize(fallback_kernel, "solver_failed")
                continue

            kernel = np.clip(kernel, 0.0, np.inf)
            kernel /= np.sum(kernel)
            metric["status"] = str(qp_result.status)
            finalize(kernel)

    @staticmethod
    def _shift_with_zeros(vec: np.ndarray, offset: int) -> np.ndarray:
        out = np.zeros_like(vec, dtype=float)
        if offset == 0:
            out[:] = vec
        elif offset > 0:
            out[offset:] = vec[:-offset]
        else:
            out[:offset] = vec[-offset:]
        return out

    def _convolve_with_kernel(self, vec: np.ndarray, kernel: np.ndarray, offsets: np.ndarray) -> np.ndarray:
        return self._convolve_rows(np.asarray(vec, float)[None, :], kernel, offsets)[0]

    def _convolve_matrix_channelwise(self, matrix: np.ndarray) -> np.ndarray:
        if matrix.size == 0:
            return matrix.copy()
        out = np.zeros_like(matrix, dtype=float)
        for channel, lo, hi in LSF_CHANNELS:
            mask = np.ones_like(self.wave, dtype=bool)
            if lo is not None:
                mask &= self.wave >= lo
            if hi is not None:
                mask &= self.wave < hi
            if not np.any(mask):
                continue
            kernel = self.lsf_kernels.get(channel, self._default_channel_kernel(mask))
            offsets = np.arange(-(kernel.size // 2), kernel.size // 2 + 1, dtype=int)
            out[:, mask] = self._convolve_rows(np.asarray(matrix[:, mask], float), kernel, offsets)
        return out

    def _default_channel_kernel(self, mask: np.ndarray) -> np.ndarray:
        offsets = np.arange(-(LSF_KERNEL_SIZE // 2), LSF_KERNEL_SIZE // 2 + 1, dtype=float)
        dw = np.gradient(self.wave[mask]) if np.count_nonzero(mask) > 1 else np.array([1.0], float)
        if np.ndim(self.lsf_sigma) > 0:
            sigma_wave = np.nanmedian(np.asarray(self.lsf_sigma[mask], float))
        else:
            sigma_wave = float(self.lsf_sigma)
        sigma_pix = sigma_wave / max(float(np.nanmedian(dw)), 1e-6)
        sigma_pix = max(float(sigma_pix), 0.5)
        kernel = np.exp(-0.5 * (offsets / sigma_pix) ** 2)
        kernel /= np.sum(kernel)
        return kernel

    def _extract_o2_chi2_red(self) -> float:
        token = "chi2_red="
        if token not in self.o2_fit_summary:
            return np.nan
        tail = self.o2_fit_summary.split(token, 1)[1]
        value = tail.split("|", 1)[0].strip()
        try:
            return float(value)
        except ValueError:
            return np.nan

    @staticmethod
    def _fmt_num(value: float) -> str:
        return f"{value:.4g}" if np.isfinite(value) else "nan"

    @staticmethod
    def _convolve_rows(matrix: np.ndarray, kernel: np.ndarray, offsets: np.ndarray) -> np.ndarray:
        matrix = np.asarray(matrix, dtype=float)
        if matrix.size == 0:
            return matrix.copy()
        kernel = np.asarray(kernel, dtype=float)
        offsets = np.asarray(offsets)
        k_size = kernel.size
        half = k_size // 2
        # Fast path when offsets are the canonical symmetric arange used everywhere internally.
        if (
            k_size > 0
            and offsets.size == k_size
            and np.array_equal(offsets.astype(int), np.arange(-half, half + 1))
        ):
            padded = np.pad(matrix, ((0, 0), (half, half)))
            windows = np.lib.stride_tricks.sliding_window_view(padded, k_size, axis=1)
            return np.einsum("riw,w->ri", windows, kernel[::-1])
        out = np.zeros_like(matrix, dtype=float)
        for weight, offset in zip(kernel, offsets):
            o = int(offset)
            if o == 0:
                out += weight * matrix
            elif o > 0:
                out[:, o:] += weight * matrix[:, :-o]
            else:
                out[:, :o] += weight * matrix[:, -o:]
        return out

    def _pmd_path(self, name: str) -> Path:
        suffix = self.palace_suffix
        if name == "pmd_popmodel_OH.dat" and self.palace_oh_suffix is not None:
            suffix = self.palace_oh_suffix
        elif name == "pmd_refcont.dat" and self.palace_diffuse_suffix is not None:
            suffix = self.palace_diffuse_suffix
        if suffix is not None and name in SUFFIXED_PMD_TABLES:
            path = Path(name)
            name = f"{path.stem}{suffix}{path.suffix}"
        return self._require_path(self.pmd_dir / name)

    @staticmethod
    def _validate_palace_suffix(suffix: str | None) -> str | None:
        if suffix is None:
            return None
        if not isinstance(suffix, str):
            raise TypeError("palace_suffix must be a string or None")

        if not suffix:
            raise ValueError("palace_suffix must not be empty")
        if suffix.endswith(".dat"):
            raise ValueError("palace_suffix must not include the .dat extension")
        if (
            re.fullmatch(r"[A-Za-z0-9._-]+", suffix) is None
            or re.search(r"[A-Za-z0-9]", suffix) is None
        ):
            raise ValueError(
                "palace_suffix must contain only letters, digits, '.', '_', or '-'"
            )
        return suffix

    @staticmethod
    def _require_path(path: Path) -> Path:
        if not path.exists():
            raise FileNotFoundError(f"Required reference file not found: {path}")
        return path


def reconstruct_component_spectra(
    wave: np.ndarray,
    coef: np.ndarray,
    lsf_sigma: np.ndarray | float,
    *,
    n_spline_knots: int = 25,
    base_dir: str | Path | None = None,
    palace_suffix: str | None = None,
    palace_oh_suffix: str | None = None,
    palace_diffuse_suffix: str | None = None,
    o2_vector: np.ndarray | None = None,
    coef_err: np.ndarray | None = None,
    moon_interline_boost: float = 0.0,
    moon_interline_red_min: float = 7454.0,
    moon_interline_exclusion_a: float = 3.0,
    moon_interline_line_flux_threshold: float = 0.0,
    split_zodi: bool = False,
    n_zodi_spline_knots: int = 3,
    moon_albedo_asset_path: str | Path | None = None,
    moon_albedo_fiducial_phase_deg: float = 30.0,
    zodi_color_exponent: float = 0.26,
) -> dict[str, np.ndarray]:
    """Reconstruct component spectra from decomposition coefficients.

    Parameters
    ----------
    wave
        Wavelength grid in Angstrom.
    coef
        Coefficient vector matching the internal design-matrix order.
    lsf_sigma
        Gaussian LSF sigma (scalar or per-pixel vector) used to build line bases.
    n_spline_knots
        Number of interior moon B-spline knots.
    base_dir
        Root path containing PALACE/PMD and solar reference files.
    palace_suffix
        Optional suffix selecting versioned ``pmd_popmodel_OH`` and
        ``pmd_refcont`` tables. For example, ``"_adam_v1"`` selects
        ``pmd_popmodel_OH_adam_v1.dat`` and ``pmd_refcont_adam_v1.dat``.
    palace_oh_suffix
        Optional exact suffix overriding ``palace_suffix`` for
        ``pmd_popmodel_OH`` only.
    palace_diffuse_suffix
        Optional exact suffix overriding ``palace_suffix`` for
        ``pmd_refcont`` only.
    o2_vector
        Optional precomputed O2 template on `wave`. If omitted, O2 is set to zero.
    coef_err
        Optional per-coefficient 1σ uncertainty (same shape as ``coef``). When
        provided, the returned dict additionally contains ``sigma`` mapping
        each component name to its per-pixel flux 1σ (first-order Jacobian
        propagation through the LSF-convolved design matrix; see
        ``SkyDecompBase._components_sigma_from_coef_err``), plus
        ``sigma_total`` (quadrature sum of the independent-component
        variances).
    moon_interline_boost
        Strength of down-weighting for masked red pixels in the moon-spline fit.
        Masked-pixel weight is `1 / (1 + moon_interline_boost)`.
    moon_interline_red_min
        Wavelength threshold (Angstrom) above which masking/down-weighting is applied.
    moon_interline_exclusion_a
        Characteristic half-width (Angstrom) used to broaden line influence when
        defining masked line neighborhoods.
    moon_interline_line_flux_threshold
        Normalized threshold in [0, 1] applied to line-signal strength before
        building exclusion windows. Only pixels above this threshold are treated
        as line centers for exclusion.
    split_zodi, n_zodi_spline_knots, moon_albedo_asset_path,
    moon_albedo_fiducial_phase_deg, zodi_color_exponent
        When ``split_zodi=True``, add a separate Zodi_bs block to the design
        matrix (moon carries ROLO albedo × solar, zodi carries a Leinert
        reddening × solar). Must match the settings used when the decomposition
        was fit; the returned dict then also contains a ``"zodi"`` key and
        ``sigma["zodi"]`` when ``coef_err`` is provided.

    Returns
    -------
    dict[str, np.ndarray]
        Component spectra with keys: `oh`, `moon`, `ho2`, `feo`, `o2ac`,
        `diffuse`, `atom`, `orc`, `o2`, and `total`. When ``split_zodi=True``,
        `zodi` is also present.
    """
    model = SkyDecomp(
        wave,
        lsf_sigma=lsf_sigma,
        n_spline_knots=n_spline_knots,
        base_dir=base_dir,
        palace_suffix=palace_suffix,
        palace_oh_suffix=palace_oh_suffix,
        palace_diffuse_suffix=palace_diffuse_suffix,
        moon_interline_boost=moon_interline_boost,
        moon_interline_red_min=moon_interline_red_min,
        moon_interline_exclusion_a=moon_interline_exclusion_a,
        moon_interline_line_flux_threshold=moon_interline_line_flux_threshold,
        split_zodi=split_zodi,
        n_zodi_spline_knots=n_zodi_spline_knots,
        moon_albedo_asset_path=moon_albedo_asset_path,
        moon_albedo_fiducial_phase_deg=moon_albedo_fiducial_phase_deg,
        zodi_color_exponent=zodi_color_exponent,
    )

    coef_arr = np.asarray(coef, float).ravel()

    if o2_vector is not None:
        o2_vec = np.asarray(o2_vector, float).ravel()
        if o2_vec.shape != model.wave.shape:
            raise ValueError(
                f"o2_vector shape mismatch: expected {model.wave.shape}, got {o2_vec.shape}"
            )
        model.vector_o2 = o2_vec
        model.matrix_o2 = o2_vec[None, :]
        model.vector_o2_stick = o2_vec.copy()
        model.matrix_o2_stick = o2_vec[None, :]

    mats = model._matrix_bundle(
        model.matrix_oh,
        model.matrix_moon,
        model.matrix_diffuse,
        model.matrix_atom,
        model.matrix_orc,
        model.matrix_o2,
        matrix_zodi=model.matrix_zodi if split_zodi else None,
    )
    n_expected = sum(m.shape[0] for m in mats.values())
    if coef_arr.size != n_expected:
        raise ValueError(
            f"Coefficient length mismatch: expected {n_expected}, got {coef_arr.size}"
        )

    comps = model._components_from_coef(coef_arr, mats)
    comps["total"] = (
        comps["oh"]
        + comps["moon"]
        + comps["zodi"]
        + comps["diffuse"]
        + comps["atom"]
        + comps["orc"]
        + comps["o2"]
    )
    if coef_err is not None:
        err_arr = np.asarray(coef_err, dtype=np.float64).ravel()
        if err_arr.size != coef_arr.size:
            raise ValueError(
                f"coef_err length mismatch: expected {coef_arr.size}, "
                f"got {err_arr.size}"
            )
        sigma_comps = model._components_sigma_from_coef_err(err_arr, mats)
        comps["sigma"] = sigma_comps
        comps["sigma_total"] = np.sqrt(
            sigma_comps["oh"] ** 2
            + sigma_comps["moon"] ** 2
            + sigma_comps["zodi"] ** 2
            + sigma_comps["diffuse"] ** 2
            + sigma_comps["atom"] ** 2
            + sigma_comps["orc"] ** 2
            + sigma_comps["o2"] ** 2
        )
    return comps


__all__ = [
    "SkyDecomp",
    "SkyDecompResult",
    "vac_to_air",
    "decode_hitran_id",
    "grp2vector",
    "reconstruct_component_spectra",
]
