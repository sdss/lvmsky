"""User-facing pipeline configuration.

The dataclasses here collect every knob the split-zodi notebook exposed as a
module-level constant.  Callers override only the fields they want to change;
domain constants that never varied (physical tables, coefficient schema)
remain hard-coded inside the modules that use them.

Every dataclass is immutable-by-convention: mutate through
``dataclasses.replace(cfg, key=value)`` rather than in-place assignment so
downstream caches (basis wavelength npz, filter mask, split indices) stay
in sync with what produced them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# Default context columns pulled from the science + sky rows.  Mirrors the
# split-zodi notebook's ``context_cols`` list (cell "starter data load").
DEFAULT_CONTEXT_COLUMNS: tuple[str, ...] = (
    "alt",
    "az_sin",
    "az_cos",
    "airmass",
    "moon_sep",
    "moon_alt",
    "moon_az_sin",
    "moon_az_cos",
    "moon_phase_sin",
    "moon_phase_cos",
    "sun_sep",
    "sun_alt",
    "sun_az_sin",
    "sun_az_cos",
    "sci_sep",
    "vanrhijn_87km",
    "vanrhijn_95km",
    "vanrhijn_285km",
    "obstime_day_sin",
    "obstime_day_cos",
    "obstime_lunation_sin",
    "obstime_lunation_cos",
    "obstime_year_sin",
    "obstime_year_cos",
    "f107",
    "f107_81d",
    "kp",
    "ew",
)


@dataclass
class DataConfig:
    """Where the decomposition corpus lives and which flavour to read.

    Attributes
    ----------
    decomp_data_root
        Directory holding the decomposition FITS files (relative to the
        notebook cwd, e.g. ``moon_zodi_spline4``).
    decomp_stem
        Stem shared by the input / sci / sky1 / sky2 FITS files
        (e.g. ``lvmsframe_median_stack_1.2.1_p40_p70``).
    use_lsf_surface_iterative
        Read the LSF-surface-iterative + split-zodi decomposition when
        ``True`` (adds suffix ``_lsf_surface_iterative_split_zodi``);
        the baseline Gaussian-LSF products when ``False``.
    palace_dir
        Base directory used by the reconstruction helpers to locate
        palace assets.
    factor
        Flux rescaling factor applied by the decomposition (1e14).
        Kept here for downstream diagnostics that need the physical
        magnitude back.
    context_columns
        Column names extracted from the META HDU into the per-row context
        vectors ``ctx_near`` / ``ctx_far`` / ``ctx_sci``.
    """

    decomp_data_root: str = "moon_zodi_spline4"
    decomp_stem: str = "lvmsframe_median_stack_1.2.1_p40_p70"
    use_lsf_surface_iterative: bool = True
    palace_dir: str = "../"
    factor: float = 1.0e14
    context_columns: tuple[str, ...] = DEFAULT_CONTEXT_COLUMNS

    @property
    def decomp_suffix(self) -> str:
        return "_lsf_surface_iterative_split_zodi" if self.use_lsf_surface_iterative else ""

    @property
    def decomp_prefix(self) -> str:
        """Full stem including the data-root prefix, e.g.
        ``moon_zodi_spline4/lvmsframe_median_stack_1.2.1_p40_p70``."""
        return f"{self.decomp_data_root}/{self.decomp_stem}"

    @property
    def input_fits_meta(self) -> str:
        return f"{self.decomp_prefix}_meta_only.fits"

    @property
    def input_fits_for_basis(self) -> str:
        """Wavelength / LSF reference used by the wavelength cache builder."""
        return f"{self.decomp_prefix}_every10.fits"

    def coef_fits(self, arm: str) -> str:
        """Path to a per-arm coef-only FITS.  ``arm`` is 'sky1', 'sky2', 'sci'."""
        return f"{self.decomp_prefix}_{arm}_meta_coef{self.decomp_suffix}.fits"

    def decomp_fits(self, arm: str) -> str:
        """Path to a per-arm full-decomp FITS (with COEF_COV_MOON / _ZODI)."""
        return f"{self.decomp_prefix}_decomp_{arm}{self.decomp_suffix}.fits"

    @property
    def wavelength_cache_path(self) -> Path:
        return Path(f"{self.decomp_data_root}/coef_wavelengths_basis_v4.npz")


@dataclass
class FilterConfig:
    """Row-level filters applied after triplet construction.

    Kept coarse for now — deep-dive filter tuning happens in the diagnostic
    cells.  Numeric thresholds match the split-zodi notebook defaults.
    """

    apply_lmc_smc_exclusion: bool = True
    kappa_sigma: float = 5.0
    kappa_mad: float = 5.0


@dataclass
class SplitConfig:
    """Train / val / test split policy."""

    strategy: str = "moon_phase"                 # 'moon_phase' | 'night' | 'random'
    train_frac: float = 0.8
    val_frac: float = 0.1
    seed: int = 42
    n_moon_phase_bins: int = 10


@dataclass
class ModelConfig:
    """Architecture for ``DualEncoderGroupHeadMLP``."""

    encoder_dims: tuple[int, ...] = (384, 192)
    ctx_dims: tuple[int, ...] = (64,)
    trunk_dims: tuple[int, ...] = (256, 128)
    head_dim: int = 128
    drop_vanrhijn_from_context: bool = False


@dataclass
class TrainConfig:
    """Trainer knobs (per-seed).  Ensemble aggregation handled separately."""

    batch_size: int = 256
    epochs: int = 400
    lr: float = 2.5e-4
    weight_decay: float = 1.0e-4
    seed: int = 42
    early_stop_patience: int = 40
    device: str = "auto"                          # 'cpu' | 'cuda' | 'auto'
    n_ensemble_seeds: int = 10
    ensemble_seed_base: int = 1000

    # Loss shaping — split-zodi 2026-08-26 default.
    flux_mse_groups: tuple[str, ...] = ("moon", "zodi")
    flux_mse_eps_frac: float = 0.0
    block_cov_loss_groups: tuple[str, ...] = ()
    relative_mse_groups: tuple[str, ...] = ()

    # Per-group weights that fold into the composite loss.
    group_weights: dict[str, float] = field(
        default_factory=lambda: {
            "moon": 1.0,
            "zodi": 1.0,
            "continuum": 1.0,
            "mesospheric": 1.0,
            "ionospheric": 1.0,
            "atomic": 1.0,
        }
    )

    # Per-group relative sigma floor (WRMSE denominator).
    sigma_floor_by_group: dict[str, float] = field(
        default_factory=lambda: {
            "moon": 0.05,
            "zodi": 0.05,
            "continuum": 0.05,
            "mesospheric": 0.05,
            "ionospheric": 0.05,
            "atomic": 0.05,
        }
    )


@dataclass
class DiagnosticsConfig:
    """Toggles + parameters shared by the diagnostic figure functions."""

    n_worst_rows: int = 15
    n_full_spectrum_sample: int = 100
    plotly_template: str = "plotly_white"


@dataclass
class PipelineConfig:
    """Top-level bundle passed to every stage.

    ``DataPipeline`` / ``Trainer`` / ``Diagnostics`` each receive the whole
    ``PipelineConfig`` and read the sub-config they need, so notebook cells
    only touch one object.
    """

    data: DataConfig = field(default_factory=DataConfig)
    filters: FilterConfig = field(default_factory=FilterConfig)
    split: SplitConfig = field(default_factory=SplitConfig)
    model: ModelConfig = field(default_factory=ModelConfig)
    train: TrainConfig = field(default_factory=TrainConfig)
    diagnostics: DiagnosticsConfig = field(default_factory=DiagnosticsConfig)


__all__ = [
    "DEFAULT_CONTEXT_COLUMNS",
    "DataConfig",
    "FilterConfig",
    "SplitConfig",
    "ModelConfig",
    "TrainConfig",
    "DiagnosticsConfig",
    "PipelineConfig",
]
