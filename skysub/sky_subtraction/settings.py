"""Operating points of the SkyML fit, recorded in provenance and reused by reconstruction."""

from dataclasses import asdict, dataclass, fields


@dataclass(frozen=True)
class Settings:
    # Gaia-clean fiber selection and medians.
    gaia_ratio_threshold: float = 0.01
    science_min_fibers: int = 50
    min_valid_fraction: float = 0.95
    # Per-fiber R/Z LSF fit.
    continuum_spacing_angstrom: float = 20.0
    max_iterations: int = 12
    amplitude_prior_sigma: float = 0.3
    variance_floor: float = 0.03
    huber_threshold: float = 2.5
    offset_roughness_fraction: float = 0.01
    fallback_prior_fraction: float = 0.001
    active_line_fraction: float = 0.01
    min_line_pixels_per_coefficient: int = 5
    min_pixels_per_basis: int = 5
    n_probes: int = 9
    convergence_width: float = 0.01
    convergence_centroid_angstrom: float = 0.015
    # Profile guards against the survey reference and the Science median.
    lsf_max_relative_width_change: float = 0.35
    lsf_max_centroid_change_angstrom: float = 0.25
    lsf_max_median_width_change: float = 0.5
    lsf_fwhm_range_angstrom: tuple = (0.9, 3.0)
    # Observed line-information gate, see lsf_information.py.
    line_bins: int = 4
    line_min_per_bin: int = 2
    line_min_snr: float = 5.0
    line_min_contrast: float = 1.0
    line_min_separation_angstrom: float = 10.0
    line_core_angstrom: float = 3.0
    line_side_angstrom: tuple = (4.0, 10.0)
    line_min_pixels: int = 6
    line_min_core_fraction: float = 0.8
    # Rendering.
    residual_max_angstrom: float = 5000.0
    fwhm_probes: int = 33

    def asdict(self):
        return asdict(self)

    @classmethod
    def from_record(cls, record):
        """Recorded values override defaults; unknown keys (older products) are ignored."""
        names = {field.name for field in fields(cls)}
        values = {
            key: tuple(value) if isinstance(value, list) else value
            for key, value in (record or {}).items()
            if key in names
        }
        return cls(**values)


SETTINGS = Settings()
