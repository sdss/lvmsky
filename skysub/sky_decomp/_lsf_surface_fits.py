"""Private FITS schema helpers for compact wavelength-dependent LSF states."""

from __future__ import annotations

import numpy as np

from .fit import LSF_CHANNELS


CHANNEL_NAMES = tuple(channel for channel, _, _ in LSF_CHANNELS)

STATE_CONFIG_COLUMNS = (
    ("line_weight", "line_weight", float),
    ("skyline_cumulative_fraction", "skyline_cumulative_fraction", float),
    ("skyline_half_width_angstrom", "skyline_half_width_angstrom", float),
    ("huber_transition_sigma", "huber_transition_sigma", float),
    ("n_refinement_cycles", "requested_cycles", int),
    ("n_basis", "config_n_basis", int),
    ("degree", "config_degree", int),
    ("roughness_fraction", "roughness_fraction", float),
    ("fallback_prior_fraction", "fallback_prior_fraction", float),
    ("information_prior_max_boost", "information_prior_max_boost", float),
    ("background_degree", "background_degree", int),
    ("blue_fit_lower", "blue_fit_lower", float),
)

METRIC_COLUMNS = (
    ("status", "", str),
    ("reason", "", str),
    ("n_pixels", 0, int),
    ("n_valid_pixels", 0, int),
    ("chi2_red", np.nan, float),
    ("rms_resid", np.nan, float),
    ("center_pix_min", np.nan, float),
    ("center_pix_median", np.nan, float),
    ("center_pix_max", np.nan, float),
    ("sigma_pix_min", np.nan, float),
    ("sigma_pix_median", np.nan, float),
    ("sigma_pix_max", np.nan, float),
    ("runtime_sec", 0.0, float),
    ("fit_lower", np.nan, float),
    ("fit_upper", np.nan, float),
    ("knots_fixed", False, bool),
    ("model", "", str),
)


def _meta_row(spectrum_index, state, channel):
    coefficient = state.coefficients[channel]
    knot_vector = state.knot_vectors[channel]
    metric = state.metrics[channel]
    lower, upper = state.channel_bounds[channel]
    row = {
        "spectrum_index": spectrum_index,
        "channel": channel,
        "available": True,
        "lower": np.nan if lower is None else lower,
        "upper": np.nan if upper is None else upper,
        "degree": state.degrees[channel],
        "n_basis": coefficient.shape[1],
        "n_knots": knot_vector.size,
        "tap_offsets": state.tap_offsets.copy(),
        "status": str(metric.get("status", "")),
        "reason": str(metric.get("reason", "")),
        "schema_version": state.schema_version,
        "requested_cycles": state.requested_cycles,
        "completed_cycles": state.completed_cycles,
        "fit_status": state.fit_status,
        "failure_reason": state.failure_reason,
        "final_continuum_status": state.final_continuum_status,
        "final_line_status": state.final_line_status,
        "knot_strategy": state.knot_strategy,
        "legacy_kernel_representation": state.legacy_kernel_representation,
        "wave_n": state.wave_n,
        "wave_min": state.wave_min,
        "wave_max": state.wave_max,
        "wave_sha256": state.wave_sha256,
    }
    for name, column, converter in STATE_CONFIG_COLUMNS:
        if name != "n_refinement_cycles":  # requested_cycles is already stored above.
            row[column] = converter(state.config[name])
    for name, default, converter in METRIC_COLUMNS[2:]:
        row[name] = converter(metric.get(name, default))
    return row


def build_lsf_hdus(states):
    """Validate compact states and return their three ordered FITS HDUs."""
    from astropy.io import fits
    from astropy.table import Table

    reference = states[0]
    for state in states:
        if state.schema_version != reference.schema_version:
            raise ValueError("All LSF states must use the same schema version")
        if not np.array_equal(state.tap_offsets, reference.tap_offsets):
            raise ValueError("All LSF states must use the same tap offsets")
        if (
            state.wave_n != reference.wave_n
            or state.wave_min != reference.wave_min
            or state.wave_max != reference.wave_max
            or state.wave_sha256 != reference.wave_sha256
        ):
            raise ValueError("All LSF states must describe the same wavelength grid")
        if set(state.coefficients) != set(CHANNEL_NAMES):
            raise ValueError("Every LSF state must contain B, R, and Z channels")
    max_basis = max(
        state.coefficients[channel].shape[1] for state in states for channel in CHANNEL_NAMES
    )
    max_knots = max(
        state.knot_vectors[channel].size for state in states for channel in CHANNEL_NAMES
    )
    coefficient_cube = np.full(
        (len(states), len(CHANNEL_NAMES), reference.tap_offsets.size, max_basis),
        np.nan,
        dtype=float,
    )
    knot_cube = np.full(
        (len(states), len(CHANNEL_NAMES), max_knots),
        np.nan,
        dtype=float,
    )
    meta_rows = []
    for spectrum_index, state in enumerate(states):
        for channel_index, channel in enumerate(CHANNEL_NAMES):
            coefficient = state.coefficients[channel]
            knot_vector = state.knot_vectors[channel]
            coefficient_cube[spectrum_index, channel_index, :, : coefficient.shape[1]] = (
                coefficient
            )
            knot_cube[spectrum_index, channel_index, : knot_vector.size] = knot_vector
            meta_rows.append(_meta_row(spectrum_index, state, channel))

    coefficient_hdu = fits.ImageHDU(coefficient_cube, name="LSF_COEF")
    coefficient_hdu.header["TAPMIN"] = int(reference.tap_offsets[0])
    coefficient_hdu.header["TAPMAX"] = int(reference.tap_offsets[-1])
    return [
        coefficient_hdu,
        fits.ImageHDU(knot_cube, name="LSF_KNOTS"),
        fits.BinTableHDU(Table(rows=meta_rows), name="LSF_META"),
    ]
