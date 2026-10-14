"""Per-coefficient effective wavelength + extinction pipeline.

Wraps the low-level ``coef_wavelengths_from_basis`` / ``resolve_coef_wavelengths_a``
/ ``fit_effective_extinction`` / ``resolve_coef_extinction_k`` helpers from
:mod:`data` into a small orchestrator so notebook cells only call a single
entry point.

Cells consolidated:
- ``infer-spline-knots``   : ``infer_spline_knots(coef_names)``.
- ``wavelength-cache-gate``: the gate around ``coef_wavelengths_from_basis``
                             (the on-disk cache it guarded was removed once the
                             design-matrix path brought the build down to ~1 s).
- ``98a37092``             : the full "populate coef wavelengths + extinction"
                             sequence (basis wavelengths, wavelength resolution,
                             physical-context assertion, fitted extinction,
                             extinction resolution).
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from sky_decomp.moon_zodi_model import LSF_FWHM_TO_SIGMA

from .data import (
    _build_group_indices,
    assert_context_is_physical,
    coef_wavelengths_from_basis,
    fit_effective_extinction,
    resolve_coef_extinction_k,
    resolve_coef_wavelengths_a,
)


def infer_spline_knots(coef_names: list[str]) -> tuple[int, int]:
    """Auto-infer B-spline knot counts from the actual coef_names.

    Returns ``(n_moon_knots, n_zodi_knots)``.  The cubic B-spline convention
    is ``n_basis = n_interior_knots + 4``.
    """

    moon_bs_count = sum(1 for n in coef_names if str(n).startswith("Moon_bs"))
    zodi_bs_count = sum(1 for n in coef_names if str(n).startswith("Zodi_bs"))
    return moon_bs_count - 4, zodi_bs_count - 4


def build_coef_wavelengths(
    *,
    input_fits_for_basis: str | Path,
    coef_names: list[str],
    n_moon_knots: int,
    n_zodi_knots: int,
    verbose: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Build the per-coefficient wavelength centroids and k_eff for the corpus.

    Nothing is cached: the centroid and the B^2-weighted ``k_eff`` are read
    straight off the decomposer's design matrix in ~1 s for the 388
    coefficients.

    The basis is divided by a PER-ROW transmission, but one set of wavelengths
    serves the whole corpus, so the row with the median sci_airmass stands in:
    the transmission enters a centroid only as a smooth multiplicative weight
    (median T = 0.9997), so the choice moves the centroids far less than the
    line grouping does.  The LSF width is the stack's first LSF_SCI row.
    """
    from . import data as _data

    t0 = time.perf_counter()
    with fits.open(str(input_fits_for_basis)) as hdul:
        wave_ref = np.asarray(hdul["WAVE"].data, dtype=np.float64)
        lsf_ref = np.asarray(hdul["LSF_SCI"].data, dtype=np.float64)
        meta = Table(hdul["META"].data)
    wave_ref = wave_ref if wave_ref.ndim == 1 else wave_ref[0]
    lsf_all = lsf_ref
    lsf_ref = lsf_ref if lsf_ref.ndim == 1 else lsf_ref[0]
    rep = _data.telluric_representative_row(input_fits_for_basis)
    tel = _data.telluric_row_kwargs(meta, rep, 'sci', wave_ref, lsf_all[rep])
    decomposer = _data.make_reconstruction_decomposer(
        wave_ref, n_spline_knots=n_moon_knots, n_zodi_spline_knots=n_zodi_knots,
        base_dir=_data._infer_base_dir_for_reconstruction(), telluric=tel,
        # AT CONSTRUCTION. Assigning `.lsf_sigma` afterwards does not rebuild
        # the design matrix, and this basis IS the design matrix: built at the
        # 1.0 A default it put the OH centroids a median 13.3 A out.
        lsf_sigma=lsf_ref / LSF_FWHM_TO_SIGMA)
    if verbose:
        print(f'[wavelength-basis] telluric basis from representative row {rep}, '
              f'sci_airmass {float(meta["sci_airmass"][rep]):.3f}, pwv {tel["pwv_mm"]:g} mm')
    result = coef_wavelengths_from_basis(coef_names, wave_ref, decomposer, verbose=verbose)
    if verbose:
        print(f"[wavelength-basis] built in {time.perf_counter() - t0:.1f} s")
    return result


@dataclass
class ExtinctionResolution:
    """Bundle returned by :func:`resolve_wavelengths_and_extinction`.

    Attributes
    ----------
    coef_wavelengths_a
        Per-coefficient effective wavelength in Angstrom (n_coef,).
    coef_extinction_k
        Per-coefficient effective extinction (mag/airmass), same shape.
    coef_wavelengths_basis
        Raw basis-derived centroids returned by
        :func:`coef_wavelengths_from_basis`, or ``None`` if unavailable.
    coef_k_eff_basis
        Raw basis-derived B²-weighted k_eff, or ``None`` if unavailable.
    coef_wavelength_source
        Per-coefficient string tag describing which fallback rule
        produced the wavelength (``basis``, ``group_default``, ...).
    coef_extinction_source
        Per-coefficient string tag for k origin (``fit_table``, ``lco_generic``,
        ...).
    extinction_fit_table
        DataFrame produced by :func:`fit_effective_extinction`, empty when
        ``use_fitted_extinction=False``.
    group_indices
        ``{group_name: ndarray}`` returned by ``_build_group_indices``.
    """

    coef_wavelengths_a: np.ndarray
    coef_extinction_k: np.ndarray
    coef_wavelengths_basis: np.ndarray | None
    coef_k_eff_basis: np.ndarray | None
    coef_wavelength_source: np.ndarray
    coef_extinction_source: np.ndarray
    extinction_fit_table: pd.DataFrame
    group_indices: dict[str, np.ndarray]


def resolve_wavelengths_and_extinction(
    filtered_triplet: dict[str, Any],
    *,
    input_fits_for_basis: str | Path,
    use_fitted_extinction: bool = True,
    n_wavelength_bins: int = 8,
    verbose: bool = True,
) -> ExtinctionResolution:
    """Populate ``filtered_triplet`` with ``coef_wavelengths_a`` + ``coef_extinction_k``.

    This is the notebook-facing entry point that mirrors cell ``98a37092``'s
    end-to-end sequence: infer knots, build the basis wavelengths, resolve
    per-coef wavelengths, verify context is physical, fit LCO-effective
    extinction, resolve per-coef extinction.

    Mutates ``filtered_triplet`` in place (adds the two arrays), and returns
    an :class:`ExtinctionResolution` bundle with all intermediates.
    """

    coef_names = list(filtered_triplet["coef_names"])
    n_moon_knots, n_zodi_knots = infer_spline_knots(coef_names)
    if verbose:
        print("Auto-inferred spline settings from filtered_triplet['coef_names']:")
        print(f"  n_spline_knots = {n_moon_knots}, n_zodi_spline_knots = {n_zodi_knots}")

    group_indices = _build_group_indices(coef_names)
    grp_sizes = {g: int(idx.size) for g, idx in group_indices.items()}
    if verbose:
        print(
            f"Coefficient group sizes: {grp_sizes} (total {sum(grp_sizes.values())})"
        )

    t_start = time.perf_counter()

    def _lap(label: str) -> None:
        nonlocal t_start
        if verbose:
            print(f"    [{label}: {time.perf_counter() - t_start:.1f} s]")
        t_start = time.perf_counter()

    try:
        basis_lam, basis_k = build_coef_wavelengths(
            input_fits_for_basis=input_fits_for_basis,
            coef_names=coef_names,
            n_moon_knots=n_moon_knots,
            n_zodi_knots=n_zodi_knots,
            verbose=verbose,
        )
    except Exception as exc:
        if verbose:
            print(
                f"Basis wavelengths unavailable ({type(exc).__name__}: {exc}); "
                f"falling back to name-token and group defaults."
            )
        basis_lam, basis_k = None, None
    _lap("basis wavelengths + B^2 k_eff")

    coef_wavelengths_a, coef_wavelength_source = resolve_coef_wavelengths_a(
        coef_names,
        group_indices=group_indices,
        basis_wavelengths_a=basis_lam,
        verbose=verbose,
    )
    filtered_triplet["coef_wavelengths_a"] = coef_wavelengths_a
    _lap("wavelength resolution")

    for label, ctx in (
        ("near", filtered_triplet["ctx_near"]),
        ("far", filtered_triplet["ctx_far"]),
        ("sci", filtered_triplet["ctx_sci"]),
    ):
        assert_context_is_physical(ctx, filtered_triplet["ctx_names"])
        if verbose:
            print(f"context[{label}] verified physical (van Rhijn columns consistent with alt)")

    extinction_fit_table = pd.DataFrame()
    if use_fitted_extinction:
        extinction_fit_table = fit_effective_extinction(
            coef_near=filtered_triplet["coef_near"],
            coef_far=filtered_triplet["coef_far"],
            ctx_near=filtered_triplet["ctx_near"],
            ctx_far=filtered_triplet["ctx_far"],
            ctx_names=filtered_triplet["ctx_names"],
            group_indices=group_indices,
            coef_wavelengths_a=coef_wavelengths_a,
            n_wavelength_bins=n_wavelength_bins,
            verbose=verbose,
        )

    coef_extinction_k, coef_extinction_source = resolve_coef_extinction_k(
        coef_names,
        coef_wavelengths_a,
        group_indices,
        fit_table=extinction_fit_table if use_fitted_extinction else None,
        clip_to_generic=True,
        coef_basis_k_generic=basis_k,
        verbose=verbose,
    )
    filtered_triplet["coef_extinction_k"] = coef_extinction_k
    _lap("effective-extinction fit")

    return ExtinctionResolution(
        coef_wavelengths_a=coef_wavelengths_a,
        coef_extinction_k=coef_extinction_k,
        coef_wavelengths_basis=basis_lam,
        coef_k_eff_basis=basis_k,
        coef_wavelength_source=coef_wavelength_source,
        coef_extinction_source=coef_extinction_source,
        extinction_fit_table=extinction_fit_table,
        group_indices=group_indices,
    )


__all__ = [
    "ExtinctionResolution",
    "build_coef_wavelengths",
    "infer_spline_knots",
    "resolve_wavelengths_and_extinction",
]
