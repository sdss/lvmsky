"""Per-coefficient effective wavelength + extinction pipeline.

Wraps the low-level ``coef_wavelengths_from_basis`` / ``resolve_coef_wavelengths_a``
/ ``fit_effective_extinction`` / ``resolve_coef_extinction_k`` helpers from
:mod:`data` into a small orchestrator so notebook cells only call a single
entry point.

Cells consolidated:
- ``infer-spline-knots``   : ``infer_spline_knots(coef_names)``.
- ``wavelength-cache-gate``: cache gate around ``coef_wavelengths_from_basis``.
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

from .data import (
    _build_group_indices,
    assert_context_is_physical,
    coef_wavelengths_from_basis,
    fit_effective_extinction,
    resolve_coef_extinction_k,
    resolve_coef_wavelengths_a,
)


def infer_spline_knots(coef_names: list[str]) -> tuple[int, bool, int]:
    """Auto-infer B-spline knot counts from the actual coef_names.

    Returns ``(n_moon_knots, split_zodi, n_zodi_knots)``.  The cubic B-spline
    convention is ``n_basis = n_interior_knots + 4``.
    """

    moon_bs_count = sum(1 for n in coef_names if str(n).startswith("Moon_bs"))
    zodi_bs_count = sum(1 for n in coef_names if str(n).startswith("Zodi_bs"))
    n_moon_knots = moon_bs_count - 4
    split_zodi = zodi_bs_count > 0
    n_zodi_knots = (zodi_bs_count - 4) if split_zodi else 0
    return n_moon_knots, split_zodi, n_zodi_knots


def wavelength_cache_matches_corpus(
    cache_path: Path,
    coef_names: list[str],
) -> bool:
    """Return True iff the cache at ``cache_path`` has the expected schema and coef names."""

    cache_path = Path(cache_path)
    if not cache_path.exists():
        return False
    with np.load(cache_path, allow_pickle=True) as c:
        if "coef_names" not in c.files or "k_eff_a" not in c.files:
            return False
        cached = [str(x) for x in c["coef_names"]]
    corpus = [str(x) for x in coef_names]
    return cached == corpus


def populate_wavelength_cache(
    *,
    cache_path: Path,
    input_fits_for_basis: str | Path,
    coef_names: list[str],
    n_moon_knots: int,
    split_zodi: bool,
    n_zodi_knots: int,
    palace_oh_suffix: str | None = None,
    palace_diffuse_suffix: str | None = None,
    verbose: bool = True,
) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    """Idempotently populate the wavelength+k_eff cache; return the two arrays.

    If the cache is stale or missing, computes it (one reconstruction call per
    coefficient); otherwise loads the existing arrays.  Returns ``(None, None)``
    only if the FITS reference is missing the required extensions and
    computation fails.
    """

    cache_path = Path(cache_path)
    if wavelength_cache_matches_corpus(cache_path, coef_names):
        with np.load(cache_path, allow_pickle=True) as c:
            n = len(c["coef_names"])
            lam_ok = int(np.isfinite(c["wavelengths_a"]).sum())
            k_ok = int(np.isfinite(c["k_eff_a"]).sum())
            wavelengths_a = np.asarray(c["wavelengths_a"], dtype=np.float64)
            k_eff_a = np.asarray(c["k_eff_a"], dtype=np.float64)
        if verbose:
            print(f"[wavelength-cache] {cache_path}:")
            print(
                f"                   {lam_ok}/{n} finite centroids, "
                f"{k_ok}/{n} finite k_eff (loaded, no compute)"
            )
        return wavelengths_a, k_eff_a

    if cache_path.exists() and verbose:
        print(
            f"[wavelength-cache] {cache_path.name} exists but does not "
            f"match current corpus; will overwrite."
        )
    if verbose:
        print(
            f"[wavelength-cache] populating for {len(coef_names)} coefs "
            f"(one reconstruction call each, ~13 min total)..."
        )
    with fits.open(str(input_fits_for_basis)) as hdul:
        ext_names = [h.name for h in hdul]
        if "WAVE" not in ext_names:
            raise KeyError(f"{input_fits_for_basis} has no WAVE extension")
        wave_ref = np.asarray(hdul["WAVE"].data, dtype=np.float64)
        lsf_name = "LSF_SCI" if "LSF_SCI" in ext_names else ("LSF" if "LSF" in ext_names else None)
        if lsf_name is None:
            raise KeyError(f"{input_fits_for_basis} has no LSF_SCI or LSF extension")
        lsf_ref = np.asarray(hdul[lsf_name].data, dtype=np.float64)
    if wave_ref.ndim > 1:
        wave_ref = wave_ref[0]
    if lsf_ref.ndim > 1:
        lsf_ref = lsf_ref[0]

    t0 = time.perf_counter()
    result = coef_wavelengths_from_basis(
        coef_names=coef_names,
        wave=wave_ref,
        lsf_sigma=lsf_ref / 2.35,
        n_spline_knots=n_moon_knots,
        split_zodi=split_zodi,
        n_zodi_spline_knots=n_zodi_knots,
        palace_oh_suffix=palace_oh_suffix,
        palace_diffuse_suffix=palace_diffuse_suffix,
        cache_path=cache_path,
        return_k_eff=True,
        verbose=verbose,
    )
    if verbose:
        print(f"[wavelength-cache] populated in {time.perf_counter() - t0:.0f} s")
    if isinstance(result, tuple):
        return result
    return result, None


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
    cache_path: Path,
    input_fits_for_basis: str | Path,
    use_fitted_extinction: bool = True,
    palace_oh_suffix: str | None = None,
    palace_diffuse_suffix: str | None = None,
    n_wavelength_bins: int = 8,
    verbose: bool = True,
) -> ExtinctionResolution:
    """Populate ``filtered_triplet`` with ``coef_wavelengths_a`` + ``coef_extinction_k``.

    This is the notebook-facing entry point that mirrors cell ``98a37092``'s
    end-to-end sequence: infer knots, populate the wavelength cache, resolve
    per-coef wavelengths, verify context is physical, fit LCO-effective
    extinction, resolve per-coef extinction.

    Mutates ``filtered_triplet`` in place (adds the two arrays), and returns
    an :class:`ExtinctionResolution` bundle with all intermediates.
    """

    coef_names = list(filtered_triplet["coef_names"])
    n_moon_knots, split_zodi, n_zodi_knots = infer_spline_knots(coef_names)
    if verbose:
        print("Auto-inferred spline settings from filtered_triplet['coef_names']:")
        print(f"  n_spline_knots = {n_moon_knots}, split_zodi = {split_zodi}, "
              f"n_zodi_spline_knots = {n_zodi_knots}")

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
        basis_lam, basis_k = populate_wavelength_cache(
            cache_path=cache_path,
            input_fits_for_basis=input_fits_for_basis,
            coef_names=coef_names,
            n_moon_knots=n_moon_knots,
            split_zodi=split_zodi,
            n_zodi_knots=n_zodi_knots,
            palace_oh_suffix=palace_oh_suffix,
            palace_diffuse_suffix=palace_diffuse_suffix,
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
    "infer_spline_knots",
    "populate_wavelength_cache",
    "resolve_wavelengths_and_extinction",
    "wavelength_cache_matches_corpus",
]
