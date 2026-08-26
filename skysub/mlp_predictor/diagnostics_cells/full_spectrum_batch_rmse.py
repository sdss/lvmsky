# Full spectrum 100 rows, RMSE stats
import numpy as np
import pandas as pd
import plotly.graph_objects as go

RUN_RMSE_SUBSET_EVAL = True  # Set True to execute this slower evaluation cell.

required = [
    "mlp_artifacts",
    "predict_sci_coefficients_default",
    "context_cols",
    "build_triplet_coef_dataset",
    "reconstruct_component_spectra",
    "reconstruct_with_lsf",
    "load_lsf_state_if_available",
    "load_o2_vector_if_available",
    "_infer_base_dir_for_reconstruction",
    "FACTOR",
]
missing = [k for k in required if k not in globals()]
if missing:
    raise RuntimeError("Run the training + residual-correction cells first. Missing: " + ", ".join(missing))

if not RUN_RMSE_SUBSET_EVAL:
    print("Cell 19 skipped. Set RUN_RMSE_SUBSET_EVAL = True to run the random-subset RMSE evaluation.")
else:
    # Use the same file inputs as the reconstruction diagnostic cell.
    _e10_stem = f"{DECOMP_DATA_ROOT}/{DECOMP_STEM}_every10"
    _e10_suffix = _DECOMP_SUFFIX  # inherited from cell 6
    EVERY10_INPUT = f"{_e10_stem}.fits"
    EVERY10_NEAR = f"{_e10_stem}_decomp_sky1{_e10_suffix}.fits"
    EVERY10_FAR = f"{_e10_stem}_decomp_sky2{_e10_suffix}.fits"
    EVERY10_SCI = f"{_e10_stem}_decomp_sci{_e10_suffix}.fits"

    n_sample = 100
    rng_seed = 42

    # 1) Build aligned triplet rows, keeping chi2 so we can apply the same
    #    quality gates the training set went through.
    e10_triplet = build_triplet_coef_dataset(
        input_fits_path=EVERY10_INPUT,
        sky_near_decomp_fits_path=EVERY10_NEAR,
        sky_far_decomp_fits_path=EVERY10_FAR,
        sci_decomp_fits_path=EVERY10_SCI,
        context_columns=context_cols,
        return_chi2=True,
    )
    # ECLIPTIC-CTX-V1: match training-time ctx layout on the e10 triplet.
    if '_augment_triplet_with_ecliptic' in globals():
        _augment_triplet_with_ecliptic(e10_triplet, meta_fits_path=EVERY10_INPUT)
    # PHYSICS-PRIORS-CTX-V1: same augment on the e10 triplet.
    if '_augment_triplet_with_physics_priors' in globals():
        _augment_triplet_with_physics_priors(e10_triplet)

    row_index_e10 = np.asarray(e10_triplet["row_index"], dtype=np.int64)
    _e10_n0 = int(e10_triplet["n_rows"])
    if _e10_n0 == 0:
        raise RuntimeError("No aligned rows available in e10_triplet")

    # 1a) Apply chi2 gating and LMC/SMC field exclusion BEFORE the random
    #     subsample, so the reconstructions we score are drawn from tiles
    #     that pass the same data-quality gates the training set went
    #     through.  Hard coefficient bounds and kappa-sigma clipping are
    #     intentionally NOT applied here -- those are training-time filters
    #     on the target that would remove the model's hardest true cases
    #     from the evaluation set.
    _e10_keep = np.ones(_e10_n0, dtype=bool)

    if "sci_ra" in e10_triplet and "sci_dec" in e10_triplet:
        _sci_ra = np.asarray(e10_triplet["sci_ra"], dtype=np.float64)
        _sci_dec = np.asarray(e10_triplet["sci_dec"], dtype=np.float64)
        _field_mask = np.ones(_e10_n0, dtype=bool)
        for _region in (LMC_EXCLUSION, SMC_EXCLUSION):
            _sep = _angular_separation_deg_vec(
                _sci_ra, _sci_dec,
                _region["ra_deg"], _region["dec_deg"])
            _inside = np.isfinite(_sep) & (_sep <= float(_region["radius_deg"]))
            print(
                f"  every10 field exclusion around {_region['name']}: "
                f"excluded {int(_inside.sum())}/{_e10_n0}"
            )
            _field_mask &= ~_inside
        _e10_keep &= _field_mask
    else:
        print("  every10 field exclusion skipped: no sci_ra/sci_dec in triplet.")

    if all(_k in e10_triplet for _k in ("chi2_near", "chi2_far", "chi2_sci")):
        # 2026-08-19: use max (not nanmax) so any-arm-NaN chi2 propagates
        # NaN and disqualifies the whole observation via the isfinite gate.
        _chi2_combined = np.max(
            np.column_stack([
                np.asarray(e10_triplet["chi2_near"], dtype=np.float64),
                np.asarray(e10_triplet["chi2_far"], dtype=np.float64),
                np.asarray(e10_triplet["chi2_sci"], dtype=np.float64),
            ]),
            axis=1,
        )
        _chi2_finite = _chi2_combined[np.isfinite(_chi2_combined)]
        _chi2_hi = (float(np.nanpercentile(_chi2_finite, 90.0))
                    if _chi2_finite.size else np.inf)
        _chi2_upper = min(10.0, _chi2_hi)
        _chi2_mask = (np.isfinite(_chi2_combined)
                      & (_chi2_combined >= 0.0)
                      & (_chi2_combined <= _chi2_upper))
        print(
            f"  every10 chi2 filter (qmax=90% -> {_chi2_hi:.3g}, "
            f"upper={_chi2_upper:.3g}): "
            f"kept {int(_chi2_mask.sum())}/{_e10_n0}"
        )
        _e10_keep &= _chi2_mask
    else:
        print("  every10 chi2 filter skipped: chi2 columns not in triplet.")

    _e10_valid_pos = np.flatnonzero(_e10_keep)
    n_rows = int(_e10_valid_pos.size)
    print(
        f"  every10 rows passing chi2 + field exclusion: "
        f"{n_rows}/{_e10_n0} ({100.0 * n_rows / max(_e10_n0, 1):.1f}%)"
    )
    if n_rows == 0:
        raise RuntimeError(
            "No aligned rows survive chi2/field filtering; relax thresholds.")

    n_use = int(min(n_sample, n_rows))
    rng = np.random.default_rng(rng_seed)

    # Stratify the draw by lunar phase so the diagnostic set spans dark -> bright
    # roughly uniformly (matches split_indices_by_moon_phase in §8.1). A plain
    # rng.choice over _e10_valid_pos would inherit the phase distribution of the
    # dataset -- weighted toward whichever quantiles happen to hold more filtered
    # rows -- and the sci_pred_vs_true RMSE would then be dominated by that
    # region rather than being representative of deployment conditions.
    _e10_moon_phase = _moon_phase_deg_from_ctx(e10_triplet)
    _valid_phase = _e10_moon_phase[_e10_valid_pos]
    if not np.isfinite(_valid_phase).all():
        raise RuntimeError('Non-finite moon_phase in the valid every10 rows; '
                           'cannot stratify by lunar phase.')

    _n_phase_bins = int(min(10, n_use))
    _phase_edges = np.quantile(_valid_phase,
                               np.linspace(0.0, 1.0, _n_phase_bins + 1))
    _phase_edges[0], _phase_edges[-1] = -np.inf, np.inf
    _bin_id = np.digitize(_valid_phase, _phase_edges[1:-1], right=False)

    # Round-robin quota with the +1s scattered randomly so no bin is systematically favored.
    _quota = np.full(_n_phase_bins, n_use // _n_phase_bins, dtype=int)
    _quota[:n_use - int(_quota.sum())] += 1
    rng.shuffle(_quota)

    _picked = []
    for _b in range(_n_phase_bins):
        _in_bin = _e10_valid_pos[_bin_id == _b]
        _take = int(min(_quota[_b], _in_bin.size))
        if _take > 0:
            _picked.append(rng.choice(_in_bin, size=_take, replace=False))
    _selected = (np.concatenate(_picked).astype(int)
                 if _picked else np.array([], dtype=int))

    # If any bin was smaller than its quota, backfill from the remaining pool.
    _shortfall = n_use - _selected.size
    if _shortfall > 0:
        _remaining = np.setdiff1d(_e10_valid_pos, _selected, assume_unique=False)
        if _remaining.size >= _shortfall:
            _selected = np.concatenate(
                [_selected, rng.choice(_remaining, size=_shortfall, replace=False)])

    sel_pos = np.sort(_selected)
    sel_rows = row_index_e10[sel_pos]

    _sel_phases = _e10_moon_phase[sel_pos]
    print(f"  phase-stratified sample: n_use={n_use} across {_n_phase_bins} "
          f"quantile bins; phase deg quartiles (min / 25 / 50 / 75 / max) = "
          f"{float(np.min(_sel_phases)):.1f} / "
          f"{float(np.percentile(_sel_phases, 25)):.1f} / "
          f"{float(np.percentile(_sel_phases, 50)):.1f} / "
          f"{float(np.percentile(_sel_phases, 75)):.1f} / "
          f"{float(np.max(_sel_phases)):.1f}")

    # 2) Load observed spectra, wavelength grid, and LSF from the same input file.
    with fits.open(EVERY10_INPUT) as hdul:
        flux_near_all = np.asarray(hdul["FLUX_SKY_NEAR"].data, dtype=np.float64)
        flux_far_all = np.asarray(hdul["FLUX_SKY_FAR"].data, dtype=np.float64)
        flux_sci_all = np.asarray(hdul["FLUX_SCI"].data, dtype=np.float64)
        wave_arr = np.asarray(hdul["WAVE"].data, dtype=np.float64)
        lsf_sci_arr = np.asarray(hdul["LSF_SCI"].data, dtype=np.float64)
        # expnum per every10 row for the per-line hover tooltip on the residual figure.
        _expnum_all = None
        if "META" in hdul:
            _meta_e10 = Table(hdul["META"].data)
            _meta_up_e10 = {c.upper(): c for c in _meta_e10.colnames}
            _expnum_col = next(
                (_meta_up_e10[k] for k in ("EXPNUM", "EXP_NUM", "EXPOSURE")
                 if k in _meta_up_e10), None)
            if _expnum_col is not None:
                _expnum_all = np.asarray(_meta_e10[_expnum_col])

    n_spec = int(flux_sci_all.shape[0])
    if np.any(sel_rows < 0) or np.any(sel_rows >= n_spec):
        raise ValueError("Selected row index is outside the valid range of EVERY10_INPUT")

    # 3) Predict SCI coefficients for sampled rows.
    coef_sci_pred = predict_sci_coefficients_default(
        mlp_artifacts,
        coef_near_phys=e10_triplet["coef_near"][sel_pos],
        coef_far_phys=e10_triplet["coef_far"][sel_pos],
        ctx_near_phys=e10_triplet["ctx_near"][sel_pos],
        ctx_far_phys=e10_triplet["ctx_far"][sel_pos],
        ctx_sci_phys=e10_triplet["ctx_sci"][sel_pos],
    ).astype(np.float64)

    coef_near_sel = np.asarray(e10_triplet["coef_near"][sel_pos], dtype=np.float64)
    coef_far_sel = np.asarray(e10_triplet["coef_far"][sel_pos], dtype=np.float64)
    coef_sci_sel = np.asarray(e10_triplet["coef_sci"][sel_pos], dtype=np.float64)

    # 4) Reconstruct and compute per-row RMSE + pixel-space WRMSE.
    base_dir_guess = _infer_base_dir_for_reconstruction()
    near_rmse = np.full(n_use, np.nan, dtype=np.float64)
    far_rmse = np.full(n_use, np.nan, dtype=np.float64)
    sci_rmse = np.full(n_use, np.nan, dtype=np.float64)
    near_wrmse = np.full(n_use, np.nan, dtype=np.float64)
    far_wrmse = np.full(n_use, np.nan, dtype=np.float64)
    sci_wrmse = np.full(n_use, np.nan, dtype=np.float64)

    # Try to load per-pixel sigma HDUs from the decomposition FITS files
    # up-front (fast path when new decompositions land).  Falls back to
    # on-the-fly propagation via coef_err inside the loop when absent.
    _pix_sigma_near_all = load_pixel_sigma_if_available(EVERY10_NEAR)
    _pix_sigma_far_all  = load_pixel_sigma_if_available(EVERY10_FAR)
    _pix_sigma_sci_all  = load_pixel_sigma_if_available(EVERY10_SCI)
    _pix_sigma_source = {
        arm: ("FITS HDU" if arr is not None else "coef_err propagation")
        for arm, arr in (("near", _pix_sigma_near_all),
                         ("far",  _pix_sigma_far_all),
                         ("sci",  _pix_sigma_sci_all))
    }
    print(f"  pixel sigma source: near={_pix_sigma_source['near']}, "
          f"far={_pix_sigma_source['far']}, sci={_pix_sigma_source['sci']}")

    # Grab coef_err arrays for the fallback path (may be all-NaN when the
    # decomposition FITS lacks a COEF_ERR HDU; the WRMSE helper falls back
    # to floor-only weighting so nothing breaks).
    _e10_coef_err_near = (np.asarray(e10_triplet.get("coef_err_near",
                                     np.full_like(coef_near_sel, np.nan)),
                                     dtype=np.float64)[sel_pos]
                          if _pix_sigma_near_all is None else None)
    _e10_coef_err_far  = (np.asarray(e10_triplet.get("coef_err_far",
                                     np.full_like(coef_far_sel, np.nan)),
                                     dtype=np.float64)[sel_pos]
                          if _pix_sigma_far_all is None else None)
    _e10_coef_err_sci  = (np.asarray(e10_triplet.get("coef_err_sci",
                                     np.full_like(coef_sci_pred, np.nan)),
                                     dtype=np.float64)[sel_pos]
                          if _pix_sigma_sci_all is None else None)
    sci_resid_rows = []
    sci_wave_rows = []
    # Per-component residuals: comps_sci[<comp>] - comps_sci_true[<comp>] per row,
    # stored in native units (like sci_resid_rows) and multiplied by FACTOR at plot time.
    sci_moon_resid_rows = []
    sci_zodi_resid_rows = []
    sci_diffuse_resid_rows = []
    sci_lines_resid_rows = []

    def _lines_sum(_c):
        return (np.asarray(_c["oh"], dtype=np.float64)
                + np.asarray(_c["atom"], dtype=np.float64)
                + np.asarray(_c["orc"], dtype=np.float64)
                + np.asarray(_c["o2"], dtype=np.float64))

    # 4a) Optimization (2026-08-11): build ONE reconstruction model outside the loop and
    #     precache per-file LSF-availability + full VECTOR_O2 cubes. Previously the loop
    #     rebuilt SkyDecompLSFSurfaceIterative (basis + solar-reference + moon spline)
    #     3 x n_use times and re-read VECTOR_O2 on every call, both of which dominated
    #     the runtime. See §12 (2026-08-11 batch-RMSE cell reconstruction hoisted).
    import time as _time
    _t_recon0 = _time.perf_counter()
    _wave_ref_recon = (wave_arr if wave_arr.ndim == 1
                       else np.asarray(wave_arr[int(sel_rows[0])], dtype=np.float64))
    if wave_arr.ndim > 1:
        _wave_probe = np.asarray(wave_arr[int(sel_rows[-1])], dtype=np.float64)
        if _wave_probe.shape != _wave_ref_recon.shape or not np.allclose(
                _wave_probe, _wave_ref_recon, rtol=0.0, atol=1e-8):
            raise RuntimeError(
                "wave_arr rows differ between sampled rows; model hoisting assumes a shared grid. "
                "Fall back to per-row reconstruct_with_lsf if this ever triggers on your dataset.")
    _lsf_model = SkyDecompLSFSurfaceIterative(
        _wave_ref_recon, lsf_sigma=1.0, n_spline_knots=N_MOON_KNOTS, base_dir=base_dir_guess,
        palace_oh_suffix='_joint_v2_updated',
        palace_diffuse_suffix='_joint_native_adam_invsky_p2_10000iter',
        split_zodi=SPLIT_ZODI, n_zodi_spline_knots=N_ZODI_KNOTS,
    )

    def _precache_decomp_state(decomp_path):
        state = {"path": Path(decomp_path), "has_lsf": False, "o2_cube": None}
        if not state["path"].exists():
            return state
        try:
            with fits.open(str(state["path"]), memmap=False) as _hdul_dec:
                _ext_names = {h.name for h in _hdul_dec}
                state["has_lsf"] = all(_e in _ext_names
                                       for _e in ("LSF_COEF", "LSF_KNOTS", "LSF_META"))
                if state["has_lsf"]:
                    # Refuse to open an inconsistent LSF cube (observed on the
                    # _p25_every10 decomp files, where LSF_COEF was rewritten
                    # to every10 size but LSF_META was inherited from every1).
                    # Loading anyway would slice cube row k with an n_basis
                    # taken from META row k -- a different fiber -- silently
                    # applying the wrong LSF to ~99% of rows and hitting NaN
                    # padding on the rest.
                    _coef_rows = int(_hdul_dec["LSF_COEF"].data.shape[0])
                    _meta_rows = int(len(_hdul_dec["LSF_META"].data))
                    _expected_meta = 3 * _coef_rows
                    if _meta_rows != _expected_meta:
                        raise RuntimeError(
                            f"Inconsistent LSF HDUs in {state['path'].name}: "
                            f"LSF_COEF has {_coef_rows} rows but LSF_META has "
                            f"{_meta_rows} rows (expected 3 x {_coef_rows} = "
                            f"{_expected_meta}). Regenerate this decomposition "
                            f"file with matching LSF_META; do NOT fall back to "
                            f"LSF_SCI sigma because the per-row LSF would be "
                            f"silently wrong on all other rows."
                        )
                if "VECTOR_O2" in _ext_names:
                    _data = np.asarray(_hdul_dec["VECTOR_O2"].data, dtype=np.float64)
                    if _data.ndim == 2:
                        state["o2_cube"] = _data
        except (KeyError, IndexError, ValueError) as _exc:
            print(f"  precache failed for {state['path'].name}: "
                  f"{type(_exc).__name__}: {_exc}")
        return state

    _state_near = _precache_decomp_state(EVERY10_NEAR)
    _state_far  = _precache_decomp_state(EVERY10_FAR)
    _state_sci  = _precache_decomp_state(EVERY10_SCI)
    print(f"  precache: has_lsf near/far/sci = "
          f"{_state_near['has_lsf']}/{_state_far['has_lsf']}/{_state_sci['has_lsf']}, "
          f"VECTOR_O2 near/far/sci = "
          f"{_state_near['o2_cube'] is not None}/"
          f"{_state_far['o2_cube'] is not None}/"
          f"{_state_sci['o2_cube'] is not None}")

    def _lsf_state_from_cache(state_dict, row_idx):
        if not state_dict["has_lsf"]:
            return None
        try:
            return load_lsf_surface_state(str(state_dict["path"]), int(row_idx))
        except (KeyError, IndexError, ValueError) as _exc:
            print(f"  LSF surface state unavailable in {state_dict['path'].name} "
                  f"row {int(row_idx)}: {type(_exc).__name__}: {_exc}")
            return None

    def _o2_vec_from_cache(state_dict, row_idx):
        cube = state_dict["o2_cube"]
        if cube is None or int(row_idx) >= cube.shape[0]:
            return None
        row = cube[int(row_idx)]
        if not np.isfinite(row).any() or float(np.nansum(np.abs(row))) == 0.0:
            return None
        return row

    def _fast_reconstruct(coef, lsf_state, o2_vec, lsf_sigma_fallback, coef_err=None):
        if isinstance(lsf_state, LSFSurfaceState):
            _lsf_model._set_lsf_state(lsf_state)
            _mats = _lsf_model._assemble_refined_matrices()
            if o2_vec is not None:
                _o2_arr = np.asarray(o2_vec, float).ravel()
                if _o2_arr.shape != _lsf_model.wave.shape:
                    raise ValueError(
                        f"o2_vector shape mismatch: expected {_lsf_model.wave.shape}, "
                        f"got {_o2_arr.shape}")
                _mats["o2"] = _o2_arr[None, :]
            _coef_arr = np.asarray(coef, float).ravel()
            _comps = _lsf_model._components_from_coef(_coef_arr, _mats)
            _comps["total"] = (_comps["oh"] + _comps["moon"] + _comps.get("zodi", 0) + _comps["diffuse"]
                                + _comps["atom"] + _comps["orc"] + _comps["o2"])
            if coef_err is not None:
                _err_arr = np.asarray(coef_err, float).ravel()
                _sigmas = _lsf_model._components_sigma_from_coef_err(_err_arr, _mats)
                _comps["sigma"] = _sigmas
                _comps["sigma_total"] = np.sqrt(
                    _sigmas["oh"] ** 2 + _sigmas["moon"] ** 2
                    + _sigmas.get("zodi", np.zeros_like(_sigmas["moon"])) ** 2
                    + _sigmas["diffuse"] ** 2 + _sigmas["atom"] ** 2
                    + _sigmas["orc"] ** 2 + _sigmas["o2"] ** 2)
            return _comps
        # Rare fallback path: no LSF surface state for this row -- take the slow route.
        return reconstruct_with_lsf(
            wave=_lsf_model.wave, coef=coef, lsf=lsf_sigma_fallback,
            n_spline_knots=N_MOON_KNOTS, base_dir=base_dir_guess, o2_vector=o2_vec,
            split_zodi=SPLIT_ZODI, n_zodi_spline_knots=N_ZODI_KNOTS,
            coef_err=coef_err)

    print(f"  recon setup: {_time.perf_counter() - _t_recon0:.2f} s "
          f"(one-time basis build + FITS precache)")

    _t_loop0 = _time.perf_counter()
    for i, r in enumerate(sel_rows):
        rr = int(r)
        wave_row = wave_arr if wave_arr.ndim == 1 else np.asarray(wave_arr[rr], dtype=np.float64)
        lsf_row = lsf_sci_arr if lsf_sci_arr.ndim == 1 else np.asarray(lsf_sci_arr[rr], dtype=np.float64)

        flux_near_true = np.asarray(flux_near_all[rr], dtype=np.float64)
        flux_far_true = np.asarray(flux_far_all[rr], dtype=np.float64)
        flux_sci_true = np.asarray(flux_sci_all[rr], dtype=np.float64)

        _lsf_state_near = _lsf_state_from_cache(_state_near, rr)
        _lsf_state_far  = _lsf_state_from_cache(_state_far,  rr)
        _lsf_state_sci  = _lsf_state_from_cache(_state_sci,  rr)
        _lsf_sigma_fallback = lsf_row / 2.35

        _o2_vec_near = _o2_vec_from_cache(_state_near, rr)
        _o2_vec_far  = _o2_vec_from_cache(_state_far,  rr)
        _o2_vec_sci  = _o2_vec_from_cache(_state_sci,  rr)

        # For each arm, only ask the reconstructor for sigma when the
        # loaded FITS sigma is absent; otherwise the propagator call would
        # duplicate work already done by the pipeline.
        _cerr_near_row = (_e10_coef_err_near[i]
                          if _e10_coef_err_near is not None else None)
        _cerr_far_row  = (_e10_coef_err_far[i]
                          if _e10_coef_err_far is not None else None)
        _cerr_sci_row  = (_e10_coef_err_sci[i]
                          if _e10_coef_err_sci is not None else None)

        comps_near = _fast_reconstruct(coef_near_sel[i], _lsf_state_near,
                                        _o2_vec_near, _lsf_sigma_fallback,
                                        coef_err=_cerr_near_row)
        comps_far  = _fast_reconstruct(coef_far_sel[i],  _lsf_state_far,
                                        _o2_vec_far,  _lsf_sigma_fallback,
                                        coef_err=_cerr_far_row)
        comps_sci  = _fast_reconstruct(coef_sci_pred[i], _lsf_state_sci,
                                        _o2_vec_sci,  _lsf_sigma_fallback,
                                        coef_err=_cerr_sci_row)

        flux_near_recon = np.asarray(comps_near["total"], dtype=np.float64) / FACTOR
        flux_far_recon = np.asarray(comps_far["total"], dtype=np.float64) / FACTOR
        flux_sci_pred = np.asarray(comps_sci["total"], dtype=np.float64) / FACTOR

        # nanmean so isolated NaN pixels (~1 pixel/row on ~40% of every10) don't
        # poison the pRMSE; a whole-row veto lives with the chi2/field filter above.
        near_rmse[i] = float(np.sqrt(np.nanmean((flux_near_recon - flux_near_true) ** 2)))
        far_rmse[i] = float(np.sqrt(np.nanmean((flux_far_recon - flux_far_true) ** 2)))

        sci_resid = flux_sci_pred - flux_sci_true
        sci_rmse[i] = float(np.sqrt(np.nanmean(sci_resid ** 2)))
        sci_resid_rows.append(np.asarray(sci_resid, dtype=np.float64))
        sci_wave_rows.append(np.asarray(wave_row, dtype=np.float64))

        # Reconstruct the sci-arm spectrum from the FITTED sci coefficients so
        # per-component residuals (pred - recon(sci_true)) can be separated in
        # the multi-panel residual plot below (matches cell 26 fig_deltas).
        comps_sci_true_batch = _fast_reconstruct(
            coef_sci_sel[i], _lsf_state_sci, _o2_vec_sci, _lsf_sigma_fallback,
            coef_err=None,
        )
        _dmoon    = (np.asarray(comps_sci["moon"], dtype=np.float64)
                     - np.asarray(comps_sci_true_batch["moon"], dtype=np.float64)) / FACTOR
        # comps_*.get('zodi') is missing on pre-split rows, so fall back to zero for backward compatibility.
        _zodi_pred_arr = np.asarray(comps_sci.get("zodi", 0.0), dtype=np.float64)
        _zodi_true_arr = np.asarray(comps_sci_true_batch.get("zodi", 0.0), dtype=np.float64)
        _dzodi    = (_zodi_pred_arr - _zodi_true_arr) / FACTOR
        _ddiffuse = (np.asarray(comps_sci["diffuse"], dtype=np.float64)
                     - np.asarray(comps_sci_true_batch["diffuse"], dtype=np.float64)) / FACTOR
        _dlines   = (_lines_sum(comps_sci) - _lines_sum(comps_sci_true_batch)) / FACTOR
        sci_moon_resid_rows.append(_dmoon)
        sci_zodi_resid_rows.append(_dzodi)
        sci_diffuse_resid_rows.append(_ddiffuse)
        sci_lines_resid_rows.append(_dlines)

        # Pixel-space WRMSE: prefer the FITS-side sigma if present, else use
        # the propagator output from _fast_reconstruct.  All three sources
        # deliver the same LSF-aware sigma; the fallback of last resort is
        # the median-floor path inside pixel_wrmse_per_row.
        _sig_near_row = (_pix_sigma_near_all[int(sel_pos[i])] * FACTOR
                         if _pix_sigma_near_all is not None
                         else comps_near.get("sigma_total"))
        _sig_far_row  = (_pix_sigma_far_all[int(sel_pos[i])] * FACTOR
                         if _pix_sigma_far_all is not None
                         else comps_far.get("sigma_total"))
        _sig_sci_row  = (_pix_sigma_sci_all[int(sel_pos[i])] * FACTOR
                         if _pix_sigma_sci_all is not None
                         else comps_sci.get("sigma_total"))
        # comps_*['sigma_total'] comes out in native units; the flux_* arrays
        # here are already divided by FACTOR, so the sigma from propagation
        # must be divided by FACTOR too for a scale match.
        if _pix_sigma_near_all is None and _sig_near_row is not None:
            _sig_near_row = np.asarray(_sig_near_row) / FACTOR
        if _pix_sigma_far_all is None and _sig_far_row is not None:
            _sig_far_row = np.asarray(_sig_far_row) / FACTOR
        if _pix_sigma_sci_all is None and _sig_sci_row is not None:
            _sig_sci_row = np.asarray(_sig_sci_row) / FACTOR

        near_wrmse[i] = float(pixel_wrmse_per_row(
            flux_near_recon, flux_near_true, _sig_near_row)[0])
        far_wrmse[i]  = float(pixel_wrmse_per_row(
            flux_far_recon,  flux_far_true,  _sig_far_row)[0])
        sci_wrmse[i]  = float(pixel_wrmse_per_row(
            flux_sci_pred,   flux_sci_true,  _sig_sci_row)[0])

    print(f"  recon loop:  {_time.perf_counter() - _t_loop0:.2f} s "
          f"({n_use} rows x 3 arms = {3 * n_use} reconstructions with hoisted basis)")

    def _rmse_stats(arr):
        x = np.asarray(arr, dtype=np.float64)
        x = x[np.isfinite(x)]
        if x.size == 0:
            return {
                "count": 0,
                "mean": np.nan,
                "median": np.nan,
                "std": np.nan,
                "min": np.nan,
                "p05": np.nan,
                "p95": np.nan,
                "max": np.nan,
            }
        return {
            "count": int(x.size),
            "mean": float(np.mean(x)),
            "median": float(np.median(x)),
            "std": float(np.std(x)),
            "min": float(np.min(x)),
            "p05": float(np.percentile(x, 5.0)),
            "p95": float(np.percentile(x, 95.0)),
            "max": float(np.max(x)),
        }

    summary_df = pd.DataFrame(
        [
            {"series": "near_self_recon_pRMSE",  **_rmse_stats(near_rmse)},
            {"series": "near_self_recon_pWRMSE", **_rmse_stats(near_wrmse)},
            {"series": "far_self_recon_pRMSE",   **_rmse_stats(far_rmse)},
            {"series": "far_self_recon_pWRMSE",  **_rmse_stats(far_wrmse)},
            {"series": "sci_pred_vs_true_pRMSE", **_rmse_stats(sci_rmse)},
            {"series": "sci_pred_vs_true_pWRMSE",**_rmse_stats(sci_wrmse)},
        ]
    )

    summary_disp_df = summary_df.copy()
    for c in ["mean", "median", "std", "min", "p05", "p95", "max"]:
        summary_disp_df[c] = summary_disp_df[c] * FACTOR

    print(f"Random per-row pRMSE / pWRMSE evaluation on {n_use} spectra from every10 inputs (seed={rng_seed})")
    print("Per-row pixel-space pRMSE / pWRMSE stats in physical flux units:")
    print(summary_df.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("")
    print(f"Per-row pRMSE / pWRMSE stats in display units (x{FACTOR:.3g}):")
    print(summary_disp_df.to_string(index=False, float_format=lambda v: f"{v:.6g}"))

    # Multi-panel residual figure (2026-08-19): row 1 keeps the historical
    # sci total residual (pred - observed); rows 2-4 show per-component
    # residuals (pred - recon(sci_true)) so a broadband deficit that lives
    # entirely in one component (e.g. moon spline) shows up separately from
    # a line-emission miss (mesospheric / atomic / ionospheric / O2).  The
    # legend is off; per-line hover shows row_idx + expnum from META.
    from plotly.subplots import make_subplots as _make_subplots_resid

    sci_resid_arr = np.vstack(sci_resid_rows) * FACTOR
    sci_moon_arr = np.vstack(sci_moon_resid_rows) * FACTOR
    sci_zodi_arr = np.vstack(sci_zodi_resid_rows) * FACTOR
    sci_diffuse_arr = np.vstack(sci_diffuse_resid_rows) * FACTOR
    sci_lines_arr = np.vstack(sci_lines_resid_rows) * FACTOR
    wave_ref = sci_wave_rows[0]
    same_grid = all(
        (w.shape == wave_ref.shape) and np.allclose(w, wave_ref, rtol=0.0, atol=1e-8)
        for w in sci_wave_rows[1:]
    )

    fig_resid = _make_subplots_resid(
        rows=5, cols=2,
        shared_xaxes=False,
        column_widths=[0.82, 0.18],
        horizontal_spacing=0.03,
        vertical_spacing=0.04,
        subplot_titles=(
            f"SCI residuals: pred - observed (n={n_use})",
            "median residual / row (±1,2,3σ)",
            "Moon component: pred - recon(sci coef)",
            "",
            "Zodi component: pred - recon(sci coef)",
            "",
            "Diffuse continuum (HO2 + FeO + O2ac): pred - recon(sci coef)",
            "",
            "Lines (OH + atom + ORC + O2): pred - recon(sci coef)",
            "",
        ),
    )
    # Share the wavelength x-axis across the four spectrum panels only; the
    # right-column histograms keep independent x-axes because per-component
    # median residuals live on very different scales.
    for _r in range(2, 6):
        fig_resid.update_xaxes(matches="x", row=_r, col=1)

    def _expnum_str(i):
        if _expnum_all is None:
            return ""
        try:
            _e = int(_expnum_all[int(sel_rows[i])])
        except (IndexError, ValueError, TypeError):
            return ""
        return f" | expnum {_e}"

    _panels = [
        ("total",   sci_resid_arr),
        ("moon",    sci_moon_arr),
        ("zodi",    sci_zodi_arr),
        ("diffuse", sci_diffuse_arr),
        ("lines",   sci_lines_arr),
    ]
    # (percent enclosed, color, dash) for the per-row median-residual histograms.
    _sigma_specs = [
        (68.27, "rgba(200, 60, 60, 0.95)", "solid"),
        (95.45, "rgba(230,140, 40, 0.85)", "dash"),
        (99.73, "rgba( 90, 90, 90, 0.75)", "dot"),
    ]
    for _row_i, (_pname, _arr) in enumerate(_panels, start=1):
        if same_grid:
            for i in range(n_use):
                _rid = int(sel_rows[i])
                _hover = (
                    f"row {_rid}{_expnum_str(i)}<br>"
                    f"λ=%{{x:.1f}} Å<br>"
                    f"Δ_{_pname}=%{{y:.4g}}"
                    "<extra></extra>"
                )
                fig_resid.add_trace(
                    go.Scattergl(
                        x=wave_ref, y=_arr[i],
                        mode="lines",
                        line=dict(width=0.8),
                        opacity=0.7,
                        hovertemplate=_hover,
                        showlegend=False,
                    ),
                    row=_row_i, col=1,
                )
            _rms_band = np.sqrt(np.mean(_arr ** 2, axis=0))
            fig_resid.add_trace(
                go.Scatter(
                    x=wave_ref, y=-_rms_band,
                    mode="lines",
                    line=dict(color="rgba(120,120,120,0.6)", width=1.0),
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=_row_i, col=1,
            )
            fig_resid.add_trace(
                go.Scatter(
                    x=wave_ref, y=_rms_band,
                    mode="lines",
                    line=dict(color="rgba(120,120,120,0.6)", width=1.0),
                    fill="tonexty",
                    fillcolor="rgba(120,120,120,0.15)",
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=_row_i, col=1,
            )
        else:
            for i in range(n_use):
                _rid = int(sel_rows[i])
                _hover = (
                    f"row {_rid}{_expnum_str(i)}<br>"
                    f"λ=%{{x:.1f}} Å<br>"
                    f"Δ_{_pname}=%{{y:.4g}}"
                    "<extra></extra>"
                )
                fig_resid.add_trace(
                    go.Scattergl(
                        x=sci_wave_rows[i], y=_arr[i],
                        mode="lines",
                        line=dict(width=0.8),
                        opacity=0.7,
                        hovertemplate=_hover,
                        showlegend=False,
                    ),
                    row=_row_i, col=1,
                )
            _global_rms = float(np.sqrt(np.mean(_arr ** 2)))
            fig_resid.add_hrect(
                y0=-_global_rms, y1=_global_rms,
                fillcolor="rgba(120,120,120,0.15)",
                line_width=0, layer="above",
                row=_row_i, col=1,
            )
        fig_resid.add_hline(
            y=0.0, line=dict(color="rgba(0,0,0,0.5)", width=0.8, dash="dash"),
            row=_row_i, col=1,
        )

        # Right-column histogram: median residual per row (median over pixels)
        # with empirical two-tailed 1σ / 2σ / 3σ percentile bars.
        _med_per_row = np.nanmedian(_arr, axis=1)
        _med_per_row = _med_per_row[np.isfinite(_med_per_row)]
        if _med_per_row.size > 0:
            _nb_hist = int(max(10, min(40, np.sqrt(_med_per_row.size) * 2.0)))
            fig_resid.add_trace(
                go.Histogram(
                    x=_med_per_row,
                    nbinsx=_nb_hist,
                    marker=dict(
                        color="rgba(70,110,170,0.65)",
                        line=dict(color="rgba(35, 60,120,1.0)", width=0.4),
                    ),
                    hovertemplate=("median Δ_" + _pname
                                   + "=%{x:.4g}<br>count=%{y}<extra></extra>"),
                    showlegend=False,
                ),
                row=_row_i, col=2,
            )
            _row_median = float(np.median(_med_per_row))
            fig_resid.add_vline(
                x=_row_median,
                line=dict(color="rgba(0,0,0,0.7)", width=0.9, dash="dashdot"),
                row=_row_i, col=2,
            )
            for _pct, _color, _dash in _sigma_specs:
                _q_lo = (100.0 - _pct) / 2.0
                _q_hi = 100.0 - _q_lo
                _v_lo = float(np.percentile(_med_per_row, _q_lo))
                _v_hi = float(np.percentile(_med_per_row, _q_hi))
                fig_resid.add_vline(
                    x=_v_lo,
                    line=dict(color=_color, width=1.0, dash=_dash),
                    row=_row_i, col=2,
                )
                fig_resid.add_vline(
                    x=_v_hi,
                    line=dict(color=_color, width=1.0, dash=_dash),
                    row=_row_i, col=2,
                )

    fig_resid.update_xaxes(title_text="Wavelength [Å]", row=4, col=1)
    fig_resid.update_xaxes(title_text="median residual / row", row=4, col=2)
    fig_resid.update_yaxes(title_text="pred - obs", row=1, col=1)
    fig_resid.update_yaxes(title_text="pred - recon(true) [moon]",    row=2, col=1)
    fig_resid.update_yaxes(title_text="pred - recon(true) [zodi]",    row=3, col=1)
    fig_resid.update_yaxes(title_text="pred - recon(true) [diffuse]", row=4, col=1)
    fig_resid.update_yaxes(title_text="pred - recon(true) [lines]",   row=5, col=1)
    fig_resid.update_layout(
        template="plotly_white",
        title=(f"SCI + per-component residuals (n={n_use} spectra, "
               f"display units x{FACTOR:.3g}). Hover any line for row_idx + expnum. "
               f"Gray band = ± RMS(λ) across rows; right-column histograms "
               f"show the per-row median residual with empirical ±1σ/2σ/3σ percentiles."),
        height=1500,
        margin=dict(l=80, r=20, t=110, b=60),
        showlegend=False,
        bargap=0.05,
    )
    fig_resid.add_annotation(
        xref="paper", yref="paper",
        x=1.0, y=1.02,
        xanchor="right", yanchor="bottom",
        showarrow=False,
        font=dict(size=11),
        text=("<span style='color:rgb(200,60,60)'>─ 1σ</span>"
              " &nbsp; <span style='color:rgb(230,140,40)'>-- 2σ</span>"
              " &nbsp; <span style='color:rgb(90,90,90)'>·· 3σ</span>"
              " &nbsp; <span style='color:rgb(50,50,50)'>-·- median</span>"),
    )
    fig_resid.show()

    rmse_subset_results = {
        "row_positions": sel_pos,
        "row_indices": sel_rows,
        "near_rmse": near_rmse,
        "far_rmse": far_rmse,
        "sci_rmse": sci_rmse,
        "near_wrmse": near_wrmse,
        "far_wrmse": far_wrmse,
        "sci_wrmse": sci_wrmse,
        "pix_sigma_source": _pix_sigma_source,
        "summary": summary_df,
        "summary_display": summary_disp_df,
        "sci_residuals": sci_resid_arr,
    }
