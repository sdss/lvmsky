"""Cached per-arm predictions from the frozen physical Moon/Zodi model.

Why a cache exists at all
-------------------------
``MoonZodiPhysicalModel.predict`` costs ~0.25 s per arm, so evaluating a whole
corpus is ~3 h of wall clock for three arms x ~14 500 rows.  That is far too
slow to sit inside the training loop, but the numbers are a pure function of
(expnum, pointing, LSF) and never change once the decomposition products are
written -- so they are computed once, in parallel, and stored beside the corpus.

What it is for
--------------
The model's moon TRANSFER RATIO, ``moon(sci) / moon(near)``, predicts the
science moon amplitude better than the trained network does.  On the gaia-stars
test split (held out by night, moon-up, n = 583), error in dex:

    flat, sci = near arm       MAD 0.0228   p90 0.0747
    model, near x r_model      MAD 0.0119   p90 0.0466
    the ML network             MAD 0.0139   p90 0.0475

and the model wins in every stratum of |log10 r_model| except one holding 5
rows.  rho(log r_model, log true ratio) = +0.794, against +0.067 for the crude
``moon_signal_proxy`` the network already receives -- so this is new
INFORMATION, not a rearrangement of existing inputs, which is what separates it
from the six refuted attempts recorded in ``ablations.RETIRED``.

NOTE this reverses an earlier refutation of the same ratio.  That test used
new-oh-2 TAIL rows only and reported the opposite ordering; a 53-row version of
the test above also pointed the wrong way.  Do not re-refute it on a small or
selected subsample.

What is cached
--------------
Per row and per arm (``near``, ``far``, ``sci``): the wavelength-integrated
moon and zodi predictions in fit units, plus the geometry the model derived --
``moon_separation_deg``, ``target_airmass``, ``moon_airmass``,
``moon_altitude_deg``, ``sun_altitude_deg``, ``zodi_b500``.  Rows the model
cannot handle store NaN and ``ok = False`` rather than being dropped, so the
cache is always aligned one-to-one with the corpus rows.

The integrals are stored, not the spectra: a per-arm spectrum cache would be
14 500 x 3 x 12 401 floats (~2 GB) and every use so far is an amplitude ratio.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from astropy.io import fits

CACHE_VERSION = 1
CACHE_BASENAME = "{stem}_moonzodi_model_cache_v{version}.npz"
ARMS = ("near", "far", "sci")
_ARM_META = {
    "near": ("sky_near", "sky_near_ra", "sky_near_dec", "LSF_SKY_NEAR"),
    "far": ("sky_far", "sky_far_ra", "sky_far_dec", "LSF_SKY_FAR"),
    "sci": ("sci", "sci_ra", "sci_dec", "LSF_SCI"),
}
# decompose_parallel's own default when no exposure column is available.  The
# corpus META carries `exposure`, which is an exposure NUMBER, not a duration.
# MoonZodiObservation validates the SOURCE string against a closed set, so it
# must be exactly 'assumed_900s' here -- anything else raises in __post_init__,
# and with a blanket except that shows up as a cache full of NaN.
DEFAULT_EXPOSURE_SECONDS = 900.0
EXPOSURE_SOURCE = "assumed_900s"
FIT_FLUX_SCALE = 1e14

_FIELDS = ("moon_total", "zodi_total", "moon_sep_deg", "target_airmass",
           "moon_airmass", "moon_alt_deg", "sun_alt_deg", "zodi_b500")

# Worker-global state, set once per process by _init_worker.  Passing the flux
# stack through the pool would pickle gigabytes per task.
_W = {}


def cache_path(corpus_prefix, version=CACHE_VERSION):
    """Cache location: beside the corpus, so it moves with it.

    The STEM is part of the filename.  The full corpus and its every10
    subsample live in the same directory, so a stem-free name would have the
    every10 build silently overwrite the corpus cache -- and since `load`
    validates on `expnum`, the symptom would be a confusing refusal on the
    corpus rather than an obvious clobber.
    """
    prefix = Path(str(corpus_prefix))
    return prefix.parent / CACHE_BASENAME.format(stem=prefix.name,
                                                 version=version)


def _sanitised_lsf(lsf):
    """Detector LSF with unusable pixels interpolated, or None if hopeless.

    Same policy as ``decompose_parallel._sanitised_lsf_row``: the gaia1over100
    selection has 9 rows carrying a single 0.0 LSF pixel at a spectrograph arm
    join, and ``predict`` rightly refuses a non-positive FWHM.  A one-pixel
    repair cannot move a band integral measurably; a row with no usable pixel
    at all is refused rather than fitted against a fabricated LSF.
    """
    lsf = np.asarray(lsf, dtype=np.float64)
    good = np.isfinite(lsf) & (lsf > 0.0)
    if good.all():
        return lsf
    if not good.any():
        return None
    idx = np.arange(lsf.size)
    out = lsf.copy()
    out[~good] = np.interp(idx[~good], idx[good], lsf[good])
    return out


def _init_worker(stack_path, exposure_seconds):
    from sky_decomp.moon_zodi_model import MoonZodiPhysicalModel
    hdul = fits.open(str(stack_path), memmap=True)
    wave = np.asarray(hdul["WAVE"].data, dtype=np.float64)
    _W["hdul"] = hdul
    _W["wave"] = wave if wave.ndim == 1 else wave[0]
    _W["meta"] = hdul["META"].data
    _W["lsf"] = {a: hdul[_ARM_META[a][3]].section for a in ARMS}
    _W["model"] = MoonZodiPhysicalModel()
    _W["exposure_seconds"] = float(exposure_seconds)


def _run_chunk(rows):
    from sky_decomp.moon_zodi_model import MoonZodiObservation
    meta, wave, model = _W["meta"], _W["wave"], _W["model"]
    out = {f"{a}_{f}": np.full(len(rows), np.nan) for a in ARMS for f in _FIELDS}
    ok = np.zeros((len(rows), len(ARMS)), dtype=bool)
    fails = []
    for i, r in enumerate(rows):
        m = meta[int(r)]
        raw = m["date_obs"]
        date_obs = (raw.decode().strip() if isinstance(raw, bytes)
                    else str(raw).strip())
        for j, arm in enumerate(ARMS):
            role, ra_col, dec_col, _ = _ARM_META[arm]
            lsf = _sanitised_lsf(np.asarray(_W["lsf"][arm][int(r)]))
            if lsf is None:
                continue
            try:
                pred = model.predict(
                    wave, lsf,
                    MoonZodiObservation(
                        expnum=int(m["expnum"]), date_obs=date_obs, role=role,
                        target_ra_deg=float(m[ra_col]),
                        target_dec_deg=float(m[dec_col]),
                        exposure_seconds=_W["exposure_seconds"],
                        exposure_seconds_source=EXPOSURE_SOURCE),
                    physical_to_fit_flux_scale=FIT_FLUX_SCALE)
            except Exception as exc:
                # Geometry the model refuses, or a bad row.  NaN + ok=False
                # keeps the cache aligned with the corpus.  The reason is kept
                # and surfaced by build(): a blanket except here once turned a
                # one-word argument error into a silently all-NaN cache.
                if len(fails) < 5:
                    fails.append((int(r), arm, f"{type(exc).__name__}: {exc}"))
                continue
            g = pred.state.geometry
            out[f"{arm}_moon_total"][i] = float(np.nansum(pred.moon))
            out[f"{arm}_zodi_total"][i] = float(np.nansum(pred.zodi))
            out[f"{arm}_moon_sep_deg"][i] = float(g.moon_separation_deg)
            out[f"{arm}_target_airmass"][i] = float(g.target_airmass)
            out[f"{arm}_moon_airmass"][i] = float(g.moon_airmass)
            out[f"{arm}_moon_alt_deg"][i] = float(g.moon_altitude_deg)
            out[f"{arm}_sun_alt_deg"][i] = float(g.sun_altitude_deg)
            out[f"{arm}_zodi_b500"][i] = float(g.zodi_b500)
            ok[i, j] = True
    return rows, out, ok, fails


def build(corpus_prefix, n_workers=8, chunk_size=32,
          exposure_seconds=DEFAULT_EXPOSURE_SECONDS, rows=None,
          overwrite=False, verbose=True):
    """Compute the cache with a process pool and write it beside the corpus.

    ``rows`` restricts the computation (for testing); a partial cache is marked
    ``complete = False`` and ``load`` refuses it, so a smoke test can never be
    mistaken for the real thing.
    """
    import multiprocessing as mp
    import time

    stack = Path(f"{corpus_prefix}.fits")
    if not stack.exists():
        raise FileNotFoundError(f"corpus stack not found: {stack}")
    out_path = cache_path(corpus_prefix)
    if out_path.exists() and not overwrite:
        raise FileExistsError(
            f"{out_path} exists; pass overwrite=True to rebuild it")
    with fits.open(stack, memmap=True) as hdul:
        n_rows = int(hdul["FLUX_SCI"].shape[0])
        _meta = hdul["META"].data
        expnum = np.asarray(_meta["expnum"], dtype=np.int64)
        # The science pointing is stored so a consumer holding only a triplet
        # can verify the cache belongs to ITS file.  A row-count check is not
        # enough: an every10 triplet's row_index (0..1446) fits happily inside
        # a 14 469-row corpus cache while meaning entirely different spectra,
        # which would attach every row's geometry to the wrong exposure.
        sci_ra = np.asarray(_meta["sci_ra"], dtype=np.float64)
        sci_dec = np.asarray(_meta["sci_dec"], dtype=np.float64)
    all_rows = np.arange(n_rows) if rows is None else np.asarray(rows, dtype=int)
    complete = rows is None
    chunks = [all_rows[i:i + int(chunk_size)]
              for i in range(0, all_rows.size, int(chunk_size))]
    res = {f"{a}_{f}": np.full(n_rows, np.nan) for a in ARMS for f in _FIELDS}
    ok = np.zeros((n_rows, len(ARMS)), dtype=bool)
    t0 = time.perf_counter()
    if verbose:
        print(f"[moon-model-cache] {all_rows.size} rows x {len(ARMS)} arms on "
              f"{n_workers} workers ({len(chunks)} chunks of {chunk_size})",
              flush=True)
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=int(n_workers), initializer=_init_worker,
                  initargs=(str(stack), float(exposure_seconds))) as pool:
        done = 0
        fails = []
        for chunk_rows, vals, chunk_ok, chunk_fails in pool.imap_unordered(
                _run_chunk, chunks):
            for key, arr in vals.items():
                res[key][chunk_rows] = arr
            ok[chunk_rows] = chunk_ok
            if len(fails) < 5:
                fails.extend(chunk_fails[:5 - len(fails)])
            done += len(chunk_rows)
            if verbose and (done % (20 * int(chunk_size)) < int(chunk_size)):
                el = time.perf_counter() - t0
                print(f"  {done}/{all_rows.size}  {el:.0f}s elapsed, "
                      f"~{el / max(done, 1) * (all_rows.size - done):.0f}s left",
                      flush=True)
    payload = dict(res)
    payload["ok"] = ok
    payload["arms"] = np.array(ARMS)
    payload["fields"] = np.array(_FIELDS)
    payload["expnum"] = expnum
    payload["sci_ra"] = sci_ra
    payload["sci_dec"] = sci_dec
    payload["n_rows"] = np.array(n_rows)
    payload["complete"] = np.array(bool(complete))
    payload["version"] = np.array(CACHE_VERSION)
    payload["exposure_seconds"] = np.array(float(exposure_seconds))
    payload["fit_flux_scale"] = np.array(FIT_FLUX_SCALE)
    np.savez_compressed(out_path, **payload)
    # Report the usable fraction over the rows actually COMPUTED, not over the
    # whole corpus -- otherwise a `rows=` smoke test always looks like a
    # near-total failure.
    usable = ok[all_rows].all(axis=1)
    if not usable.any():
        raise RuntimeError(
            "the moon-model cache came out empty: no row produced all three "
            "arms. First failures: "
            + "; ".join(f"row {r} {a}: {m}" for r, a, m in fails[:5]))
    if verbose:
        print(f"[moon-model-cache] wrote {out_path} in "
              f"{time.perf_counter() - t0:.0f}s; all three arms usable on "
              f"{usable.mean() * 100.0:.2f}% of the {all_rows.size} rows "
              f"computed"
              + ("" if complete else "  (PARTIAL -- load() will refuse it)"))
        if fails:
            print("  first failures: "
                  + "; ".join(f"row {r} {a}: {m}" for r, a, m in fails[:3]))
    return out_path


def load(corpus_prefix, expnum=None, require_complete=True):
    """Load the cache, validating it against the corpus it claims to describe.

    ``expnum`` (the corpus META column) is checked element-wise.  That is the
    self-invalidation: a cache built for a different selection, or for a corpus
    whose rows were reordered, is refused rather than silently misaligned --
    which would attach every row's geometry to the wrong spectrum.
    """
    path = cache_path(corpus_prefix)
    if not path.exists():
        raise FileNotFoundError(
            f"no moon-model cache at {path}; build it with "
            f"`python -m mlp_predictor.moon_model_cache {corpus_prefix}`")
    z = np.load(path, allow_pickle=False)
    if int(z["version"]) != CACHE_VERSION:
        raise RuntimeError(
            f"{path} is cache version {int(z['version'])}, expected "
            f"{CACHE_VERSION}; rebuild it")
    if require_complete and not bool(z["complete"]):
        raise RuntimeError(
            f"{path} is a PARTIAL cache (built with rows=...); rebuild it "
            f"without `rows` before using it for training")
    if expnum is not None:
        want = np.asarray(expnum, dtype=np.int64)
        have = np.asarray(z["expnum"], dtype=np.int64)
        if want.shape != have.shape or not np.array_equal(want, have):
            raise RuntimeError(
                f"{path} does not match this corpus: cache has {have.size} "
                f"rows, corpus has {want.size}"
                + ("" if want.shape != have.shape else
                   f", and {int((want != have).sum())} expnum values differ"))
    return {k: z[k] for k in z.files}


def default_workers():
    """Worker count for an implicit build: 8, or fewer on a smaller machine."""
    return max(1, min(8, os.cpu_count() or 1))


def load_or_build(corpus_prefix, n_workers=None, expnum=None, verbose=True,
                  chunk_size=32):
    """Load the cache, building it first if it does not exist yet.

    Mirrors the wavelength cache's contract: a corpus without a cache simply
    grows one on first use.  ABSENCE triggers a build; a cache that exists but
    does not VALIDATE does not -- it raises instead.  That asymmetry is
    deliberate.  A validation failure means either a stale cache or a caller
    that passed the wrong prefix, and silently spending ~26 minutes of CPU to
    find out which is not a reasonable default; the message says what to run.
    """
    path = cache_path(corpus_prefix)
    if not path.exists():
        workers = default_workers() if n_workers is None else int(n_workers)
        if verbose:
            print(f"  [moon-model-cache] no cache at {path}; building it now on "
                  f"{workers} worker(s).  This is a one-off: ~26 min for a full "
                  f"corpus, ~3 min for an every10 subsample, then free.")
        build(corpus_prefix, n_workers=workers, chunk_size=chunk_size,
              overwrite=False, verbose=verbose)
    return load(corpus_prefix, expnum=expnum)


def transfer_ratio(cache, reference="near", eps=1e-30):
    """``log10(moon_sci / moon_reference)`` -- the quantity that beat the network.

    NaN where either arm is unusable or non-positive, which includes every
    moon-down row: the model's moon there is ~1e-2 in fit units against ~1e4
    with the moon up, so the ratio is numerically meaningless and the consumer
    must gate on ``moon_alt_deg`` rather than trusting a finite value.
    """
    num = np.asarray(cache["sci_moon_total"], dtype=np.float64)
    den = np.asarray(cache[f"{reference}_moon_total"], dtype=np.float64)
    good = np.isfinite(num) & np.isfinite(den) & (num > eps) & (den > eps)
    out = np.full(num.shape, np.nan)
    out[good] = np.log10(num[good] / den[good])
    return out


def _main(argv=None):
    import argparse
    p = argparse.ArgumentParser(
        description="Build the physical Moon/Zodi model cache for one corpus.")
    p.add_argument("corpus_prefix",
                   help="e.g. gaia-stars/lvmsframe_median_stack_1.2.1_gaia1over100 "
                        "(without the .fits)")
    p.add_argument("--n-workers", type=int, default=8)
    p.add_argument("--chunk-size", type=int, default=32)
    p.add_argument("--rows", type=int, default=None,
                   help="only the first N rows, for a smoke test; the result is "
                        "marked PARTIAL and load() refuses it")
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args(argv)
    build(a.corpus_prefix, n_workers=a.n_workers, chunk_size=a.chunk_size,
          rows=(None if a.rows is None else np.arange(a.rows)),
          overwrite=a.overwrite)


if __name__ == "__main__":
    _main()
