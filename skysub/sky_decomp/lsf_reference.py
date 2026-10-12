"""Read empirical LSF references without fitting or resampling their coefficients."""

from pathlib import Path
from collections.abc import Mapping
from contextlib import ExitStack
from copy import deepcopy
import re

from astropy.time import Time

import numpy as np
from astropy.io import fits

from .lsf_surface_iterative import LSFSurfaceState, _wave_fingerprint
from .lsf_spline2d import LSF_SPLINE2D_REPRESENTATION, evaluate_bspline_basis
from .result_io import STATE_CONFIG_COLUMNS, _text_value


def validate_coefficients(coefficients, tolerance=2e-7):
    """Validate the convex constraints, including the common central component."""
    value = np.asarray(coefficients, dtype=float)
    if value.ndim < 2 or not np.isfinite(value).all():
        raise ValueError("LSF coefficients must be finite matrices")
    if np.any(value < -tolerance) or not np.allclose(value.sum(axis=-2), 1, atol=tolerance, rtol=0):
        raise ValueError("LSF coefficients must be nonnegative and unit-sum")
    center = value.shape[-2] // 2
    if (np.any(np.diff(value[..., :center + 1, :], axis=-2) < -tolerance)
            or np.any(np.diff(value[..., center:, :], axis=-2) > tolerance)):
        raise ValueError("LSF coefficients must obey the central-component shape constraints")


def load_reference_lsf(path: str | Path | Mapping, fiberid: int, *, allow_epoch_pooling=False):
    """Return a B/R/Z LSFSurfaceState for one fiber in one reference FITS file.

    A path uses one file for B/R/Z. A mapping {"B": path_b, "R": path_r,
    "Z": path_z} supports independent camera epochs or sparse coverage.
    Each selected channel must have exactly one reference. Survey-wide
    descriptive pooling requires opt-in.
    """
    channels = ("B", "R", "Z")
    paths = dict(path) if isinstance(path, Mapping) else dict.fromkeys(channels, path)
    if set(paths) != set(channels):
        raise ValueError("Channel paths must identify B, R and Z")
    with ExitStack() as stack:
        opened = {str(p): stack.enter_context(fits.open(p, memmap=True)) for p in set(map(str, paths.values()))}
        coefficients, knots, degrees, bounds, metrics = {}, {}, {}, {}, {}
        config, reference_fingerprint, offsets = {}, None, None
        for channel in channels:
            source = paths[channel]
            hdus = opened[str(source)]
            if hdus[0].header.get("REFSCHEM") != 1:
                raise ValueError("Unsupported LSF reference schema")
            if not hdus[0].header["APPLYOK"] and not allow_epoch_pooling:
                raise ValueError("Survey pooling spans epochs; use an epoch-specific reference or explicitly allow pooling")
            table = hdus["REFERENCES"].data
            chosen = np.flatnonzero(table["fiberid"] == fiberid)
            config_table = hdus["LSF_CONFIG"].data
            wave = np.asarray(hdus["WAVE"].data, dtype=float)
            rows = [i for i in chosen if _text_value(table[i]["channel"]) == channel]
            if len(rows) != 1:
                raise ValueError(f"Expected one {channel} reference for fiber {fiberid}; found {len(rows)}")
            i = rows[0]
            configs = [r for r in config_table if _text_value(r["channel"]) == channel]
            if len(configs) != 1:
                raise ValueError(f"Expected one configuration for channel {channel}")
            cf = configs[0]
            if _text_value(cf["legacy_kernel_representation"]) != LSF_SPLINE2D_REPRESENTATION:
                raise ValueError("Reference must contain continuous M-spline LSF densities")
            coefficient = np.asarray(hdus["LSF_COEF"].data[i], dtype=float).copy()
            knot = np.asarray(hdus["LSF_KNOTS"].data[i], dtype=float).copy()
            degree = int(cf["degree"])
            if coefficient.shape != (len(cf["tap_offsets"]), int(cf["n_basis"])):
                raise ValueError("LSF coefficient dimensions disagree with configuration")
            if len(knot) != int(cf["n_knots"]) or len(knot) - degree - 1 != coefficient.shape[1]:
                raise ValueError("LSF knot dimensions disagree with configuration")
            basis = evaluate_bspline_basis(np.linspace(knot[degree], knot[-degree - 1], 31), knot, degree)
            if not np.isfinite(basis).all() or not np.allclose(basis.sum(axis=1), 1):
                raise ValueError("LSF wavelength basis must form a partition of unity")
            validate_coefficients(coefficient)
            coefficients[channel], knots[channel], degrees[channel] = coefficient, knot, degree
            bounds[channel] = tuple(None if not np.isfinite(cf[k]) else float(cf[k]) for k in ("lower", "upper"))
            local_config = {key: converter(cf[column]) for key, column, converter in STATE_CONFIG_COLUMNS}
            local_config.update(offset_half_width_angstrom=float(cf["offset_half_width_angstrom"]),
                                offset_roughness_fraction=float(cf["offset_roughness_fraction"]))
            if config and local_config != config:
                raise ValueError("Channel configurations disagree")
            config = local_config
            fingerprint = _wave_fingerprint(wave)
            expected = (int(cf["wave_n"]), float(cf["wave_min"]), float(cf["wave_max"]), _text_value(cf["wave_sha256"]))
            if fingerprint != expected:
                raise ValueError("Reference wavelength-grid fingerprint is inconsistent")
            if reference_fingerprint is not None and reference_fingerprint != fingerprint:
                raise ValueError("Selected channels use different native wavelength grids")
            reference_fingerprint = fingerprint
            channel_offsets = np.asarray(cf["tap_offsets"], dtype=int)
            if offsets is not None and not np.array_equal(offsets, channel_offsets):
                raise ValueError("Selected channels use different displacement components")
            offsets = channel_offsets.copy()
            metrics[channel] = dict(status="Solved", model="dispersion_weighted_reference",
                                    n_exposures=int(table[i]["n_exposures"]), n_nights=int(table[i]["n_nights"]),
                                    calibration_epoch=int(table[i]["calibration_epoch"]),
                                    reference_file=str(Path(source).resolve()), reference_id=int(table[i]["reference_id"]))
        return LSFSurfaceState(coefficients=coefficients, knot_vectors=knots, degrees=degrees,
                               channel_bounds=bounds, tap_offsets=offsets,
                               config=config, metrics=metrics, requested_cycles=0, completed_cycles=0,
                               wave_n=fingerprint[0], wave_min=fingerprint[1], wave_max=fingerprint[2],
                               wave_sha256=fingerprint[3], fit_status="empirical_reference",
                               knot_strategy="inherited_fixed_basis", legacy_kernel_representation=LSF_SPLINE2D_REPRESENTATION)


class SurveyReferenceSet:
    """Every per-fiber reference of one file, for means over any set of fibers.

    For spectra stacked over many fibers (the ML corpus rows).  All fibers of a
    channel share one knot vector, and a mean of convex LSF coefficients is
    itself convex, so every mean is a valid B/R/Z state.  Only fibers with B,
    R and Z references enter a mean.
    """

    def __init__(self, path: str | Path, *, allow_epoch_pooling=False):
        self.path = Path(path)
        with fits.open(self.path, memmap=False) as hdus:
            table = hdus["REFERENCES"].data
            telescopes = np.char.strip(np.asarray(table["telescope"]).astype(str))
            channels = np.char.strip(np.asarray(table["channel"]).astype(str))
            fiberids = np.asarray(table["fiberid"], dtype=np.int64)
            n_exposures = np.asarray(table["n_exposures"], dtype=np.int64)
            coefficients = np.asarray(hdus["LSF_COEF"].data, dtype=float)
            knots = np.asarray(hdus["LSF_KNOTS"].data, dtype=float)
        self._channel = {}
        complete = None
        for channel in ("B", "R", "Z"):
            rows = np.flatnonzero(channels == channel)
            if rows.size == 0:
                raise ValueError(f"No {channel} references in {path}")
            rows = rows[np.argsort(fiberids[rows], kind="stable")]
            if np.any(np.diff(fiberids[rows]) == 0):
                raise ValueError(f"Duplicate {channel} fiber references in {path}")
            if not np.all(knots[rows] == knots[rows[0]]):
                raise ValueError(f"{channel} references use different knot vectors")
            self._channel[channel] = (fiberids[rows], coefficients[rows], n_exposures[rows])
            complete = set(fiberids[rows]) if complete is None else complete & set(fiberids[rows])
        self.fiberids = np.array(sorted(complete), dtype=np.int64)
        if self.fiberids.size == 0:
            raise ValueError(f"No fiber has B, R and Z references in {path}")
        by_fiber = dict(zip(fiberids, telescopes))
        self.telescope = np.array([by_fiber[f] for f in self.fiberids])
        # One member supplies (and validates) knots, configuration and grid.
        self._template = load_reference_lsf(self.path, int(self.fiberids[0]),
                                            allow_epoch_pooling=allow_epoch_pooling)

    def mean(self, fiberids, *, model="fiber_mean_reference", label=""):
        """Equal-weight mean over those of ``fiberids`` that have references."""
        requested = np.unique(np.asarray(fiberids, dtype=np.int64))
        used = np.intersect1d(requested, self.fiberids)
        if used.size == 0:
            raise ValueError(f"None of the {requested.size} fibers has a reference in {self.path}")
        state = deepcopy(self._template)
        for channel, (ids, coefficients, n_exposures) in self._channel.items():
            rows = np.searchsorted(ids, used)
            mean = coefficients[rows].mean(axis=0)
            validate_coefficients(mean)
            state.coefficients[channel] = mean
            metric = state.metrics[channel]
            metric.pop("n_nights", None)   # the template fiber's, not the mean's
            metric.update(model=model, label=label, n_fibers=int(used.size),
                          n_fibers_requested=int(requested.size),
                          n_exposures=int(n_exposures[rows].sum()), reference_id=-1)
        state.fit_status = f"empirical_reference_{model}"
        return state

    def telescope_mean(self, telescope):
        """Mean over every referenced fiber of one telescope (Sci, SkyE or SkyW)."""
        return self.mean(self.fiberids[self.telescope == telescope],
                         model="telescope_mean", label=telescope)


def telescope_mean_reference_lsf(path: str | Path, telescope: str, *, allow_epoch_pooling=False):
    """Equal-weight mean over every fiber reference of one telescope in one file."""
    return SurveyReferenceSet(path, allow_epoch_pooling=allow_epoch_pooling).telescope_mean(telescope)


def select_reference_lsf(directory, fiberid, camera_epochs, observation_mjd):
    """Choose the nearest measured lunation independently within each camera epoch."""
    if set(camera_epochs) != {"B", "R", "Z"} or not np.isfinite(observation_mjd):
        raise ValueError("Provide B/R/Z calibration epochs and a finite observation MJD")
    candidates = {c: [] for c in camera_epochs}
    for path in sorted((Path(directory) / "lunation").glob("*.fits")):
        with fits.open(path, memmap=True) as hdus:
            for row in hdus["REFERENCES"].data:
                channel = _text_value(row["channel"])
                if (int(row["fiberid"]) == fiberid and channel in candidates
                        and int(row["calibration_epoch"]) == int(camera_epochs[channel])):
                    date = float(row["weighted_mjd"])
                    if not np.isfinite(date):
                        raise ValueError(f"Nonfinite reference timestamp: {path}")
                    candidates[channel].append((abs(date-observation_mjd), str(path)))
    paths = {}
    for channel, choices in candidates.items():
        if not choices:
            raise ValueError(f"No {channel} LSF for fiber {fiberid}, camera epoch {camera_epochs[channel]}")
        paths[channel] = min(choices)[1]
    return load_reference_lsf(paths, fiberid)


def load_lsf_from_sframe(sframe, fiberid, directory=None, *, scope="lunation"):
    """Load a nearest epoch–lunation or survey mean for the target physical fiber."""
    if scope not in ("lunation", "survey"):
        raise ValueError("LSF scope must be lunation or survey")
    directory = Path(directory) if directory is not None else Path(__file__).parent
    with fits.open(sframe, memmap=True) as hdus:
        header = hdus[0].header
        table = hdus["SLITMAP"].data
        rows = table[table["fiberid"] == fiberid]
        if len(rows) != 1:
            raise ValueError(f"Expected one physical fiber {fiberid} in the target SLITMAP")
        if scope == "lunation":
            sp = int(rows[0]["spectrographid"])
            comments = "".join(header.get("COMMENT", []))
            epochs = {}
            for epoch, camera in re.findall(r"/calib/(\d+)/lvm-mlsf-([brz][123])\.fits", comments):
                if int(camera[1]) == sp:
                    channel = camera[0].upper()
                    if channel in epochs and epochs[channel] != int(epoch):
                        raise ValueError(f"Conflicting calibration provenance for camera {camera}")
                    epochs[channel] = int(epoch)
            if set(epochs) != {"B", "R", "Z"}:
                raise ValueError("Target SFrame must identify the three camera LSF calibration files")
            if "OBSTIME" not in header:
                raise ValueError("Target SFrame must identify its UTC observation timestamp")
            date = Time(header["OBSTIME"], scale="utc").mjd
        fingerprint = _wave_fingerprint(np.asarray(hdus["WAVE"].data, dtype=float))
    state = (load_reference_lsf(directory / "survey/lsf_survey.fits", fiberid, allow_epoch_pooling=True)
             if scope == "survey" else select_reference_lsf(directory, fiberid, epochs, date))
    if fingerprint != (state.wave_n, state.wave_min, state.wave_max, state.wave_sha256):
        raise ValueError("Target wavelength grid differs from the stored native calibration grid")
    return state


def reference_covariance(path: str | Path, reference_id: int, *, kind="scatter"):
    """Reconstruct a full 44x44 covariance from stored factors in C order.

    ``scatter`` is empirical exposure scatter. ``conditional_mean`` is a
    fixed-weight, night-cluster estimate, not a calibrated prediction error.
    """
    path = Path(path)
    with fits.open(path, memmap=True) as probe:
        expected_build = probe[0].header.get("BUILDFP")
        archive = probe[0].header.get("AUDIT")
        if "SCATTER_FACTORS" not in probe and archive:
            path = path.parent / archive
    with fits.open(path, memmap=True) as hdus:
        if hdus[0].header.get("BUILDFP") != expected_build:
            raise ValueError("Reference and calculation archive fingerprints disagree")
        refs = hdus["REFERENCES"].data
        ids = np.flatnonzero(refs["reference_id"] == reference_id)
        if len(ids) != 1:
            raise ValueError("Expected one reference_id")
        row = refs[ids[0]]
        if kind == "scatter":
            mapping, extension = "CONTRIBUTORS", "SCATTER_FACTORS"
            denominator = float(row["scatter_denominator"])
        elif kind == "conditional_mean":
            mapping, extension = "NIGHTS", "NIGHT_FACTORS"
            denominator = (row["n_nights"] - 1) / row["n_nights"]
        else:
            raise ValueError("kind must be scatter or conditional_mean")
        factors = np.asarray(hdus[extension].data[hdus[mapping].data["reference_id"] == reference_id], dtype=float)
        if denominator <= 0:
            return np.full((factors.shape[1], factors.shape[1]), np.nan)
        return factors.T @ factors / denominator
