"""Deterministic exposure selection, DRP paths and verified product reuse."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path, PurePosixPath
import hashlib
import json
import math
import os
import re

from .options import DEFAULTS, validate_model


@dataclass
class Exposure:
    input_path: Path
    expnum: int
    mjd: int
    tileid: int
    version: str
    tilegroup: str
    output: Path | None = None
    seed: Path | None = None
    state_mode: str = ""
    fingerprint: str = ""
    input_sha256: str = ""
    dependencies: dict | None = None
    status: str = "queued"
    error: str = ""
    log: str = ""
    pid: int | None = None
    elapsed_seconds: float = 0.0

    @property
    def identity(self):
        return self.version, self.mjd, self.expnum, self.tileid

    def record(self):
        return asdict(self)


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def read_metadata(hdus):
    return json.loads(bytes(hdus["SKYML_META"].data).decode("utf-8"))


def _text(value):
    return (value.decode() if isinstance(value, bytes) else str(value)).strip()


def validate_version(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value) or value in (".", ".."):
        raise ValueError("DRP version must be one directory name, e.g. 1.3.2")
    return value


def paths(sas_base_dir=None, drp_version="1.3.2", cframe_root=None, drpall=None):
    validate_version(drp_version)
    sas = Path(sas_base_dir).expanduser() if sas_base_dir else None
    if cframe_root:
        root = Path(cframe_root).expanduser()
    elif sas:
        root = sas / "sdsswork/lvm/spectro/redux" / drp_version
    elif os.environ.get("LVM_SPECTRO_REDUX"):
        root = Path(os.environ["LVM_SPECTRO_REDUX"]).expanduser() / drp_version
    elif os.environ.get("SAS_BASE_DIR"):
        root = Path(os.environ["SAS_BASE_DIR"]).expanduser() / "sdsswork/lvm/spectro/redux" / drp_version
    elif drpall:
        raise ValueError("A relocated --drpall needs --cframe-root or a SAS environment")
    else:
        raise ValueError("Set SAS_BASE_DIR, --sas-base-dir or --cframe-root for catalogue selection")
    root = root.resolve()
    catalogue = Path(drpall).expanduser().resolve() if drpall else root / f"drpall-{drp_version}.fits"
    return root, catalogue


def parse_values(values):
    """Represent inclusive/open integer ranges without expanding large ranges."""
    intervals = []
    for value in values or []:
        for token in value.split(","):
            token = token.strip()
            if not re.fullmatch(r"\d+|\d*:\d*", token) or token == ":":
                raise ValueError(f"Invalid selector {token!r}; use 11017,11020:11025 or 11112:")
            if ":" in token:
                lo, hi = token.split(":")
                lo, hi = int(lo) if lo else None, int(hi) if hi else None
            else:
                lo = hi = int(token)
            if lo is not None and hi is not None and hi < lo:
                raise ValueError(f"Range {token!r} must be increasing")
            intervals.append((lo, hi))
    return intervals


def includes(value, intervals):
    return any((lo is None or value >= lo) and (hi is None or value <= hi) for lo, hi in intervals)


def exposure_from_file(path):
    from astropy.io import fits

    path = Path(path).expanduser().resolve()
    if not re.fullmatch(r"lvmCFrame-\d+\.fits(?:\.gz)?", path.name):
        raise ValueError(f"Input must be an lvmCFrame FITS, got {path}")
    header = fits.getheader(path)
    expnum = int(header["EXPOSURE"])
    if int(path.name.split("-")[-1].split(".")[0]) != expnum:
        raise ValueError(f"Input basename and EXPOSURE disagree: {path}")
    version = validate_version(str(header["DRPVER"]).strip())
    tileid = int(header.get("TILE_ID", header.get("TILEID", -1)))
    mjd = int(header["MJD"])
    if min(expnum, tileid, mjd) < 0:
        raise ValueError(f"Invalid exposure/tile/MJD identity in {path}")
    return Exposure(path, expnum, mjd, tileid, version, f"{tileid // 1000:04d}XX")


def _catalogue_exposures(catalogue, root, version, selectors, exptime, fluxcal, all_exposures):
    import numpy as np
    from astropy.table import Table

    table = Table.read(catalogue)
    columns = {name.lower(): name for name in table.colnames}
    for name in ("expnum", "mjd", "tileid"):
        if name not in columns:
            raise ValueError(f"Catalogue requires column {name}: {catalogue}")
    mask = np.ones(len(table), dtype=bool)
    intervals = {name: parse_values(value) for name, value in selectors.items()}
    for name, ranges in intervals.items():
        if ranges:
            mask &= np.array([includes(int(v), ranges) for v in table[columns[name]]])
    if exptime is not None:
        if not math.isfinite(exptime) or exptime < 0:
            raise ValueError("--exptime must be finite and >= 0")
        if "exptime" not in columns:
            raise ValueError("--exptime needs the catalogue exptime column")
        mask &= np.asarray(table[columns["exptime"]], float) == exptime
    if fluxcal is not None:
        if "fluxcal" not in columns:
            raise ValueError("--catalog-fluxcal needs the catalogue fluxcal column")
        mask &= np.array([_text(v) == fluxcal for v in table[columns["fluxcal"]]])
    # Explicit exposure requests must exist even when another selector would filter them.
    if intervals["expnum"]:
        available = set(map(int, table[columns["expnum"]]))
        for lo, hi in intervals["expnum"]:
            if lo is None or hi is None:
                continue
            count = sum(lo <= v <= hi for v in available)
            if count != hi - lo + 1:
                raise ValueError(f"Requested exposure(s) in {lo}:{hi} are absent from {catalogue}")
    result = []
    for row in table[mask]:
        expnum, mjd, tileid = (int(row[columns[n]]) for n in ("expnum", "mjd", "tileid"))
        row_version = _text(row[columns["drpver"]]) if "drpver" in columns else version
        if row_version != version:
            raise ValueError(f"Catalogue DRP version {row_version} differs from requested {version}")
        group = _text(row[columns["tilegrp"]]) if "tilegrp" in columns else f"{tileid // 1000:04d}XX"
        if not re.fullmatch(r"\d+XX", group) or group != f"{tileid // 1000:04d}XX":
            raise ValueError(f"Invalid tile group {group!r} for tile {tileid}")
        expected = (group, str(tileid), str(mjd))
        if "location" in columns:
            location = PurePosixPath(_text(row[columns["location"]]))
            if ".." in location.parts:
                raise ValueError(f"Unsafe catalogue location for exposure {expnum}")
            if location.name and location.parent.parts[-3:] != expected:
                raise ValueError(f"Catalogue location disagrees with tile/MJD for exposure {expnum}")
            if location.name:
                match = re.fullmatch(r"lvm[CS]Frame-(\d+)\.fits(?:\.gz)?", location.name)
                if not match or int(match.group(1)) != expnum:
                    raise ValueError(f"Catalogue location disagrees with exposure {expnum}")
                if len(location.parts) >= 5 and location.parts[-5] != version:
                    raise ValueError(f"Catalogue location has a different DRP version: {location}")
        path = root.joinpath(*expected, f"lvmCFrame-{expnum:08d}.fits")
        if not path.exists() and path.with_suffix(".fits.gz").exists():
            path = path.with_suffix(".fits.gz")
        result.append(Exposure(path.resolve(), expnum, mjd, tileid, version, group))
    return result


def select_exposures(*, inputs=None, input_list=None, drpall=None, sas_base_dir=None,
                     cframe_root=None, drp_version="1.3.2", expnum=None, mjd=None,
                     tileid=None, exptime=None, catalog_fluxcal=None, all_exposures=False,
                     max_exposures=None, missing_input="error"):
    explicit = bool(inputs) or input_list is not None
    selectors = dict(expnum=expnum, mjd=mjd, tileid=tileid)
    filters = any(selectors.values()) or exptime is not None or catalog_fluxcal is not None
    if inputs and input_list:
        raise ValueError("Choose positional INPUTS or --input-list, not both")
    if explicit and (filters or all_exposures or drpall is not None):
        raise ValueError("Choose explicit paths or catalogue selection, not both")
    if all_exposures and filters:
        raise ValueError("Use --all-exposures instead of catalogue filters")
    if not explicit and not (filters or all_exposures):
        raise ValueError("Select --expnum/--mjd/--tileid/other filters, --all-exposures, or explicit INPUTS")
    if max_exposures is not None and max_exposures < 1:
        raise ValueError("--max-exposures must be >= 1")
    context = {}
    if explicit:
        paths_in = list(inputs or [])
        if input_list is not None:
            manifest = Path(input_list).expanduser().resolve()
            for line in manifest.read_text().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    p = Path(line).expanduser()
                    paths_in.append(p if p.is_absolute() else manifest.parent / p)
        exposures = []
        for path in paths_in:
            path = Path(path).expanduser().resolve()
            if not path.is_file() and missing_input == "skip":
                # Unknown header identities use -1; retain each skipped path in run records.
                match = re.fullmatch(r"lvmCFrame-(\d+)\.fits(?:\.gz)?", path.name)
                exposures.append(Exposure(path, int(match[1]) if match else -1, -1, -1,
                                          drp_version, "", status="skipped_input",
                                          error=f"Local CFrame unavailable: {path}"))
            else:
                exposures.append(exposure_from_file(path))
    else:
        root, catalogue = paths(sas_base_dir, drp_version, cframe_root, drpall)
        context = dict(cframe_root=str(root), drpall=str(catalogue))
        exposures = _catalogue_exposures(catalogue, root, drp_version, selectors,
                                         exptime, catalog_fluxcal, all_exposures)
    unique = {}
    for item in exposures:
        key = ("skipped", str(item.input_path)) if item.status == "skipped_input" else item.identity
        previous = unique.get(key)
        if previous is not None and previous.input_path != item.input_path:
            raise ValueError(f"Conflicting paths for exposure identity {item.identity}")
        unique[key] = item
    selected = sorted(unique.values(), key=lambda v: (*v.identity, str(v.input_path)))
    if not selected:
        raise ValueError("The selection contains no exposures")
    return selected[:max_exposures] if max_exposures else selected, context


def scientific_implementation():
    base = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    for path in sorted(base.rglob("*.py")):
        if "tests" in path.parts or "__pycache__" in path.parts:
            continue
        if path.parent.name == "sky_subtraction_block" and path.name in (
            "cli.py", "batch.py", "planning.py", "resources.py", "options.py", "_worker.py", "__main__.py",
        ):
            continue
        digest.update(str(path.relative_to(base)).encode())
        digest.update(path.read_bytes())
    manifest = base / "sky_decomp/data/lsf/manifest.json"
    digest.update(manifest.read_bytes())
    return digest.hexdigest()


def _check_state(seed, item, input_hash):
    import numpy as np
    from astropy.io import fits

    required = {"MED_COEF", "MED_LABELS", "ML_COEF", "MED_LSF_COEF", "RENDER_META",
                "LSF_FIBERS", "SLITMAP", "WAVE", "SKYML_META"}
    with fits.open(seed, memmap=True) as saved, fits.open(item.input_path, memmap=True) as source:
        if not required.issubset({h.name for h in saved}):
            raise ValueError(f"Saved state needs compact SkyML median/ML extensions: {seed}")
        if read_metadata(saved)["input_sha256"] != input_hash:
            raise ValueError(f"Saved state belongs to a different CFrame: {seed}")
        if (not np.array_equal(saved["WAVE"].data, source["WAVE"].data)
            or not np.array_equal(saved["SLITMAP"].data["fiberid"], source["SLITMAP"].data["fiberid"])):
            raise ValueError(f"Saved-state wavelength/fiber identities differ: {seed}")


def verify_product(item):
    from astropy.io import fits

    from ..sky_subtraction.product import SHARED_ORDER

    optional = {"GAIA_SOURCES", "FIT_LSF_HISTORY"}
    required = (set(SHARED_ORDER) - optional) | {
        "BLOCK_LSF_COEF", "BLOCK_LSF_META", "BLOCK_SLIT_KNOTS", "BLOCK_WAVE_KNOTS",
        "SKY_FIBER_COEF", "SKYML_META",
    }
    with fits.open(item.output, memmap=True) as hdus:
        if not required.issubset({h.name for h in hdus}) or hdus[0].header.get("SKYMETH") != "skyml-block":
            raise ValueError(f"Existing destination is not a complete block product: {item.output}")
        if read_metadata(hdus).get("cli_run", {}).get("fingerprint") != item.fingerprint:
            raise ValueError(f"Existing product has a different/absent CLI scientific fingerprint: {item.output}")
        shape = hdus["FLUX"].shape
        if any(hdus[n].shape != shape for n in ("IVAR", "MASK", "LSF", "SKY")):
            raise ValueError(f"Existing spectral array dimensions disagree: {item.output}")
        if any(h.verify_checksum() != 1 or h.verify_datasum() != 1 for h in hdus):
            raise ValueError(f"Existing FITS checksum validation failed: {item.output}")


def prepare_exposure(item, model, *, seed=None, implementation=None, on_existing="error"):
    """Resolve actual state mode and fingerprint without loading numerical models."""
    from ..sky_subtraction.settings import SETTINGS

    validate_model(model)
    if item.output == item.input_path or not item.output.name.startswith("lvmSFrame-skyml-block-"):
        raise ValueError("Output must use lvmSFrame-skyml-block-* and protect the input")
    if item.output.exists() and on_existing == "error":
        raise FileExistsError(f"Output exists: {item.output}; choose a new directory or --on-existing skip")
    if not item.input_path.is_file():
        raise FileNotFoundError(item.input_path)
    hash_cache = {}
    def hashed(path):
        path = Path(path).expanduser().resolve()
        if path not in hash_cache:
            hash_cache[path] = file_sha256(path)
        return hash_cache[path]

    item.input_sha256 = hashed(item.input_path)
    adjacent = item.output.parent / ("lvmSFrame-skyml-preparation-" +
                                      item.input_path.name.removeprefix("lvmCFrame-").removesuffix(".gz"))
    item.seed = Path(seed).expanduser().resolve() if seed else (adjacent if adjacent.exists() else None)
    item.state_mode = "reuse explicitly selected state" if seed else (
        "reuse adjacent saved state" if item.seed else "compute state in memory")
    dependencies = {str(item.input_path): item.input_sha256}
    state_hash = None
    ensemble_hash = None
    same_computed_snapshot = False
    if item.seed is not None:
        if item.seed == item.output:
            raise ValueError("Saved state and output must differ")
        _check_state(item.seed, item, item.input_sha256)
        state_hash = hashed(item.seed)
        dependencies[str(item.seed)] = state_hash
        if seed is None and item.output.exists() and on_existing == "skip" and model.get("ensemble"):
            from astropy.io import fits
            with fits.open(item.output, memmap=True) as previous:
                provenance = read_metadata(previous)
            same_computed_snapshot = (
                provenance.get("preparation_sha256") == state_hash
                and provenance.get("cli_run", {}).get("state_mode") == "compute state in memory"
            )
            if same_computed_snapshot:
                # --keep-preparation must not prevent resuming its own completed calculation.
                state_hash = None
                ensemble = Path(model["ensemble"]).expanduser().resolve()
                ensemble_hash = hashed(ensemble)
                dependencies[str(ensemble)] = ensemble_hash
    else:
        if not model.get("ensemble"):
            raise ValueError(f"--ml-model-file (alias --ensemble) is required for exposure {item.expnum} computing state anew")
        ensemble = Path(model["ensemble"]).expanduser().resolve()
        ensemble_hash = hashed(ensemble)
        dependencies[str(ensemble)] = ensemble_hash
    gaia_hash = None
    if (item.seed is None or same_computed_snapshot) and model.get("gaia_dir"):
        gaia = Path(model["gaia_dir"]).expanduser().resolve() / f"lvmGAIA-fibers-{item.expnum:08d}.fits"
        if gaia.is_file():
            gaia_hash = hashed(gaia)
            dependencies[str(gaia)] = gaia_hash
    scientific = {name: model[name] for name in DEFAULTS if name != "keep_preparation"}
    fingerprint = dict(input_sha256=item.input_sha256, state_sha256=state_hash,
                       ensemble_sha256=ensemble_hash, gaia_cache_sha256=gaia_hash,
                       block_settings=scientific, settings=SETTINGS.asdict(),
                       implementation_sha256=implementation or scientific_implementation())
    item.fingerprint = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()
    item.dependencies = dependencies
    if item.output.exists():
        verify_product(item)
        item.status = "reused"
    return item


def make_plan(exposures, model, *, output_dir, flat_output=False, reuse_state_dir=None,
              on_existing="error", missing_input="error"):
    output_dir = Path(output_dir).expanduser().resolve()
    implementation = scientific_implementation()
    destinations = set()
    for item in exposures:
        if item.status == "skipped_input":
            continue
        relative = Path() if flat_output else Path(item.version, item.tilegroup, str(item.tileid), str(item.mjd))
        item.output = output_dir / relative / f"lvmSFrame-skyml-block-{item.expnum:08d}.fits"
        if item.output in destinations:
            raise ValueError(f"Multiple exposures map to the same output: {item.output}")
        destinations.add(item.output)
    for item in exposures:
        if item.status == "skipped_input":
            continue
        if not item.input_path.is_file():
            if missing_input == "skip":
                item.status, item.error = "skipped_input", f"Local CFrame unavailable: {item.input_path}"
                continue
            raise FileNotFoundError(item.input_path)
        actual = exposure_from_file(item.input_path)
        if actual.identity != item.identity:
            raise ValueError(f"Catalogue and CFrame identity disagree: {item.input_path}")
        seed = None
        if reuse_state_dir:
            relative = Path() if flat_output else Path(item.version, item.tilegroup, str(item.tileid), str(item.mjd))
            seed = Path(reuse_state_dir).expanduser().resolve() / relative / f"lvmSFrame-skyml-preparation-{item.expnum:08d}.fits"
        prepare_exposure(item, model, seed=seed, implementation=implementation, on_existing=on_existing)
    return exposures
