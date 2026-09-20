from importlib.util import find_spec
from pathlib import Path

from astropy.io import fits
from astropy.table import Table
from typer.testing import CliRunner

from lvm_medians.cli import app
from lvm_medians.stack import read_manifest


def test_scan_can_write_a_small_sample(tmp_path: Path) -> None:
    for expnum in range(1, 5):
        (tmp_path / f"lvmSFrame-{expnum:08d}.fits").touch()

    work_dir = tmp_path / "work"
    result = CliRunner().invoke(
        app,
        [
            "--work-dir",
            str(work_dir),
            "--no-progress",
            "scan",
            "--sframes-root",
            str(tmp_path),
            "--every-nth",
            "2",
            "--limit",
            "2",
        ],
    )

    assert result.exit_code == 0, result.output
    selected = [path.name for _, path in read_manifest(work_dir / "sframes.txt")]
    assert selected == ["lvmSFrame-00000001.fits", "lvmSFrame-00000003.fits"]

    status = CliRunner().invoke(app, ["--work-dir", str(work_dir), "status"])
    assert status.exit_code == 0, status.output
    assert "SFrames in list: 2" in status.output
    assert "awaiting download: 2" in status.output
    assert "manifest" not in status.output.lower()


def test_combine_gaia_writes_one_fits_table(tmp_path: Path) -> None:
    cache_dir = tmp_path / "gaia"
    fibers_dir = cache_dir / "fibers"
    fibers_dir.mkdir(parents=True)
    inputs = []
    for expnum in (1, 2):
        sframe = tmp_path / f"lvmSFrame-{expnum:08d}.fits"
        sframe.touch()
        inputs.append(sframe)
        Table(
            {
                "expnum": [expnum, expnum],
                "fiberid": [1, 2],
                "gaia_g_flux": [0.1, 0.2],
            }
        ).write(fibers_dir / f"lvmGAIA-fibers-{expnum:08d}.fits")

    manifest = tmp_path / "sframes.txt"
    manifest.write_text("".join(f"{path}\n" for path in inputs), encoding="utf-8")
    output = tmp_path / "all-fibers.fits"
    result = CliRunner().invoke(
        app,
        [
            "--work-dir",
            str(tmp_path / "work"),
            "--no-progress",
            "combine-gaia",
            "--sframe-list",
            str(manifest),
            "--cache-dir",
            str(cache_dir),
            "--table",
            "fibers",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    with fits.open(output, checksum=True) as hdul:
        assert len(hdul["FIBERS"].data) == 4
        assert hdul[0].header["NEXP"] == 2
        assert hdul[0].header["COMPLETE"]

    parquet = tmp_path / "all-fibers.parquet"
    parquet_result = CliRunner().invoke(
        app,
        [
            "--work-dir",
            str(tmp_path / "work"),
            "--no-progress",
            "combine-gaia",
            "--sframe-list",
            str(manifest),
            "--cache-dir",
            str(cache_dir),
            "--table",
            "fibers",
            "--output",
            str(parquet),
        ],
    )
    if find_spec("pyarrow") is None:
        assert parquet_result.exit_code == 1
        assert "Parquet output requires" in parquet_result.output
    else:
        assert parquet_result.exit_code == 0, parquet_result.output
        assert len(Table.read(parquet)) == 4


def test_help_is_concise() -> None:
    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "--install-completion" not in result.output
    assert "--show-completion" not in result.output
    assert "10 MiB" not in result.output
    assert "older copies" not in result.output
