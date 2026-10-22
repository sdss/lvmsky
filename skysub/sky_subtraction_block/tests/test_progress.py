"""Batch display counters, honest ETA and terminal/file logging separation."""
from io import StringIO
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from rich.console import Console

from skysub.sky_subtraction_block import batch
from skysub.sky_subtraction_block.batch import _progress_panel, _timing_statistics
from skysub.sky_subtraction_block.planning import Exposure


def render(statuses, elapsed=60, width=120):
    entries = [Exposure(None, 29088 + i, 60623, 11111, "1.3.2", "0011XX", status=status)
               for i, status in enumerate(statuses)]
    stream = StringIO()
    Console(file=stream, width=width, color_system=None).print(_progress_panel(entries, elapsed))
    return stream.getvalue()


def test_panel_separates_success_errors_and_active_work():
    text = render(["complete", "error", "running", "queued"])
    assert "2/4 processed" in text
    assert "OK 1 · Failed 1 · Running 1 · Queued 1" in text
    assert "Elapsed 00:01:00" in text and "1.0 successful exp/min" in text
    assert "ETA ~00:02:00" in text
    assert "Running: 29090" in text


def test_fast_failures_do_not_create_a_success_rate_or_eta():
    text = render(["error", "error", "running", "queued"], elapsed=1)
    assert "2/4 processed" in text and "OK 0 · Failed 2" in text
    assert "Rate 0.0 successful exp/min · ETA --" in text


def test_reuse_and_skips_are_processed_but_not_new_successes():
    text = render(["complete", "reused", "skipped_input", "error"])
    assert "4/4 processed" in text and "Reused 1 · Skipped 1" in text
    assert "Rate 1.0 successful exp/min · ETA 00:00:00" in text
    assert "Running: none" in text


def test_unstarted_and_interrupted_work_is_not_reported_as_processed():
    text = render(["complete", "interrupted", "not_started"])
    assert "1/3 processed" in text
    assert "Not started 1 · Interrupted 1" in text


def test_panel_fits_a_standard_terminal_without_dropping_counts():
    text = render(["complete", "error", "running", "queued"], width=80)
    assert "2/4 processed" in text and "Failed 1" in text
    assert max(map(len, text.splitlines())) <= 80


def test_success_duration_statistics_exclude_failures_and_reuse():
    entries = [Exposure(None, 100 + i, 60623, 11111, "1.3.2", "0011XX",
                        status="complete", elapsed_seconds=value)
               for i, value in enumerate((60, 120, 180, 240))]
    entries.extend(Exposure(None, 200 + i, 60623, 11111, "1.3.2", "0011XX",
                            status=status, elapsed_seconds=1)
                   for i, status in enumerate(("error", "reused", "skipped_input")))
    stats = _timing_statistics(entries)
    assert stats["n"] == 4
    assert stats["mean"] == stats["median"] == 150
    assert stats["p10"] == pytest.approx(78)
    assert stats["p90"] == pytest.approx(222)
    assert stats["p95"] == pytest.approx(231)
    assert (stats["min"], stats["fastest_expnum"]) == (60, 100)
    assert (stats["max"], stats["slowest_expnum"]) == (240, 103)


def test_timing_statistics_before_and_after_first_success():
    assert _timing_statistics([]) == {"n": 0}
    entry = Exposure(None, 100, 60623, 11111, "1.3.2", "0011XX", status="complete", elapsed_seconds=42)
    stats = _timing_statistics([entry])
    assert all(stats[key] == 42 for key in ("mean", "median", "p10", "p90", "p95", "min", "max"))
    for invalid in (0, float("nan"), float("inf")):
        entry.elapsed_seconds = invalid
        assert _timing_statistics([entry]) == {"n": 0}


def test_twenty_thousand_exposures_render_as_one_compact_summary():
    entries = [Exposure(None, i, 60623, 11111, "1.3.2", "0011XX", status="complete", elapsed_seconds=120)
               for i in range(20000)]
    stream = StringIO()
    Console(file=stream, width=120, color_system=None).print(_progress_panel(entries, 2400))
    text = stream.getvalue()
    assert "20000/20000 processed" in text and "n=20000" in text
    assert "mean 00:02:00" in text and "P95 00:02:00" in text
    assert len(text.splitlines()) <= 8


@pytest.mark.parametrize("failed,level", [(False, "INFO"), (False, "DEBUG"),
                                        (True, "INFO"), (True, "DEBUG")])
def test_terminal_events_are_bounded_and_detailed_logs_retained(tmp_path, monkeypatch, capsys, failed, level):
    entries = [Exposure(tmp_path / f"input-{i}.fits", i, 60623, 11111, "1.3.2", "0011XX",
                        output=tmp_path / f"output-{i}.fits") for i in range(10)]
    def launch(command, **kwargs):
        request = json.loads(Path(command[-1]).read_text())
        entry = request["exposure"]
        if not failed:
            Path(entry["output"]).touch()
        Path(request["result_path"]).write_text(json.dumps(dict(
            state="error" if failed else "complete", error="test failure", elapsed_seconds=10)))
        return SimpleNamespace(pid=100 + entry["expnum"], poll=lambda: int(failed))
    monkeypatch.setattr(batch.subprocess, "Popen", launch)
    monkeypatch.setattr(batch, "apply_affinity", lambda _: None)
    result = batch.run_batch(entries, {}, dict(processes=1, jobs=1, cpus=1, cpu_ids=[0], subsets=[[0]]),
                             output_dir=tmp_path, log_level=level, no_progress=True)
    terminal = capsys.readouterr().err
    details = (Path(result["run_dir"]) / "batch.log").read_text()
    if failed:
        assert terminal.count("failed: test failure") == (10 if level == "DEBUG" else 3)
        assert details.count("test failure") == 10
        assert result["timing_seconds"] == {"n": 0}
    else:
        assert terminal.count("completed in") == (10 if level == "DEBUG" else 0)
        assert details.count("complete in") == 10
        assert result["timing_seconds"]["n"] == 10
