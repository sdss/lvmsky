"""Batch display counters, honest ETA and terminal/file logging separation."""
from io import StringIO

from rich.console import Console

from skysub.sky_subtraction_block.batch import _progress_panel
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
