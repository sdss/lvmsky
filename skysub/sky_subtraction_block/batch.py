"""Bounded subprocess scheduling and durable batch run records."""
from __future__ import annotations

from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path
import json
import logging
import os
import shlex
import signal
import subprocess
import sys
import tempfile
import time

from .resources import apply_affinity


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(value, stream, indent=2, default=str)
            stream.write("\n")
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def new_run_dir(output_dir):
    return Path(output_dir) / ".batch" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")


def write_table(path, entries):
    from astropy.table import Table

    rows = []
    for entry in entries:
        rows.append(dict(expnum=entry.expnum, mjd=entry.mjd, tileid=entry.tileid, version=entry.version,
                         input=str(entry.input_path), output=str(entry.output), seed=str(entry.seed or ""),
                         state_mode=entry.state_mode, status=entry.status, error=entry.error,
                         fingerprint=entry.fingerprint, pid=entry.pid or -1,
                         elapsed_seconds=entry.elapsed_seconds, log=entry.log))
    # Only the coordinator writes these records.
    temporary = path.with_name("." + path.name)
    Table(rows=rows).write(temporary, format="ascii.ecsv", overwrite=True)
    os.replace(temporary, path)


def _terminate(active):
    for process, stream, result_path, slot, entry in active.values():
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
    for process, stream, result_path, slot, entry in active.values():
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
        stream.close()
        entry.status = "interrupted"


def _duration(seconds):
    minutes, seconds = divmod(max(0, int(seconds)), 60)
    hours, minutes = divmod(minutes, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def _progress_panel(entries, elapsed):
    from rich.console import Group
    from rich.progress_bar import ProgressBar
    from rich.table import Table
    from rich.text import Text

    counts = Counter(entry.status for entry in entries)
    completed = sum(counts[name] for name in ("complete", "error", "reused", "skipped_input"))
    remaining = counts["running"] + counts["queued"]
    rate = counts["complete"] / elapsed if elapsed > 0 else 0
    eta = "00:00:00" if not remaining else ("~" + _duration(remaining / rate) if rate else "--")
    bar = Table.grid(expand=True, padding=(0, 1))
    bar.add_column(no_wrap=True)
    bar.add_column(ratio=1)
    bar.add_column(no_wrap=True)
    bar.add_row("Sky subtraction", ProgressBar(total=len(entries), completed=completed),
                f"{completed}/{len(entries)} processed")
    summary = (f"OK {counts['complete']} · Failed {counts['error']} · "
               f"Running {counts['running']} · Queued {counts['queued']}")
    for status, label in (("reused", "Reused"), ("skipped_input", "Skipped"),
                          ("not_started", "Not started"), ("interrupted", "Interrupted")):
        if counts[status]:
            summary += f" · {label} {counts[status]}"
    running = ", ".join(str(entry.expnum) for entry in entries if entry.status == "running") or "none"
    return Group(bar, Text(summary),
                 Text(f"Elapsed {_duration(elapsed)} · Rate {rate * 60:.1f} successful exp/min · ETA {eta}"),
                 Text(f"Running: {running}"))


def run_batch(entries, model, budget, *, output_dir, run_dir=None, context=None,
              on_error="continue", log_level="INFO", no_progress=False):
    run_dir = Path(run_dir).expanduser().resolve() if run_dir else new_run_dir(output_dir).resolve()
    if run_dir.exists() and any(run_dir.iterdir()):
        raise ValueError(f"Run directory must be new/empty: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    state = dict(state="running", started=timestamp(), command=shlex.join([sys.executable, "-m", "skysub.sky_subtraction_block.cli", *sys.argv[1:]]),
                 environment={name: os.environ.get(name) for name in ("SAS_BASE_DIR", "LVM_SPECTRO_REDUX")},
                 budget=budget, paths=context or {}, run_dir=str(run_dir))
    queue = deque(entry for entry in entries if entry.status == "queued")
    active = {}
    slots = deque(range(budget["processes"]))
    apply_affinity(budget["cpu_ids"])
    write_table(run_dir / "plan.ecsv", entries)

    def save():
        state.update(updated=timestamp(), elapsed_seconds=time.monotonic() - started,
                     counts=dict(Counter(e.status for e in entries)),
                     active=[dict(expnum=e.expnum, pid=p.pid, cpu_ids=budget["subsets"][slot])
                             for p, stream, result, slot, e in active.values()])
        atomic_json(run_dir / "status.json", state)
        write_table(run_dir / "results.ecsv", entries)

    from rich.console import Console
    from rich.text import Text

    console = Console(stderr=True)
    log = logging.Logger(__name__, level=logging.DEBUG)
    handler = logging.FileHandler(run_dir / "batch.log")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(handler)
    log.info("Command: %s; paths=%s; CPU budget=%s", state["command"], context, budget)
    terminal_level = logging.getLevelName(log_level.upper())
    def event(message, level=logging.INFO, style=None):
        if level >= terminal_level:
            console.print(Text(message, style=style))

    live = None
    if not no_progress and sys.stderr.isatty():
        from rich.live import Live
        live = Live(console=console, refresh_per_second=4,
                    get_renderable=lambda: _progress_panel(entries, time.monotonic() - started))
    try:
        save()
        console.print(Text(f"Run records: {run_dir}"))
        if live:
            live.start()
        stop = False
        while queue or active:
            while queue and slots and not stop:
                slot = slots.popleft()
                entry = queue.popleft()
                stem = f"{entry.expnum:08d}-{entries.index(entry):05d}"
                request_path = run_dir / (stem + ".request.json")
                result_path = run_dir / (stem + ".result.json")
                entry.log = str(run_dir / (stem + ".log"))
                request = dict(exposure=entry.record(), model=model,
                               cpu_ids=budget["subsets"][slot], budget=budget,
                               result_path=str(result_path), log_level=log_level,
                               command=state["command"])
                atomic_json(request_path, request)
                stream = Path(entry.log).open("w")
                command = [sys.executable, "-m", "skysub.sky_subtraction_block._worker", str(request_path)]
                environment = os.environ.copy()
                repository = str(Path(__file__).resolve().parents[2])
                environment["PYTHONPATH"] = os.pathsep.join(filter(None, (repository, environment.get("PYTHONPATH"))))
                try:
                    process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT,
                                               start_new_session=True, env=environment)
                except BaseException:
                    stream.close()
                    raise
                entry.pid, entry.status = process.pid, "running"
                active[process.pid] = (process, stream, result_path, slot, entry)
                log.info("Exposure %s started pid=%s CPUs=%s: %s -> %s; log=%s",
                         entry.expnum, process.pid, budget["subsets"][slot],
                         entry.input_path, entry.output, entry.log)
                event(f"Exposure {entry.expnum} started", level=logging.DEBUG)
                save()
            for pid, record in list(active.items()):
                process, stream, result_path, slot, entry = record
                code = process.poll()
                if code is None:
                    continue
                stream.close()
                try:
                    result = json.loads(result_path.read_text())
                except (OSError, ValueError):
                    result = dict(error=f"Worker exited {code}; inspect {entry.log}")
                entry.elapsed_seconds = float(result.get("elapsed_seconds", 0))
                if code == 0 and result.get("state") == "complete" and entry.output.is_file():
                    entry.status = "complete"
                    log.info("Exposure %s complete in %.1f s: %s", entry.expnum, entry.elapsed_seconds, entry.output)
                    event(f"✓ {entry.expnum} completed in {_duration(entry.elapsed_seconds)}", style="green")
                else:
                    entry.status, entry.error = "error", result.get("error", f"Worker exited {code}")
                    log.error("Exposure %s: %s; log=%s", entry.expnum, entry.error, entry.log)
                    reason = " ".join(str(entry.error).split())
                    if len(reason) > 160:
                        reason = reason[:157] + "..."
                    event(f"✗ {entry.expnum} failed: {reason} — log: {Path(entry.log).name}",
                          level=logging.ERROR, style="red")
                    stop = on_error == "stop"
                del active[pid]
                slots.append(slot)
                save()
            if stop:
                for entry in queue:
                    entry.status = "not_started"
                    entry.error = "Submission stopped by --on-error stop"
                queue.clear()
            if active:
                time.sleep(0.1)
        state["state"] = "complete_with_errors" if any(e.status == "error" for e in entries) else "complete"
    except KeyboardInterrupt:
        log.warning("Batch interrupted by user")
        state["state"] = "interrupted"
        _terminate(active)
        active.clear()
        for entry in queue:
            entry.status = "not_started"
            entry.error = "Interrupted before submission"
        state["state"] = "interrupted"
        raise
    except BaseException:
        log.exception("Batch failed")
        _terminate(active)
        active.clear()
        state["state"] = "error"
        raise
    finally:
        try:
            if live:
                live.stop()
            save()
        finally:
            log.removeHandler(handler)
            handler.close()
    return state
