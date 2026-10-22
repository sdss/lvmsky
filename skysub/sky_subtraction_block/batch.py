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

LOG = logging.getLogger(__name__)


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

    progress = None
    if not no_progress and sys.stderr.isatty():
        from rich.progress import Progress
        progress = Progress()
        progress.start()
        bar = progress.add_task("Sky subtraction", total=len(entries),
                                completed=sum(e.status != "queued" for e in entries))
    save()
    print(f"Run records: {run_dir}", flush=True)
    try:
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
                LOG.info("Exposure %s started pid=%s CPUs=%s: %s -> %s; log=%s",
                         entry.expnum, process.pid, budget["subsets"][slot],
                         entry.input_path, entry.output, entry.log)
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
                    LOG.info("Exposure %s complete: %s", entry.expnum, entry.output)
                else:
                    entry.status, entry.error = "error", result.get("error", f"Worker exited {code}")
                    LOG.error("Exposure %s: %s; log=%s", entry.expnum, entry.error, entry.log)
                    stop = on_error == "stop"
                del active[pid]
                slots.append(slot)
                if progress:
                    progress.advance(bar)
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
        state["state"] = "interrupted"
        _terminate(active)
        active.clear()
        for entry in queue:
            entry.status = "not_started"
            entry.error = "Interrupted before submission"
        state["state"] = "interrupted"
        raise
    except BaseException:
        _terminate(active)
        active.clear()
        state["state"] = "error"
        raise
    finally:
        if progress:
            progress.stop()
        save()
    return state
