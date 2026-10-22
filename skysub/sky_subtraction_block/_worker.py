"""Fresh-interpreter batch worker; numerical imports follow CPU/thread setup."""
from __future__ import annotations

import json
import logging
from pathlib import Path
import resource
import signal
import sys
import time
import traceback

from .resources import apply_affinity, limit_environment, limit_numerics, runtime_resources

limit_environment()


def main(request_path):
    from .batch import atomic_json
    from .planning import file_sha256

    request = json.loads(Path(request_path).read_text())
    result_path = Path(request["result_path"])
    entry = request["exposure"]
    started = time.monotonic()
    logging.basicConfig(level=request["log_level"], format="%(asctime)s %(levelname)s %(message)s")
    apply_affinity(request["cpu_ids"])

    def interrupted(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupted)
    result = dict(state="running", pid=__import__("os").getpid(), cpu_ids=request["cpu_ids"])
    try:
        for path, expected in entry["dependencies"].items():
            if file_sha256(path) != expected:
                raise ValueError(f"Input dependency changed after planning: {path}")
        limit_numerics()
        from .pipeline import run
        result["resources_start"] = runtime_resources()
        atomic_json(result_path, result)
        metadata = dict(fingerprint=entry["fingerprint"], command=request["command"],
                        state_mode=entry["state_mode"], budget=request["budget"],
                        resources=result["resources_start"])
        output = run(entry["input_path"], output=entry["output"], seed=entry["seed"],
                     execution_metadata=metadata, **request["model"])
        result.update(state="complete", output=str(output), resources_end=runtime_resources())
        from astropy.io import fits
        with fits.open(output, memmap=True) as hdus:
            result["qa"] = dict(skyml_reliability=int(hdus[0].header.get("SKYMLREL", -1)),
                                reference_blocks=int(sum(not bool(v) for v in hdus["BLOCK_LSF_META"].data["adopted"])))
        return 0
    except KeyboardInterrupt:
        result.update(state="interrupted", error="User interruption")
        return 130
    except Exception as error:
        traceback.print_exc()
        result.update(state="error", error=f"{type(error).__name__}: {error}")
        return 1
    finally:
        own = resource.getrusage(resource.RUSAGE_SELF)
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        result.update(elapsed_seconds=time.monotonic() - started,
                      cpu_seconds=own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime,
                      peak_rss_kib=own.ru_maxrss)
        atomic_json(result_path, result)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
