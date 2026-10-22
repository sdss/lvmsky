"""Explicit native-thread and Linux CPU-affinity limits, without numerical imports."""
from __future__ import annotations

import os
import re

THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    "NUMEXPR_MAX_THREADS", "RAYON_NUM_THREADS", "POLARS_MAX_THREADS",
    "NUMBA_NUM_THREADS", "TBB_NUM_THREADS",
)
_LIMITER = None


def limit_environment():
    for name in THREAD_VARIABLES:
        os.environ[name] = "1"
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["MKL_DYNAMIC"] = "FALSE"


def parse_cpu_list(value):
    cpus = []
    for token in value.split(","):
        if not re.fullmatch(r"\d+(?:-\d+)?", token.strip()):
            raise ValueError("--cpu-list accepts CPU IDs/ranges, e.g. 0,2,8-11")
        ends = [int(v) for v in token.strip().split("-")]
        lower, upper = ends[0], ends[-1]
        if upper < lower:
            raise ValueError("--cpu-list ranges must be increasing")
        if upper - lower + 1 > len(os.sched_getaffinity(0)):
            raise ValueError("--cpu-list range exceeds the allowed CPU count")
        cpus.extend(range(lower, upper + 1))
    if len(set(cpus)) != len(cpus):
        raise ValueError("--cpu-list must contain distinct CPU IDs")
    return sorted(cpus)


def cpu_budget(processes=1, jobs=1, cpus=None, cpu_list=None, *, strict=True):
    if processes < 1 or jobs < 1 or (cpus is not None and cpus < 1):
        raise ValueError("--processes, --jobs and --cpus must be positive integers")
    required = processes * jobs
    if cpu_list is not None and cpus is None:
        raise ValueError("Supply --cpus COUNT with --cpu-list IDS")
    if not hasattr(os, "sched_getaffinity"):
        if strict or cpus is not None:
            raise ValueError("Strict CPU control requires Linux CPU affinity")
        return dict(processes=processes, jobs=jobs, cpus=None, cpu_ids=None, subsets=[None])
    allowed = sorted(os.sched_getaffinity(0))
    if not strict and cpus is None:
        return dict(processes=1, jobs=jobs, cpus=len(allowed), cpu_ids=None, subsets=[None])
    count = required if cpus is None else cpus
    if required > count:
        raise ValueError(f"Need processes * jobs = {required} CPUs; --cpus supplies {count}")
    ids = parse_cpu_list(cpu_list) if cpu_list is not None else allowed[:count]
    if len(ids) != count or not set(ids).issubset(allowed):
        raise ValueError(f"Requested {count} CPUs/IDs must fit the allowed set {allowed}")
    return dict(processes=processes, jobs=jobs, cpus=count, cpu_ids=ids,
                subsets=[ids[i * jobs:(i + 1) * jobs] for i in range(processes)])


def apply_affinity(cpu_ids):
    if cpu_ids is not None:
        os.sched_setaffinity(0, set(cpu_ids))


def limit_numerics():
    """Call once before inference; keep the native-pool controller alive."""
    global _LIMITER
    limit_environment()
    import torch
    from threadpoolctl import threadpool_limits
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    _LIMITER = threadpool_limits(limits=1)


def runtime_resources():
    from threadpoolctl import threadpool_info
    import torch
    pools = threadpool_info()
    if (torch.get_num_threads() != 1 or torch.get_num_interop_threads() != 1
            or any(p["num_threads"] != 1 for p in pools)):
        raise RuntimeError("A supported native numerical pool exceeds one thread")
    return dict(pid=os.getpid(), cpu_ids=sorted(os.sched_getaffinity(0)),
                torch_threads=torch.get_num_threads(),
                torch_interop_threads=torch.get_num_interop_threads(), native_pools=pools,
                environment={n: os.environ[n] for n in THREAD_VARIABLES})
