"""Single-machine scheduler for the SevenNet atlas run (infrastructure only, no physics): one GPU and one CPU slot.
Machine settings in configs/compute.yaml (section scheduler); jobs run `python -m filteratlas.dynamics.mlip_atlas_run`.
Written for Windows; works on POSIX too (detached sessions).

Two slots over stage3_mlip_set.json, one mlip_atlas_run process per compound:
- GPU slot (python_gpu, CUDA build of torch): smallest supercells first (largest-first spent most of its time on
  out-of-memory tries). A CUDA out-of-memory sends the compound to the CPU overflow queue. When adaptive_failures of
  the last adaptive_window GPU tries ran out of memory, the GPU cap drops below the smallest of those failures.
- CPU slot (python_cpu, cpu_threads OMP threads): smallest supercells first; when the regular queue is empty it takes
  the overflow queue (largest that fits in RAM).
RAM admission for every start: available RAM - expected peak of the new job - not-yet-reached peaks of running jobs
>= RESERVE_GB. Expected host peak: CPU max(min_gb, gb_per_atom n + intercept_gb) (fit to the MLIP-benchmark runs +
max residual); GPU job gpu_host_gb host RAM.
A compound that fits neither slot is logged as not_run_memory and stays flagged.
Jobs are started detached (own process group), so stopping the scheduler does not stop running jobs; a restarted
scheduler skips compounds that have an output or a final CPU error, and adopts nothing (it waits for running jobs).
STOP: create the stop file (path in the constant STOP; no new starts, running jobs finish). Stopping jobs: by PID only.
-> D("logs")/laptop_scheduler_stage3.jsonl
Usage (MLIP environment): python pipelines/03_dynamics_mlip.py --with mlip, or python -m filteratlas.dynamics.scheduler
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import psutil

from filteratlas.common.config import REPO, cfg
from filteratlas.common.paths import D, F

# Machine settings: configs/compute.yaml, section "scheduler" (values used for the atlas run).
_SC = cfg("compute")["scheduler"]
OUT = D("mlip_atlas") / cfg("mlip")["atlas_model"]
LOG = D("logs") / "mlip_atlas.jsonl"
SLOG = D("logs") / "laptop_scheduler_stage3.jsonl"
PY_CPU = Path(os.environ.get("FA_PY_CPU", REPO / _SC["python_cpu"]))
PY_GPU = Path(os.environ.get("FA_PY_GPU", REPO / _SC["python_gpu"]))
STOP = D("logs") / "STOP_STAGE3_LAPTOP"
RESERVE_GB, GPU_HOST_GB, CPU_THREADS = _SC["ram_reserve_gb"], _SC["gpu_host_gb"], str(_SC["cpu_threads"])
GPU_MAX_ATOMS = _SC["gpu_max_atoms"]
# second pass: the adaptive cap can skip compounds <= GPU_MAX_ATOMS untried; FA_NO_ADAPTIVE=1 tries each of them
# once, Tier A first
NO_ADAPTIVE = os.environ.get("FA_NO_ADAPTIVE") == "1"
if os.name == "nt":   # detached, below-normal priority: stopping the scheduler does not stop running jobs
    FLAGS = subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS | subprocess.BELOW_NORMAL_PRIORITY_CLASS
    POPEN_KW = {"creationflags": FLAGS}
else:
    POPEN_KW = {"start_new_session": True}
JOB_MODULE = "filteratlas.dynamics.mlip_atlas_run"


def need_cpu(n):
    m = _SC["cpu_memory_model"]   # expected peak RSS of a CPU job [GB] for n supercell atoms (fit to the MLIP-benchmark runs)
    return max(m["min_gb"], m["gb_per_atom"] * n + m["intercept_gb"])


def last_log():
    out = {}
    if LOG.exists():
        for line in open(LOG, encoding="utf-8"):
            if line.strip():
                r = json.loads(line)
                out[r["material_id"]] = r
    return out


def gpu_oom(rec):
    e = (rec or {}).get("error", "")
    return rec is not None and rec.get("device") == "cuda" and ("CUDA" in e or "OutOfMemory" in e or "cuda" in e)


def gpu_mapping_error(rec):
    """'device not ready / Failed to create GPU mapping' came from expandable_segments under the Windows GPU driver:
    retry on the GPU."""
    e = (rec or {}).get("error", "")
    return gpu_oom(rec) and ("device not ready" in e or "GPU mapping" in e)


def adopt():
    """{material_id: (psutil.Process, device)} for mlip_atlas_run jobs already running (started by an earlier
    scheduler); the venv python.exe is a launcher, so the interpreter child is the one that is tracked."""
    out = {}
    for pr in psutil.process_iter(["cmdline"]):
        cl = pr.info.get("cmdline") or []
        if JOB_MODULE in cl and "--ids" in cl:
            mid = cl[cl.index("--ids") + 1]
            dev = cl[cl.index("--device") + 1] if "--device" in cl else "cpu"
            if mid not in out or pr.ppid() == out[mid][0].pid:
                out[mid] = (pr, dev)
    return out


def avail_gb():
    return psutil.virtual_memory().available / 2 ** 30


def rss_gb(pid):
    """RSS of the job: the launcher process plus its children (the real interpreter)."""
    try:
        p = psutil.Process(pid)
        return sum(q.memory_info().rss for q in [p] + p.children(recursive=True)) / 2 ** 30
    except psutil.Error:
        return 0.0


def start(mid, device):
    env = {**os.environ, "OMP_NUM_THREADS": "2" if device == "cuda" else CPU_THREADS, "PYTHONIOENCODING": "utf-8"}
    if device != "cuda":   # expandable_segments is not used: under the Windows GPU driver it failed with "device not ready"
        env["CUDA_VISIBLE_DEVICES"] = ""
    py = PY_GPU if device == "cuda" else PY_CPU
    env["PYTHONPATH"] = str(REPO / "src")
    err = open(D("logs") / "laptop_scheduler_stage3.err", "a")
    return subprocess.Popen([str(py), "-W", "ignore", "-m", JOB_MODULE, "--threads",
                             "2" if device == "cuda" else CPU_THREADS, "--device", device, "--ids", mid],
                            cwd=REPO, env=env, stdout=subprocess.DEVNULL, stderr=err, **POPEN_KW)


def main():
    comp = json.loads(F("stage3_mlip_set.json").read_text())["compounds"]
    size = {c["material_id"]: (c["supercell_atoms"], c["formula"]) for c in comp}
    logged = last_log()
    busy = adopt()
    todo, overflow = [], []
    for c in comp:
        mid = c["material_id"]
        if (OUT / f"{mid}.json.gz").exists() or mid in busy:
            continue
        if gpu_mapping_error(logged.get(mid)) and not any(
                json.loads(x).get("material_id") == mid and json.loads(x).get("retry_after_mapping_error")
                for x in open(SLOG, encoding="utf-8") if x.strip()):
            todo.append(mid)                                # one GPU retry without expandable_segments
            continue
        r = logged.get(mid)
        if r is not None and r.get("status") == "error" and not gpu_oom(r):
            continue                                        # final CPU failure: recorded, not retried
        (overflow if gpu_oom(r) else todo).append(mid)
    # supercells above GPU_MAX_ATOMS go to the CPU overflow queue; those that cannot fit the RAM end not_run_memory
    overflow += [m for m in todo if size[m][0] > GPU_MAX_ATOMS]
    todo = [m for m in todo if size[m][0] <= GPU_MAX_ATOMS]
    tier_a = {c["material_id"] for c in comp if c["arm"].startswith("atlas_a")}
    rank = (lambda m: (m not in tier_a, size[m][0], m)) if NO_ADAPTIVE else (lambda m: (size[m][0], m))
    todo.sort(key=rank)                                     # ascending (second pass: Tier A first); both lanes take the front
    gpu_cap, gpu_hist = GPU_MAX_ATOMS, []                   # gpu_hist: (atoms, ran out of memory) per finished GPU try
    slog = open(SLOG, "a")
    running = {}
    for mid, (pr, dev) in busy.items():
        running[pr.pid] = {"proc": None, "ps": pr, "mid": mid, "n": size.get(mid, (0, ""))[0], "device": dev,
                           "t0": time.time(), "peak": 0.0}
    print(f"{len(todo)} to run, {len(overflow)} in CPU overflow, adopted {sorted(busy)}", flush=True)
    while todo or overflow or running:
        for pid, j in list(running.items()):
            j["peak"] = max(j["peak"], rss_gb(pid))
            alive = (j["proc"].poll() is None) if j["proc"] is not None else j["ps"].is_running()
            if alive:
                continue
            mid = j["mid"]
            rec = {"material_id": mid, "formula": size[mid][1], "supercell_atoms": j["n"], "device": j["device"],
                   "rc": j["proc"].returncode if j["proc"] is not None else None, "peak_rss_gb": round(j["peak"], 2), "wall_s": round(time.time() - j["t0"], 1),
                   "output": (OUT / f"{mid}.json.gz").exists(), "t_end": time.strftime("%Y-%m-%d %H:%M:%S")}
            if j["device"] == "cuda" and not rec["output"]:
                lr = last_log().get(mid)
                if j.get("retry"):
                    rec["retry_after_mapping_error"] = True
                if gpu_mapping_error(lr) and not j.get("retry"):
                    todo.append(mid)                        # one GPU retry without expandable_segments
                    todo.sort(key=rank)
                    rec["requeued_gpu_retry"] = True
                elif gpu_oom(lr) or lr is None:
                    overflow.append(mid)
                    rec["to_cpu_overflow"] = True
            if j["device"] == "cuda":
                gpu_hist.append((j["n"], not rec["output"] and (gpu_oom(last_log().get(mid)) or last_log().get(mid) is None)))
                last = gpu_hist[-_SC["adaptive_window"]:]
                if not NO_ADAPTIVE and len(last) == _SC["adaptive_window"] and sum(f for _, f in last) >= _SC["adaptive_failures"]:
                    gpu_cap = min(n for n, f in last if f) - 1   # adaptive cap after repeated out-of-memory
                    overflow += [m for m in todo if size[m][0] > gpu_cap]
                    todo = [m for m in todo if size[m][0] <= gpu_cap]
                    gpu_hist.clear()
                    rec["gpu_cap_now"] = gpu_cap
            slog.write(json.dumps(rec) + "\n")
            slog.flush()
            print(json.dumps(rec), flush=True)
            del running[pid]
        if STOP.exists():
            if not running:
                break
            time.sleep(15)
            continue
        pending = sum(max(0.0, (GPU_HOST_GB if j["device"] == "cuda" else need_cpu(j["n"])) - rss_gb(p))
                      for p, j in running.items())
        free = avail_gb() - pending - RESERVE_GB
        lanes = {j["device"] for j in running.values()}
        if "cuda" not in lanes and todo and free > GPU_HOST_GB:
            mid = todo.pop(0)                               # smallest remaining
            p = start(mid, "cuda")
            running[p.pid] = {"proc": p, "mid": mid, "n": size[mid][0], "device": "cuda", "t0": time.time(), "peak": 0.0,
                              "retry": gpu_mapping_error(logged.get(mid))}
            print(f"start GPU {size[mid][1]} ({size[mid][0]} atoms) avail {avail_gb():.1f} GB", flush=True)
            continue
        if "cpu" not in lanes:
            pick = None
            if todo and need_cpu(size[todo[0]][0]) < free:
                pick = todo.pop(0)                          # smallest remaining
            elif not todo and overflow:
                fit = [m for m in sorted(overflow, key=lambda m: -size[m][0]) if need_cpu(size[m][0]) < free]
                if fit:
                    pick = fit[0]
                    overflow.remove(pick)
                elif not running:                           # nothing fits even with the machine idle
                    for m in overflow:
                        rec = {"material_id": m, "formula": size[m][1], "supercell_atoms": size[m][0],
                               "status": "not_run_memory", "need_gb": round(need_cpu(size[m][0]), 1),
                               "avail_gb": round(avail_gb(), 1)}
                        slog.write(json.dumps(rec) + "\n")
                        print(json.dumps(rec), flush=True)
                    overflow = []
            if pick:
                p = start(pick, "cpu")
                running[p.pid] = {"proc": p, "mid": pick, "n": size[pick][0], "device": "cpu", "t0": time.time(), "peak": 0.0}
                print(f"start CPU {size[pick][1]} ({size[pick][0]} atoms) avail {avail_gb():.1f} GB", flush=True)
        time.sleep(2)                                       # short poll: the GPU idles between jobs otherwise
    print("SCHEDULER_DONE", flush=True)


if __name__ == "__main__":
    main()
