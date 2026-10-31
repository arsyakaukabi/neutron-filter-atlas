"""Tiny step runner shared by the pipelines/ entry points.

A pipeline is an ordered list of Step(name, function, env). `env` names where a step can run:
  "main"   the main environment (environment.yml)            -> run by default
  "openmc" OpenMC + ENDF data (OPENMC_CROSS_SECTIONS set)     -> run only with --with openmc
  "mlip"   SevenNet / MACE / ORB environments                 -> run only with --with mlip
  "network" needs internet / API keys (Materials Project)     -> run only with --with network
Command line: --list (show steps), --only a,b (run these), --from a (start at a), --with openmc,mlip,network.
"""
import argparse
import time
from dataclasses import dataclass, field
from typing import Callable


@dataclass
class Step:
    name: str
    fn: Callable
    env: str = "main"
    note: str = ""
    kwargs: dict = field(default_factory=dict)


def run(steps, description):
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--from", dest="start", default="")
    ap.add_argument("--with", dest="extra", default="")
    a = ap.parse_args()
    if a.list:
        for s in steps:
            print(f"{s.name:28s} [{s.env}] {s.note}")
        return
    only = set(filter(None, a.only.split(",")))
    extra = set(filter(None, a.extra.split(","))) | {"main"}
    started = not a.start
    for s in steps:
        started = started or s.name == a.start
        if not started or (only and s.name not in only):
            continue
        if s.env not in extra and not (only and s.name in only):
            print(f"-- skip {s.name} (needs environment '{s.env}'; use --with {s.env})", flush=True)
            continue
        t0 = time.time()
        print(f"== {s.name}", flush=True)
        s.fn(**s.kwargs)
        print(f"== {s.name} done in {time.time() - t0:.0f} s", flush=True)
