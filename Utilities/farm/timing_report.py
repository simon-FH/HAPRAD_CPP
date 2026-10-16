#!/usr/bin/env python3
"""Benchmark of an RC workflow from the timing records its jobs wrote.

    ./timing_report.py <workflow dir> [--runs N] [--chunk-minutes M]

Reads <workflow dir>/jobs.tsv (every job and what it waits for, written by
submit_rc_workflow.sh) and timing/*.tsv (one line per finished job, written by
rc_job.sh), and prints:

  1. the timeline: submission, first start, last end, and each iteration's span;
  2. per job type: jobs done / failed / missing, wall time, CPU, memory, and the
     wait between a job becoming ready (its dependencies done) and starting --
     the farm's queue latency;
  3. throughput: CPU per million events (fill, apply) and per grid node;
  4. with --runs N: the CPU and the wall time at full parallelism expected for
     N runs. fill and apply scale with the number of runs; fit, grid and merge
     do not (the grid's cost is set by its nodes, not by the data), so those are
     valid only if this workflow used the production config;
  5. suggested swif2 resources (3x the longest wall time, 2x the largest
     memory) and the number of grid chunks for chunks of --chunk-minutes.
"""

import argparse
import glob
import math
import os
import statistics
import sys
import time

STEPS = ["fill", "fit", "grid", "merge", "apply"]
PER_RUN = {"fill", "apply"}


def read_jobs(d):
    """jobs.tsv as {name: dict(step, iteration, arg, after=[names])}, in order."""
    jobs = {}
    with open(os.path.join(d, "jobs.tsv")) as f:
        next(f)
        for line in f:
            name, step, it, arg, after = (line.rstrip("\n").split("\t") + [""])[:5]
            jobs[name] = dict(step=step, iteration=int(it), arg=arg, after=[a for a in after.split(",") if a])
    return jobs


def num(x):
    """A float from a timing field, None for NA or empty."""
    try:
        return float(x)
    except ValueError:
        return None


def read_timing(d):
    """timing/*.tsv as {name: dict(...)}."""
    out = {}
    for path in glob.glob(os.path.join(d, "timing", "*.tsv")):
        with open(path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 14:
                    continue
                out[p[0]] = dict(step=p[1], iteration=int(p[2]), arg=p[3], host=p[4], start=float(p[5]),
                                 end=float(p[6]), wall=num(p[7]), user=num(p[8]), sys=num(p[9]),
                                 rss=num(p[10]), status=int(p[11]), w1=num(p[12]), w2=num(p[13]))
    return out


def settings(d):
    """KEY=value lines of the workflow's farm.cfg copy (last one wins)."""
    s = {}
    with open(os.path.join(d, "farm.cfg")) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if "=" in line:
                k, v = line.split("=", 1)
                s[k.strip()] = v.strip().strip('"')
    return s


def hms(sec):
    """Seconds as e.g. 1h05m, 12m30s or 45s."""
    if sec is None:
        return "--"
    sec = int(round(sec))
    if sec >= 3600:
        return "%dh%02dm" % (sec // 3600, sec % 3600 // 60)
    if sec >= 60:
        return "%dm%02ds" % (sec // 60, sec % 60)
    return "%ds" % sec


def cpu(t):
    """User + system CPU of a timing record, or its wall time if CPU is unknown."""
    if t["user"] is not None and t["sys"] is not None:
        return t["user"] + t["sys"]
    return t["wall"]


def swif_time(sec):
    """A swif2 -time value covering sec: whole hours, or 30-minute steps below 2h."""
    if sec <= 1800:
        return "30m"
    if sec <= 3600:
        return "1h"
    return "%dh" % math.ceil(sec / 3600.)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workdir")
    ap.add_argument("--runs", type=int, help="extrapolate to this many runs")
    ap.add_argument("--chunk-minutes", type=float, default=30., help="target wall time of one grid chunk")
    args = ap.parse_args()
    d = os.path.realpath(args.workdir)

    jobs, tim, cfg = read_jobs(d), read_timing(d), settings(d)
    with open(os.path.join(d, "timing", "submitted")) as f:
        submitted = float(f.read())
    with open(os.path.join(d, "runs.txt")) as f:
        runs = [l.strip() for l in f if l.strip()]
    iterations = int(cfg.get("ITERATIONS", 1))
    nchunks = int(cfg.get("NCHUNKS", 1))

    # ---- When each job became ready: its last dependency finished (or submission) ----
    for name, t in tim.items():
        deps = jobs.get(name, {}).get("after", [])
        ends = [tim[a]["end"] for a in deps if a in tim]
        t["ready"] = max(ends) if ends and len(ends) == len(deps) else submitted
        t["wait"] = t["start"] - t["ready"]

    # ==== 1. Timeline ====
    done = [t for t in tim.values()]
    print("workflow %s: %d runs, %d iteration(s), %d grid chunks; %d of %d jobs have a timing record"
          % (os.path.basename(d), len(runs), iterations, nchunks, len(tim), len(jobs)))
    print("  submitted   %s" % time.strftime("%F %T", time.localtime(submitted)))
    if done:
        first, last = min(t["start"] for t in done), max(t["end"] for t in done)
        print("  first start %s  (%s after submission)" % (time.strftime("%F %T", time.localtime(first)), hms(first - submitted)))
        print("  last end    %s  (%s after submission)" % (time.strftime("%F %T", time.localtime(last)), hms(last - submitted)))
    for k in range(iterations):
        ts = [t for t in done if t["iteration"] == k]
        if ts:
            s, e = min(t["start"] for t in ts), max(t["end"] for t in ts)
            print("  iteration %d %s  (first start to last end)" % (k, hms(e - s)))

    # ==== 2. Per job type ====
    print("\n%-6s %5s %5s %5s  %9s %9s  %9s  %8s  %9s %9s" % ("step", "done", "fail", "miss", "wall mean", "wall max",
                                                            "CPU total", "max RSS", "wait mean", "wait max"))
    per = {}
    for step in STEPS:
        names = [n for n, j in jobs.items() if j["step"] == step]
        ts = [tim[n] for n in names if n in tim]
        ok = [t for t in ts if t["status"] == 0]
        per[step] = ok
        if not names:
            continue
        walls = [t["wall"] for t in ok if t["wall"] is not None]
        rss = [t["rss"] for t in ok if t["rss"] is not None]
        waits = [t["wait"] for t in ok]
        print("%-6s %5d %5d %5d  %9s %9s  %9s  %8s  %9s %9s" % (
            step, len(ok), len(ts) - len(ok), len(names) - len(ts),
            hms(statistics.mean(walls)) if walls else "--", hms(max(walls)) if walls else "--",
            hms(sum(cpu(t) for t in ok)) if ok else "--",
            "%.0f MB" % (max(rss) / 1024.) if rss else "--",
            hms(statistics.mean(waits)) if waits else "--", hms(max(waits)) if waits else "--"))
    failed = [n for n, t in tim.items() if t["status"] != 0]
    if failed:
        print("  failed: %s  (see logs/<name>.log)" % ", ".join(sorted(failed)[:8]))

    # ==== 3. Throughput ====
    print("\nthroughput")
    fill, grid, apply_ = per["fill"], per["grid"], per["apply"]
    ev = sum(t["w1"] or 0 for t in fill)
    if ev:
        print("  fill   %8.1f CPU s per million entries read   (%.3g entries per run on average)"
              % (sum(cpu(t) for t in fill) / ev * 1e6, ev / len(fill)))
    nodes = sum(t["w1"] or 0 for t in grid)
    valid = sum(t["w2"] or 0 for t in grid)
    if nodes:
        print("  grid   %8.2f CPU s per node, %.2f per valid node   (%d nodes per iteration, %d valid)"
              % (sum(cpu(t) for t in grid) / nodes, sum(cpu(t) for t in grid) / max(valid, 1),
                 nodes / max(1, len({t["iteration"] for t in grid})), valid / max(1, len({t["iteration"] for t in grid}))))
    ent = sum(t["w1"] or 0 for t in apply_)
    if ent:
        print("  apply  %8.1f CPU s per million entries" % (sum(cpu(t) for t in apply_) / ent * 1e6))

    # ==== 4. Extrapolation ====
    if args.runs:
        n = args.runs
        print("\nexpected for %d runs, %d iteration(s), %d grid chunks (from this workflow's averages)" % (n, iterations, nchunks))
        print("  fill and apply scale with the runs (assumed the size of these); fit, grid and merge are")
        print("  per iteration and valid only for this workflow's config -- benchmark with the production grid")
        cpu_total, path = 0., 0.
        for step in STEPS:
            ok = per[step]
            if not ok:
                print("  %-6s no timing records -- cannot extrapolate" % step)
                continue
            its = max(1, len({t["iteration"] for t in ok}))
            c_iter = sum(cpu(t) for t in ok) / its
            if step in PER_RUN:
                c_iter *= n / float(len(runs))
            cpu_total += c_iter * iterations
            path += max(t["wall"] for t in ok) * iterations
            print("  %-6s %9s CPU per iteration%s" % (step, hms(c_iter), "  (scales with runs)" if step in PER_RUN else ""))
        print("  total  %9s CPU; %s wall at full parallelism (longest job of each step, summed), plus queue waits"
              % (hms(cpu_total), hms(path)))

    # ==== 5. Suggested resources ====
    print("\nsuggested resources (3x the longest wall, 2x the largest RSS, at least 1g)")
    for step in STEPS:
        ok = per[step]
        walls = [t["wall"] for t in ok if t["wall"] is not None]
        rss = [t["rss"] for t in ok if t["rss"] is not None]
        if not walls:
            continue
        ram = max(1, math.ceil(2 * max(rss) / 1024. / 1024.)) if rss else None
        print("  RES_%-6s -time %-4s -ram %s   (now: %s)" % (step.upper(), swif_time(3 * max(walls)),
                                                          "%dg" % ram if ram else "?", cfg.get("RES_" + step.upper(), "?")))
    if grid:
        its = max(1, len({t["iteration"] for t in grid}))
        per_iter = sum(cpu(t) for t in grid) / its
        print("  NCHUNKS = %d for grid chunks of about %g minutes (now %d: %s each)"
              % (max(1, math.ceil(per_iter / (60. * args.chunk_minutes))), args.chunk_minutes, nchunks,
                 hms(statistics.mean(t["wall"] for t in grid))))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
