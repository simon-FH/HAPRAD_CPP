#!/usr/bin/env python3
"""The radiative-correction iteration (PLAN.md Phase 4).

    table_k -> RC grid_k -> weights_k -> table_k+1 from the RC-corrected data

HAPRAD needs Born-level amplitudes, the data give observed ones; iterating
removes the difference. Iteration 0 fills the table from the data as they are
(times `weight_branch`, e.g. acceptance, if the config sets one); iteration k
multiplies in the RC weights of iteration k-1 (MakePhiTable fill --rc-friend).
It stops when the RC factor (without the exclusive tail) changes by less than
--tol at 95% of the grid nodes, or after --max-iter iterations.

    ./rc_iterate.py <config> <workdir> <ntuple files or globs ...>
        [--max-iter 4] [--tol 0.002] [--jobs N] [--resume]

Layout: <workdir>/iter_<k>/{counts,table,grid}.root and friends/<name>_rc.root
(one RC friend per input file -- the last iteration's are the weights for the
analysis). Run per target: each target's own data give its own table.

Cost per iteration is dominated by the RC grid (run_rc_grid.sh, ~30 min on 14
cores for the default grid); on the cluster, replace run_grid() with a job
array of `MakeRCGrid run` chunks and an hadd.
"""

import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")


def run(cmd, log, **kw):
    with open(log, "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True, **kw)
    if r.returncode:
        sys.exit("failed (%d): %s\n  see %s" % (r.returncode, " ".join(cmd), log))


def friend_name(path):
    base = os.path.basename(path)
    return (base[:-5] if base.endswith(".root") else base) + "_rc.root"


def compare(a, b):
    r = subprocess.run([os.path.join(BIN, "MakeRCGrid"), "compare", a, b], capture_output=True, text=True)
    line = [l for l in r.stdout.splitlines() if l.startswith("@compare")]
    if r.returncode or not line:
        sys.exit("grid comparison failed:\n" + r.stderr)
    f = line[0].split()
    return dict(nodes=int(f[1]), max=float(f[2]), median=float(f[3]), p95=float(f[4]), mean=float(f[5]),
                one_only=int(f[6]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("workdir")
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--max-iter", type=int, default=4)
    ap.add_argument("--tol", type=float, default=0.002, help="95%% of nodes must change by less than this")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--resume", action="store_true", help="keep outputs that already exist")
    args = ap.parse_args()

    config = os.path.realpath(args.config)
    work = os.path.realpath(args.workdir)
    files = sorted({os.path.realpath(f) for p in args.inputs for f in (glob.glob(p) or [p])})
    if not files:
        sys.exit("no input files")
    names = [friend_name(f) for f in files]
    if len(set(names)) != len(names):
        sys.exit("two input files share a name; their RC friends would collide")
    os.makedirs(work, exist_ok=True)

    def need(path):
        return not (args.resume and os.path.exists(path))

    history = []
    prev = None
    for k in range(args.max_iter):
        d = os.path.join(work, "iter_%d" % k)
        fr = os.path.join(d, "friends")
        os.makedirs(fr, exist_ok=True)
        counts, table, grid = (os.path.join(d, x) for x in ("counts.root", "table.root", "grid.root"))
        print("iteration %d" % k, flush=True)

        if need(counts):
            extra = ["--rc-friend", os.path.join(prev, "friends")] if prev else []
            run([os.path.join(BIN, "MakePhiTable"), "fill", config, counts] + extra + files, os.path.join(d, "fill.log"))
        if need(table):
            run([os.path.join(BIN, "MakePhiTable"), "fit", config, counts, table], os.path.join(d, "fit.log"))
        if need(grid):
            run([os.path.join(HERE, "run_rc_grid.sh"), config, table, grid, str(args.jobs)], os.path.join(d, "grid.log"))

        def apply(i):
            out = os.path.join(fr, names[i])
            if need(out):
                run([os.path.join(BIN, "ApplyRC"), config, grid, files[i], out], out[:-5] + ".log")

        with ThreadPoolExecutor(args.jobs) as pool:
            list(pool.map(apply, range(len(files))))

        if prev:
            c = compare(os.path.join(prev, "grid.root"), grid)
            history.append((k, c))
            print("  RC change vs iteration %d over %d nodes: 95%% %.2e, median %.2e, max %.2e, mean %+.2e"
                  % (k - 1, c["nodes"], c["p95"], c["median"], c["max"], c["mean"]), flush=True)
            if c["p95"] < args.tol:
                print("converged after %d iterations: the weights are in %s" % (k + 1, fr))
                return 0
        prev = d
    print("NOT converged within %d iterations (tolerance %g); last weights in %s"
          % (args.max_iter, args.tol, os.path.join(prev, "friends")))
    return 1


if __name__ == "__main__":
    sys.exit(main())
