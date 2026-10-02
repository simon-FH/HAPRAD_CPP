#!/usr/bin/env python3
"""Interpolation error budget of an RC grid (PLAN.md Phase 3).

Samples selected events, calls HAPRAD directly at each one's own kinematics
(MakeRCGrid check, in parallel chunks) and compares with the grid's
interpolated RC factor. Reports the relative error of

    rc_noex = r_vr + r_in              (what ApplyRC uses by default)
    rc_ex   = r_vr + r_in + naz r_ex   (with the exclusive tail)

overall, for events with full and partial grid support, and per bin of each
grid variable -- which says which axis needs more nodes.

    ./check_rc_grid.py <config> <grid.root> <table.root> <n_events> <ntuple files ...>
        [--jobs N] [--naz 0.5] [--save points.txt]

The grid is fine enough when these errors sit well below the statistical
errors of the analysis bins.
"""

import argparse
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
MAKERCGRID = os.path.join(HERE, "bin", "MakeRCGrid")


def quantile(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * (len(v) - 1) + 0.5))] if v else float("nan")


def describe(rel):
    a = [abs(x) for x in rel]
    if not a:
        return "     --"
    mean = sum(rel) / len(rel)
    return "%6d   %+8.4f%%  %8.4f%%  %8.4f%%  %8.4f%%  %8.4f%%" % (
        len(a), 100 * mean, 100 * quantile(a, 0.5), 100 * quantile(a, 0.68), 100 * quantile(a, 0.95), 100 * max(a))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("grid")
    ap.add_argument("table")
    ap.add_argument("n", type=int)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--naz", type=float, default=0.5)
    ap.add_argument("--save", help="write the per-event comparison here")
    args = ap.parse_args()

    data = os.path.realpath(os.environ.get("HAPRAD_DATA_DIR", os.path.join(HERE, "..", "haprad2")))
    if not os.path.isfile(os.path.join(data, "pi_n_maid.dat")):
        sys.exit("no pi_n_maid.dat in %s (set HAPRAD_DATA_DIR)" % data)
    absp = lambda p: os.path.realpath(p)
    files = [absp(f) for f in args.files]

    def chunk(i):
        cmd = [MAKERCGRID, "check", absp(args.config), absp(args.grid), absp(args.table), str(args.n), str(i),
               str(args.jobs)] + files
        r = subprocess.run(cmd, cwd=data, capture_output=True, text=True)
        if r.returncode:
            sys.exit("chunk %d failed:\n%s" % (i, r.stderr[-2000:]))
        return [l for l in r.stdout.splitlines() if l.startswith("@ ")]

    with ThreadPoolExecutor(args.jobs) as pool:
        lines = [l for part in pool.map(chunk, range(args.jobs)) for l in part]

    pts, nodirect = [], 0
    for l in lines:
        f = l.split()
        if "no direct result" in l:
            nodirect += 1
            continue
        v = [float(x) for x in f[2:7]]
        status, cover = int(f[7]), float(f[8])
        gi = [float(x) for x in f[9:12]]
        gd = [float(x) for x in f[12:15]]
        fi, fd = gi[0] + gi[1], gd[0] + gd[1]
        xi, xd = fi + args.naz * gi[2], fd + args.naz * gd[2]
        pts.append(dict(v=v, status=status, cover=cover, noex=(fi - fd) / fd, ex=(xi - xd) / xd, fd=fd, fi=fi))
    if args.save:
        with open(args.save, "w") as out:
            out.write("# Q2 nu z pt phi status cover rc_noex_direct rc_noex_interp rel_err_noex rel_err_ex\n")
            for p in pts:
                out.write("%s %d %.4f %.6f %.6f %+.6f %+.6f\n" % (" ".join("%.5g" % x for x in p["v"]), p["status"],
                                                                 p["cover"], p["fd"], p["fi"], p["noex"], p["ex"]))

    print("events compared: %d  (%d had no direct HAPRAD result)" % (len(pts), nodirect))
    print("\nrelative error, interpolated vs direct")
    print("%-28s %6s   %9s  %9s  %9s  %9s  %9s" % ("", "n", "mean", "median|.|", "68%|.|", "95%|.|", "max|.|"))
    for label, sel in (("rc_noex, all", lambda p: True), ("rc_noex, full support", lambda p: p["status"] == 0),
                       ("rc_noex, partial support", lambda p: p["status"] == 1)):
        print("%-28s %s" % (label, describe([p["noex"] for p in pts if sel(p)])))
    print("%-28s %s" % ("rc_ex (naz=%g), all" % args.naz, describe([p["ex"] for p in pts])))

    names = ["Q2", "nu", "z", "pt", "phi"]
    print("\nrc_noex, full support: mean / 68%|.| of the relative error per bin of each variable")
    full = [p for p in pts if p["status"] == 0]
    for d, name in enumerate(names):
        vals = [p["v"][d] for p in full]
        if not vals:
            continue
        lo, hi = min(vals), max(vals)
        nb = 6
        row = []
        for b in range(nb):
            a, c = lo + (hi - lo) * b / nb, lo + (hi - lo) * (b + 1) / nb
            sel = [p["noex"] for p in full if (a <= p["v"][d] < c) or (b == nb - 1 and p["v"][d] == hi)]
            row.append("[%5.3g,%5.3g) %+6.2f/%5.2f%%" % (a, c, 100 * sum(sel) / len(sel), 100 * quantile([abs(x) for x in sel], .68))
                       if sel else "[%5.3g,%5.3g)      --" % (a, c))
        print("  %-4s " % name + "  ".join(row))

    print("\nlargest errors (rc_noex):")
    for p in sorted(pts, key=lambda p: -abs(p["noex"]))[:8]:
        print("  Q2 %5.2f nu %5.2f z %4.2f pt %5.3f phi %5.1f  status %d cover %.2f  direct %.4f interp %.4f  %+.2f%%" % (
            *p["v"], p["status"], p["cover"], p["fd"], p["fi"], 100 * p["noex"]))


if __name__ == "__main__":
    main()
