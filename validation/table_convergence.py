#!/usr/bin/env python3
"""How fine must the structure-function table be?  (PLAN.md Phase 2.3)

Builds tables from HAPRAD 2.0's own PDF x FF model at several resolutions,
runs the C++ with each, and compares against the full FORTRAN calculation at
the RG-E test points. With the same physics on both sides, what is left is the
error the table introduces: discretisation and interpolation, propagated
through the Born cross section and the tail integrals.

    ./table_convergence.py                       # default set of grids
    ./table_convergence.py --grid "18:0.5:10,16:1:10.5,20:0:1,12:0:1.5" ...

Needs bin/make_model_table, bin/born_harmonics, bin/rc_point,
bin/haprad2_point (so CERN_LIB must have been set at build time).
"""

import argparse
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")
DATA = os.path.join(HERE, "data")
E = 10.5473

POINTS = [  # (x, Q2, z, pt, phi) at E -- the RG-E grid of compare.py
    (0.234, 2.5, 0.34, 0.424, 60.0),
    (0.234, 2.5, 0.34, 0.424, 180.0),
    (0.150, 2.0, 0.30, 0.350, 90.0),
    (0.300, 4.0, 0.50, 0.500, 30.0),
    (0.234, 2.5, 0.70, 0.300, 60.0),
]

GRIDS = [  # label, "nQ2:lo:hi,nnu:lo:hi,nz:lo:hi,npt2:lo:hi"
    ("base        9x 8x10x 6", "9:0.5:10,8:1:10.5,10:0:1,6:0:1.5"),
    ("pt2 x4      9x 8x10x24", "9:0.5:10,8:1:10.5,10:0:1,24:0:1.5"),
    ("all x2     18x16x20x12", "18:0.5:10,16:1:10.5,20:0:1,12:0:1.5"),
    ("all x2,pt2x4 18x16x20x24", "18:0.5:10,16:1:10.5,20:0:1,24:0:1.5"),
]

KEYS = ("sib", "tai_in", "f1", "f2")


def at(binary, point, env=None):
    cmd = [os.path.join(BIN, binary), "%.10g" % E] + ["%.10g" % v for v in point]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=DATA,
                         env=env, timeout=3600).stdout
    vals, block = {}, None
    for line in out.splitlines():
        if line.startswith("### BEGIN "):
            block = line[10:].strip()
        elif line.startswith("@") and block == "RES":
            k, _, v = line[1:].partition(" ")
            vals[k.strip()] = float(v.replace("D", "e").replace("E", "e"))
    return vals


def build_table(spec, path):
    cells = subprocess.run([os.path.join(BIN, "make_model_table"), "points", spec, str(E)],
                           capture_output=True, text=True, check=True).stdout
    harm = subprocess.run([os.path.join(BIN, "born_harmonics")], input=cells,
                          capture_output=True, text=True, cwd=DATA, check=True).stdout
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(harm)
        hpath = f.name
    r = subprocess.run([os.path.join(BIN, "make_model_table"), "build", spec, str(E), hpath, path],
                       capture_output=True, text=True, cwd=DATA)
    os.unlink(hpath)
    if r.returncode != 0:
        sys.exit("table build failed for %s:\n%s" % (spec, r.stdout + r.stderr))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", action="append", help="extra grid spec(s) to try")
    args = ap.parse_args()
    grids = GRIDS + [("custom %s" % g, g) for g in (args.grid or [])]

    print("FORTRAN reference (full calculation, its own model) ...")
    ref = [at("haprad2_point", p) for p in POINTS]

    tmpdir = tempfile.TemporaryDirectory(prefix="tblconv_")  # removed when main() returns
    tmp = tmpdir.name
    rows = []
    for label, spec in grids:
        path = os.path.join(tmp, "t.root")
        build_table(spec, path)
        env = dict(os.environ, HAPRAD_SI_TABLE=path)
        worst = {k: 0.0 for k in KEYS}
        for p, r in zip(POINTS, ref):
            c = at("rc_point", p, env)
            for k in KEYS:
                if k in c and k in r and r[k] != 0:
                    worst[k] = max(worst[k], abs(c[k] - r[k]) / abs(r[k]))
        rows.append((label, worst))
        print("  done: %s" % label)

    print()
    print("max |C++ - FORTRAN| / |FORTRAN| over %d RG-E points, same model both sides" % len(POINTS))
    print("  %-26s" % "grid (Q2 x nu x z x pt2)" + "".join("%11s" % k for k in KEYS))
    for label, w in rows:
        print("  %-26s" % label + "".join("%11.2e" % w[k] for k in KEYS))
    print()
    print("sib is pure table interpolation (the codes agree to 1.6e-6 at cell centres).")
    print("f1/f2 also carry the H11 floor -- exp(delta_inf) vs the linear form --")
    print("measured at 0.2-0.4% of sigma_obs/sigma_Born at these points.")


if __name__ == "__main__":
    main()
