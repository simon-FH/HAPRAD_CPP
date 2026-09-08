#!/usr/bin/env python3
"""Validation harness for HAPRAD_CPP.

Three layers, because they fail differently and in a useful order:

  Tier 0  self-consistency of the C++ code alone. Properties that must hold
          under ANY structure-function model, so they stay meaningful while the
          model is still missing (hurdle H1).

  Tier 1  kinematics vs the original HAPRAD 2.0 Fortran. Lorentz invariants and
          hadron kinematics are model-independent, so the two codes should agree
          to machine precision. A failure here invalidates everything below it.

  Tier 2  the exclusive radiative tail vs HAPRAD 2.0. This one is comparable
          RIGHT NOW: it is built from the MAID exclusive model (TSffun) and the
          theta matrix only, and never touches the semi-inclusive structure
          functions, so H1 does not affect it. Any disagreement here is a real
          porting or integration defect.

  Tier 3  cross sections and RC factors vs HAPRAD 2.0. These will NOT agree
          while H1 is open, and even afterwards will differ at the third digit
          because the C++ exponentiates delta_inf where the Fortran does not
          (hurdle H11).

Exit status is 0 only if every check that is currently expected to pass does.
"""

import argparse
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CPP = os.path.join(HERE, "bin", "rc_point")
F77 = os.path.join(HERE, "bin", "haprad2_point")
# Both drivers open pi_n_maid.dat and the *.grid files from the cwd.
DATADIR = os.path.join(HERE, "data")

# ---------------------------------------------------------------- point grids
# (E, x, Q2, z, pt, phi_deg)
GRIDS = {
    # The point rcdat.f ships with, plus neighbours: 6 GeV CLAS/EG2-era.
    "eg2": [
        (6.0000, 0.32, 2.5, 0.50, 0.35,  60.0),
        (6.0000, 0.32, 2.5, 0.50, 0.35, 140.0),
        (6.0000, 0.24, 1.5, 0.40, 0.30,  60.0),
        (5.0150, 0.24, 1.5, 0.50, 0.40,  60.0),
    ],
    # RG-E: beam energy and centroids from run 020026 pi+ DIS sample.
    "rge": [
        (10.5473, 0.234, 2.5, 0.34, 0.424,  60.0),
        (10.5473, 0.234, 2.5, 0.34, 0.424, 180.0),
        (10.5473, 0.150, 2.0, 0.30, 0.350,  90.0),
        (10.5473, 0.300, 4.0, 0.50, 0.500,  30.0),
        (10.5473, 0.234, 2.5, 0.70, 0.300,  60.0),
    ],
}

KIN_TOL = 1e-10   # kinematics: same formulas, both double precision
EXC_TOL = 5e-2    # exclusive tail: same model both sides, should be close
RES_TOL = 5e-2    # cross sections: see H11


def run(binary, point, extra=()):
    """Run a point driver and return {block: {key: float}}."""
    cmd = [binary] + ["%.10g" % v for v in point] + list(extra)
    try:
        out = subprocess.run(cmd, capture_output=True, text=True,
                             timeout=1800, cwd=DATADIR).stdout
    except subprocess.TimeoutExpired:
        return None
    blocks, cur = {}, None
    for line in out.splitlines():
        if line.startswith("### BEGIN "):
            cur = line[10:].strip()
            blocks[cur] = {}
        elif line.startswith("### END "):
            cur = None
        elif cur is not None and line.startswith("@"):
            name, _, val = line[1:].partition(" ")
            try:
                blocks[cur][name.strip()] = float(val)
            except ValueError:
                pass
    return blocks


def reldiff(a, b):
    scale = max(abs(a), abs(b))
    return 0.0 if scale == 0 else abs(a - b) / scale


def model_dead(res):
    """True if the C++ semi-inclusive model is uninitialised (hurdle H1)."""
    sib = res.get("sib", 0.0)
    return (not math.isfinite(sib)) or abs(sib) < 1e-100


def fmt(v):
    return "  n/a  " if v is None else ("%9.3g" % v)


# --------------------------------------------------------------------- tier 0
def tier0(points):
    print("=" * 78)
    print("TIER 0  self-consistency of the C++ code (model-independent)")
    print("=" * 78)
    results = []

    E, x, Q2, z, pt, _ = points[0]

    # (a) phi symmetry: the Born cross section depends on phi only through
    #     cos(phi) and cos(2 phi), so sigma(phi) == sigma(360 - phi) exactly.
    a = run(CPP, (E, x, Q2, z, pt, 60.0))
    b = run(CPP, (E, x, Q2, z, pt, 300.0))
    if a and b and not model_dead(a["RES"]) and not model_dead(b["RES"]):
        d = reldiff(a["RES"]["sib"], b["RES"]["sib"])
        results.append(("phi symmetry  sigma(60) == sigma(300)", d < 1e-8, "reldiff %.2e" % d))
    else:
        results.append(("phi symmetry  sigma(60) == sigma(300)", None, "model uninitialised (H1)"))

    # (b) the Born cross section MUST depend on z and on pt. With A/Ac/Acc held
    #     constant every h_i carries a leading tldZ that cancels the 1/tldZ in
    #     coef, so a flat model gives a z-independent sigma_Born -- which is the
    #     signature of H1 and the reason this check exists.
    for label, p1, p2 in (
        ("sigma_Born depends on z ", (E, x, Q2, 0.35, pt, 60.0), (E, x, Q2, 0.65, pt, 60.0)),
        ("sigma_Born depends on pt", (E, x, Q2, z, 0.20, 60.0), (E, x, Q2, z, 0.60, 60.0)),
    ):
        r1, r2 = run(CPP, p1), run(CPP, p2)
        if r1 and r2 and not model_dead(r1["RES"]) and not model_dead(r2["RES"]):
            d = reldiff(r1["RES"]["sib"], r2["RES"]["sib"])
            results.append((label, d > 1e-6, "reldiff %.2e" % d))
        else:
            results.append((label, None, "model uninitialised (H1)"))

    # (c) nothing may be NaN or infinite anywhere on the grid (hurdle H4).
    bad = []
    for p in points:
        r = run(CPP, p)
        if not r:
            bad.append("%s TIMEOUT" % (p,))
            continue
        for blk in ("KIN", "RES"):
            for k, v in r.get(blk, {}).items():
                if not math.isfinite(v):
                    bad.append("phi=%g %s=%s" % (p[5], k, v))
    results.append(("all outputs finite over grid", not bad,
                    "OK" if not bad else "; ".join(bad[:3])))

    for name, ok, note in results:
        tag = "SKIP" if ok is None else ("PASS" if ok else "FAIL")
        print("  [%s] %-38s %s" % (tag, name, note))
    return results


# --------------------------------------------------------------------- tier 1
def tier_compare(points, block, tol, title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)

    worst = {}
    n_ok = 0
    for p in points:
        c, f = run(CPP, p), run(F77, p)
        if not c or not f:
            print("  point %s: driver failed" % (p,))
            continue
        cb, fb = c.get(block, {}), f.get(block, {})
        if not cb.get("kin_ok", 0) or not fb.get("kin_ok", 0):
            print("  point E=%.4g x=%.3g Q2=%.3g z=%.3g pt=%.3g phi=%g : "
                  "rejected (cpp=%g f77=%g)"
                  % (p + (cb.get("kin_ok", -1), fb.get("kin_ok", -1))))
            continue
        n_ok += 1
        for k in sorted(set(cb) & set(fb)):
            if k == "kin_ok":
                continue
            d = reldiff(cb[k], fb[k])
            if k not in worst or d > worst[k][0]:
                worst[k] = (d, cb[k], fb[k], p)

    if not worst:
        print("  no comparable points")
        return False

    print("  %-10s %-11s %-24s %-24s" % ("quantity", "max reldiff", "C++", "HAPRAD 2.0"))
    print("  " + "-" * 72)
    allpass = True
    for k, (d, cv, fv, _) in sorted(worst.items(), key=lambda kv: -kv[1][0]):
        ok = d <= tol
        allpass &= ok
        print("  %-10s %-11.3e % -24.15g % -24.15g %s"
              % (k, d, cv, fv, "" if ok else "  <-- exceeds %.0e" % tol))
    print("  " + "-" * 72)
    print("  %d point(s) compared, tolerance %.0e : %s"
          % (n_ok, tol, "PASS" if allpass else "FAIL"))
    return allpass


def tier2(points):
    """Compare only the exclusive radiative tail.

    tai_ex is independent of the (currently missing) semi-inclusive model, so
    unlike Tier 3 this is expected to agree today.
    """
    print()
    print("=" * 78)
    print("TIER 2  exclusive radiative tail: C++ vs HAPRAD 2.0")
    print("        (independent of H1 -- should agree NOW)")
    print("=" * 78)

    rows, allpass = [], True
    for p in points:
        c, f = run(CPP, p), run(F77, p)
        if not c or not f:
            continue
        cr, fr = c.get("RES", {}), f.get("RES", {})
        if not cr.get("kin_ok", 0) or not fr.get("kin_ok", 0):
            continue
        cv, fv = cr.get("tai_ex"), fr.get("tai_ex")
        if cv is None or fv is None:
            continue
        d = reldiff(cv, fv)
        allpass &= d <= EXC_TOL
        rows.append((p, cv, fv, d))

    if not rows:
        print("  no comparable points")
        return False

    print("  %-34s %-14s %-14s %-10s" % ("point (x,Q2,z,pt,phi)", "C++", "HAPRAD 2.0", "reldiff"))
    print("  " + "-" * 72)
    for p, cv, fv, d in rows:
        print("  x=%-5.3g Q2=%-4.3g z=%-4.3g pt=%-5.3g phi=%-5.4g % -13.6g % -13.6g %-9.3e %s"
              % (p[1], p[2], p[3], p[4], p[5], cv, fv, d,
                 "" if d <= EXC_TOL else "<-- FAIL"))
    print("  " + "-" * 72)
    print("  tolerance %.0e : %s" % (EXC_TOL, "PASS" if allpass else "FAIL"))
    return allpass


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", choices=sorted(GRIDS), default="rge")
    ap.add_argument("--tier0-only", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(CPP):
        sys.exit("missing %s -- run `make` first" % CPP)
    points = GRIDS[args.grid]
    print("grid '%s' : %d point(s)\n" % (args.grid, len(points)))

    t0 = tier0(points)
    hard_fail = any(ok is False for _, ok, _ in t0)

    if args.tier0_only:
        return 1 if hard_fail else 0

    if not os.path.exists(F77):
        print("\nHAPRAD 2.0 reference not built (CERN_LIB unset) -- "
              "skipping Tier 1 and Tier 3.")
        return 1 if hard_fail else 0

    t1 = tier_compare(points, "KIN", KIN_TOL,
                      "TIER 1  kinematics: C++ vs HAPRAD 2.0 Fortran")
    t2 = tier2(points)
    t3 = tier_compare(points, "RES", RES_TOL,
                      "TIER 3  cross sections: C++ vs HAPRAD 2.0 Fortran")

    print()
    print("=" * 78)
    print("  Tier 0 : %s" % ("FAIL" if hard_fail else "pass"))
    print("  Tier 1 : %s   <-- must hold; everything else rests on it"
          % ("pass" if t1 else "FAIL"))
    print("  Tier 2 : %s   <-- should hold today; H1 does not excuse it"
          % ("pass" if t2 else "FAIL"))
    print("  Tier 3 : %s   (expected to fail while H1 is open)"
          % ("pass" if t3 else "FAIL"))
    print("=" * 78)
    return 0 if (t1 and t2 and not hard_fail) else 1


if __name__ == "__main__":
    sys.exit(main())
