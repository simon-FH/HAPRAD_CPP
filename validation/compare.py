#!/usr/bin/env python3
"""Validation harness for HAPRAD_CPP.

Runs kinematic points through the C++ driver (bin/rc_point) and the HAPRAD
2.0 Fortran driver (bin/haprad2_point), and checks them in four tiers, which
fail differently and in a useful order:

  Tier 0  self-consistency of the C++ code alone. Properties that must hold
          under ANY structure-function model. Checks that need a non-zero
          sigma_Born are skipped when no table is given (--table).

  Tier 1  kinematics vs the original HAPRAD 2.0 Fortran. Lorentz invariants and
          hadron kinematics are model-independent, so the two codes should agree
          to machine precision. A failure here invalidates everything below it.

  Tier 2  the exclusive radiative tail vs HAPRAD 2.0. It is built from the MAID
          exclusive model (TSffun) and the theta matrix only, and never touches
          the semi-inclusive structure functions, so it must agree whatever
          table is used. Any disagreement here is a porting or integration
          defect.

  Tier 3  cross sections and RC factors vs HAPRAD 2.0. These agree only when the
          C++ is given a table of the Fortran's own model (make_model_table),
          and even then differ at the third digit because the C++ exponentiates
          delta_inf where the Fortran does not.

With --table, a table-coverage summary is also printed: how often the
integrals' structure-function lookups fell outside the table, into
second-neighbour cells, or into empty cells.

Exit status is 0 only if Tiers 0, 1 and 2 pass.
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
EXC_TOL = 1e-3    # exclusive tail: same model and same integrand both sides
RES_TOL = 5e-2    # cross sections: exp(delta_inf) vs the Fortran's linear form


# Structure-function table handed to the C++ via HAPRAD_SI_TABLE (--table).
TABLE = None


def run(binary, point, extra=()):
    """Run a point driver from the data directory and return its output as
    {block: {key: float}}, parsed from the `### BEGIN <block>` / `### END`
    sections and their `@key value` lines; None on timeout. Only the C++
    driver is given the table (through HAPRAD_SI_TABLE)."""
    cmd = [binary] + ["%.10g" % v for v in point] + list(extra)
    env = dict(os.environ)
    env.pop("HAPRAD_SI_TABLE", None)
    if TABLE and binary == CPP:
        env["HAPRAD_SI_TABLE"] = TABLE
    try:
        out = subprocess.run(cmd, capture_output=True, text=True,
                             timeout=1800, cwd=DATADIR, env=env).stdout
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
    """|a - b| relative to the larger of |a|, |b| (0 if both are 0)."""
    scale = max(abs(a), abs(b))
    return 0.0 if scale == 0 else abs(a - b) / scale


def model_dead(res):
    """True if the C++ result has no Born cross section (zero or not finite),
    as happens without a structure-function table."""
    sib = res.get("sib", 0.0)
    return (not math.isfinite(sib)) or abs(sib) < 1e-100


def fmt(v):
    """A number for the tables, or n/a for None."""
    return "  n/a  " if v is None else ("%9.3g" % v)


# --------------------------------------------------------------------- tier 0
def tier0(points):
    """Self-consistency checks of the C++ alone; returns (name, ok, note) per
    check, ok None for a skipped one."""
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
        results.append(("phi symmetry  sigma(60) == sigma(300)", None, "no sigma_Born (no table)"))

    # (b) the Born cross section MUST depend on z and on pt. With A/Ac/Acc held
    #     constant every h_i carries a leading tldZ that cancels the 1/tldZ in
    #     coef, so a flat model gives a z-independent sigma_Born; this check
    #     catches such a model.
    for label, p1, p2 in (
        ("sigma_Born depends on z ", (E, x, Q2, 0.35, pt, 60.0), (E, x, Q2, 0.65, pt, 60.0)),
        ("sigma_Born depends on pt", (E, x, Q2, z, 0.20, 60.0), (E, x, Q2, z, 0.60, 60.0)),
    ):
        r1, r2 = run(CPP, p1), run(CPP, p2)
        if r1 and r2 and not model_dead(r1["RES"]) and not model_dead(r2["RES"]):
            d = reldiff(r1["RES"]["sib"], r2["RES"]["sib"])
            results.append((label, d > 1e-6, "reldiff %.2e" % d))
        else:
            results.append((label, None, "no sigma_Born (no table)"))

    # (c) no spikes in phi. The Born and tail cross sections vary smoothly with
    #     phi_h, so a large jump between neighbouring phi is a numerical
    #     pathology rather than physics. This catches the phi=180 divergence
    #     of 1/tldPt2 in the A_cc -> H4z inversion, which blows up unless the
    #     azimuthal amplitudes carry their proper p_t dependence.
    #     Needs a non-zero sigma_Born -- without one the tail is noise.
    scan = []
    for ph in (170.0, 175.0, 179.0, 180.0, 181.0, 185.0):
        r = run(CPP, (E, x, Q2, z, pt, ph))
        scan.append((ph, r, (r or {}).get("RES", {}).get("tai_in", float("nan"))))
    if any(r is None or model_dead(r["RES"]) for _, r, _ in scan):
        results.append(("no spikes in phi (tail)", None, "no sigma_Born (no table)"))
    else:
        jump, at = 1.0, None
        for (p1, _, v1), (p2, _, v2) in zip(scan, scan[1:]):
            if v1 and v2 and math.isfinite(v1) and math.isfinite(v2):
                r_ = max(abs(v1 / v2), abs(v2 / v1))
                if r_ > jump:
                    jump, at = r_, "%g->%g" % (p1, p2)
        results.append(("no spikes in phi (tail)", jump < 5.0,
                        "max adjacent ratio %.1f at phi %s" % (jump, at)))

    # (d) nothing may be NaN or infinite. Split by what a missing table
    #     invalidates: KIN and tai_ex are model-independent and must always be
    #     finite; sib/sig_obs/tai_in/f1..f3 depend on the structure functions,
    #     and are reported separately.
    MODEL_DEP = {"sib", "sig_obs", "tai_in", "f1", "f2", "f3"}
    bad, bad_model = [], []
    for p in points:
        r = run(CPP, p)
        if not r:
            bad.append("%s TIMEOUT" % (p,))
            continue
        dead = model_dead(r.get("RES", {}))
        for blk in ("KIN", "RES"):
            for k, v in r.get(blk, {}).items():
                if math.isfinite(v):
                    continue
                if k in MODEL_DEP:
                    bad_model.append("phi=%g %s=%s" % (p[5], k, v))
                else:
                    bad.append("phi=%g %s=%s" % (p[5], k, v))
    note = "OK" if not bad else "; ".join(bad[:3])
    if not bad and bad_model:
        note = "OK (model-independent); %d non-finite in sib/f*" % len(bad_model)
    results.append(("all outputs finite over grid", not bad, note))

    for name, ok, note in results:
        tag = "SKIP" if ok is None else ("PASS" if ok else "FAIL")
        print("  [%s] %-38s %s" % (tag, name, note))
    return results


# --------------------------------------------------------------------- tier 1
def tier_compare(points, block, tol, title):
    """Compare every quantity of output block `block` between the C++ and the
    Fortran over the points; print the worst relative difference of each and
    return True if all are within tol. Points either code rejects are listed
    and skipped."""
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

    tai_ex is independent of the semi-inclusive model, so unlike Tier 3 it
    must agree whatever table the C++ is given.
    """
    print()
    print("=" * 78)
    print("TIER 2  exclusive radiative tail: C++ vs HAPRAD 2.0")
    print("        (independent of the structure-function table)")
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


def table_coverage(points):
    """How often the integrals looked the structure functions up outside the
    table, into second-neighbour cells, into empty cells, or at unphysical
    kinematics; per point and in total."""
    print()
    print("=" * 78)
    print("TABLE COVERAGE  %s" % TABLE)
    print("=" * 78)
    print("  %-40s %12s %9s %9s %9s %9s" % ("point (x,Q2,z,pt,phi)", "lookups", "outside", "2nd nbr", "empty", "unphys"))
    print("  " + "-" * 82)
    tot = {"sf_lookups": 0, "sf_oor": 0, "sf_second": 0, "sf_empty": 0, "sf_unphys": 0}
    for p in points:
        r = run(CPP, p)
        res = (r or {}).get("RES", {})
        n = res.get("sf_lookups", 0)
        if not n:
            continue
        for k in tot:
            tot[k] += res.get(k, 0)
        print("  x=%-5.3g Q2=%-4.3g z=%-4.3g pt=%-5.3g phi=%-5.4g %12d %8.2f%% %8.2f%% %8.2f%% %8.2f%%"
              % (p[1], p[2], p[3], p[4], p[5], n, 100. * res.get("sf_oor", 0) / n, 100. * res.get("sf_second", 0) / n,
                 100. * res.get("sf_empty", 0) / n, 100. * res.get("sf_unphys", 0) / n))
    if tot["sf_lookups"]:
        n = tot["sf_lookups"]
        print("  " + "-" * 82)
        print("  %-40s %12d %8.2f%% %8.2f%% %8.2f%% %8.2f%%" % ("all points", n, 100. * tot["sf_oor"] / n,
                                                               100. * tot["sf_second"] / n, 100. * tot["sf_empty"] / n,
                                                               100. * tot["sf_unphys"] / n))


def main():
    """Parse the options, run the tiers on the chosen grid of points and print
    the summary. Exit status 0 only if Tiers 0-2 pass."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", choices=sorted(GRIDS), default="rge")
    ap.add_argument("--tier0-only", action="store_true")
    ap.add_argument("--table", help="structure-function table for the C++ "
                    "(TSemiInclusiveModel format); without one sigma_Born is zero")
    args = ap.parse_args()

    global TABLE
    if args.table:
        TABLE = os.path.abspath(args.table)
        if not os.path.exists(TABLE):
            sys.exit("no such table: %s" % TABLE)

    if not os.path.exists(CPP):
        sys.exit("missing %s -- run `make` first" % CPP)
    # ---- Tier 0 always; the FORTRAN tiers only if haprad2_point was built ----
    points = GRIDS[args.grid]
    print("grid '%s' : %d point(s)\n" % (args.grid, len(points)))

    t0 = tier0(points)
    hard_fail = any(ok is False for _, ok, _ in t0)
    if TABLE:
        table_coverage(points)

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

    # ---- Summary ----
    print()
    print("=" * 78)
    print("  Tier 0 : %s" % ("FAIL" if hard_fail else "pass"))
    print("  Tier 1 : %s   <-- must hold; everything else rests on it"
          % ("pass" if t1 else "FAIL"))
    print("  Tier 2 : %s   <-- must hold for any table"
          % ("pass" if t2 else "FAIL"))
    print("  Tier 3 : %s   (fails unless both codes see the same model -- make_model_table)"
          % ("pass" if t3 else "FAIL"))
    print("=" * 78)
    return 0 if (t1 and t2 and not hard_fail) else 1


if __name__ == "__main__":
    sys.exit(main())
