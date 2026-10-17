#!/usr/bin/env python3
"""Figure: the radiative correction as a function of phi_h, for presentations.

Runs bin/rc_point (the C++ HAPRAD) at fixed Q2, z, p_t for several y over
phi_h = 0..180 deg, with a structure-function table given by --table
(default: HAPRAD 2.0's own PDF x FF model, the table make_model_table builds),
and draws two panels:

  left   eta = sigma_obs / sigma_Born (no exclusive tail) against phi_h, one
         curve per y;
  right  the parts of eta - 1 at the middle y: virtual + soft-photon factor,
         inelastic radiative tail, exclusive radiative tail (times Z/A = 0.5).

    ./plot_rc_phi.py --table model_table.root [--out figures/rc_vs_phi]

Writes <out>.pdf, <out>.png and <out>.tsv (the numbers plotted).
"""

import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
RC_POINT = os.path.join(HERE, "bin", "rc_point")
DATA = os.path.join(HERE, "data")

E = 10.5473            # RG-E beam energy [GeV]
M = 0.938272           # proton mass [GeV]
Q2, Z, PT = 2.5, 0.4, 0.4
YS = [0.3, 0.5, 0.65, 0.8]
Y_PARTS = 0.5          # the y whose decomposition is drawn
PHIS = list(range(0, 181, 10))
NAZ = 0.5              # Z/A of deuterium and carbon

# Reference palette: blue ramp steps 250/400/550/700 for the ordered y values,
# categorical slots 1-3 for the three parts; chart chrome and ink.
Y_COLORS = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
PART_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"


def point(table, y, phi):
    """One HAPRAD evaluation: (y, phi, sigma_Born, r_vr, r_in, r_ex)."""
    x = Q2 / (2. * M * E * y)
    env = dict(os.environ, HAPRAD_SI_TABLE=table)
    out = subprocess.run([RC_POINT] + ["%.10g" % v for v in (E, x, Q2, Z, PT, phi)], capture_output=True,
                         text=True, cwd=DATA, env=env, timeout=600).stdout
    vals, block = {}, None
    for line in out.splitlines():
        if line.startswith("### BEGIN "):
            block = line[10:].strip()
        elif line.startswith("@") and block == "RES":
            k, _, v = line[1:].partition(" ")
            vals[k.strip()] = float(v)
    sib = vals.get("sib", 0.)
    if not vals.get("kin_ok") or sib <= 0.:
        return (y, phi, 0., float("nan"), float("nan"), float("nan"))
    return (y, phi, sib, vals["sig_obs"] / sib, vals["tai_in"] / sib, vals["tai_ex"] / sib)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="structure-function table (TSemiInclusiveModel format)")
    ap.add_argument("--out", default=os.path.join(HERE, "figures", "rc_vs_phi"))
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--note", default="Structure functions: HAPRAD 2.0's PDF x FF model (GRV94 LO x PKH), tabulated. "
                    "\u03b7 without the exclusive tail; Z/A = 0.5.", help="footnote naming the table's source")
    args = ap.parse_args()
    table = os.path.realpath(args.table)
    os.makedirs(os.path.dirname(os.path.realpath(args.out)), exist_ok=True)

    # ---- HAPRAD at every (y, phi) ----
    with ThreadPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(lambda a: point(table, *a), [(y, p) for y in YS for p in PHIS]))
    with open(args.out + ".tsv", "w") as f:
        f.write("# E %g  Q2 %g  z %g  pt %g  table %s\n# y\tphi\tsigma_Born\tr_vr\tr_in\tr_ex\n" % (E, Q2, Z, PT, table))
        for r in rows:
            f.write("%g\t%g\t%.6g\t%.6g\t%.6g\t%.6g\n" % r)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.facecolor": SURFACE,
                         "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    for ax in (a1, a2):
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_xlim(0, 180)
        ax.set_xticks(range(0, 181, 30))
        ax.set_xlabel(r"$\phi_h$ [deg]")

    # ---- Left: eta(phi) per y ----
    for y, c in zip(YS, Y_COLORS):
        pts = [r for r in rows if r[0] == y and r[2] > 0]
        phi = [r[1] for r in pts]
        eta = [r[3] + r[4] for r in pts]
        a1.plot(phi, eta, color=c, linewidth=2, marker="o", markersize=4, label="y = %.2f  (x = %.3f)" % (y, Q2 / (2 * M * E * y)))
        if eta:
            a1.annotate("y = %.2f" % y, (phi[-1], eta[-1]), xytext=(6, 0), textcoords="offset points", va="center",
                        color=INK2, fontsize=10, annotation_clip=False)
    a1.axhline(1.0, color=AXIS, linewidth=1)
    a1.set_ylabel(r"$\eta = \sigma_{\mathrm{obs}} / \sigma_{\mathrm{Born}}$")
    a1.set_title("Radiative correction factor", color=INK, loc="left", fontsize=12)
    # Room below the lowest curve for the legend.
    lo = min(r[3] + r[4] for r in rows if r[2] > 0)
    a1.set_ylim(lo - 0.05, a1.get_ylim()[1])
    a1.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK2)

    # ---- Right: the parts of eta - 1 at y = Y_PARTS ----
    pts = [r for r in rows if r[0] == Y_PARTS and r[2] > 0]
    phi = [r[1] for r in pts]
    parts = [("virtual + soft photons", [r[3] - 1. for r in pts]), ("inelastic tail", [r[4] for r in pts]),
             ("exclusive tail $\\times$ Z/A", [NAZ * r[5] for r in pts])]
    for (label, v), c in zip(parts, PART_COLORS):
        a2.plot(phi, v, color=c, linewidth=2, marker="o", markersize=4, label=label)
    a2.axhline(0.0, color=AXIS, linewidth=1)
    a2.set_ylabel(r"contribution to $\eta - 1$")
    a2.set_title("Its parts at y = %.2f" % Y_PARTS, color=INK, loc="left", fontsize=12)
    a2.legend(frameon=False, fontsize=9, loc="best", labelcolor=INK2)

    fig.suptitle(r"HAPRAD, $\pi^+$ SIDIS at E = %.4g GeV:  Q$^2$ = %g GeV$^2$, z = %g, p$_t$ = %g GeV" % (E, Q2, Z, PT),
                 color=INK, fontsize=12, x=0.01, ha="left")
    fig.get_layout_engine().set(rect=(0, 0.05, 1, 0.95))
    fig.text(0.01, 0.01, args.note, color=MUTED, fontsize=9, ha="left", va="bottom")
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (args.out, ext), dpi=200)
    print("wrote %s.{pdf,png,tsv}" % args.out)


if __name__ == "__main__":
    main()
