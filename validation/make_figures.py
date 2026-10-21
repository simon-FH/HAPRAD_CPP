#!/usr/bin/env python3
"""Figures for presentations: what HAPRAD's radiative correction does, how the
C++ port and the data-driven method were validated, and the RG-E LD2+C results.

    ./make_figures.py --table <model_table.root> [--only fig02,fig05] [--recompute]

Each figure is written to figures/<name>.pdf and .png, with the numbers it
shows in figures/<name>.tsv. Inputs:

  bin/rc_point, bin/haprad2_point   HAPRAD (C++ and FORTRAN) at single points:
                                    kinematic scans and the port comparison,
                                    cached in figures/data/scan_*.tsv
  --table                           structure-function table for the C++: the
                                    one make_model_table builds from HAPRAD
                                    2.0's own PDF x FF model, so the scans show
                                    physical corrections and the comparison
                                    has the same physics on both sides
  figures/data/closure_*.txt        outputs of closure_normalisation.py and
                                    closure_iteration.py
  figures/data/shapes_*.tsv         table_shapes.C on a data table and on the
                                    model table
  figures/data/farm_*.json          the full LD2+C farm run, summarised
  figures/data/broad_comparison.tsv C++ and FORTRAN at 504 grid points and 17
                                    nodes (load_broad; a few minutes to recompute)

Figures:
  fig01_rc_vs_phi        eta against phi_h at four y, and its parts
  fig02_rc_kinematics    phi-averaged eta and its range over phi, against y, z, p_t
  fig03_rc_parts         the parts of eta - 1 against y, z, p_t
  fig04_fake_modulations the change of the cos(phi), cos(2 phi) amplitudes
  fig05_port_validation  C++ against HAPRAD 2.0 FORTRAN, point by point
  fig06_normalisation    the absolute-normalisation closure of the table
  fig07_convergence      the RC iteration: toy closure against real data
  fig08_data_rc          RC factors of the selected LD2+C events, every iteration
  fig09_table_shapes     why: phi shapes of the data table's cells
  fig10_farm_cost        CPU and wall time of the farm run
  fig11_where_rc_gt2     where in the kinematics events get a factor above 2
  fig12_coverage         table cells fitted / filled / empty, extreme grid nodes
  fig13_data_vs_model    the data-driven eta(phi) against the model's
  fig14_output_phi       the RC factor against phi_h at 12 kinematic points, every iteration
  fig15_output_maps      the phi-averaged RC factor over Q2 and nu, every iteration
  fig16_output_parts     its two parts (virtual + soft, inelastic tail), every iteration
  fig17_output_table     the weights' effect on the table: yields and modulations
  fig18_output_cells     single table cells' phi shapes, every iteration
  fig19_port_broad       C++ against FORTRAN at 504 RG-E points and the full run's extreme nodes
  fig20_rc_phi           pages: eta against phi, one page per (Q2, nu) node, panels z x p_t
  fig21_xsec_phi         pages: sigma_obs and sigma_Born against phi, same layout
  fig22_rc_pt2           pages: phi-integrated eta against p_t^2, one page per Q2
  fig23_rc_z             pages: phi-integrated eta against z, one page per Q2
  fig24_convergence_tail distribution of convergence: the bulk of the grid converges, a tail runs away

fig07-09 and fig11-18 read the full LD2+C run from figures/farm/tsv/, made by
export_results.sh on the farm and export_farm.C here. fig20-23 overlay the
grid of the HAPRAD 2.0 model table, figures/farm/tsv/nodes_model.tsv: built by
Utilities/run_rc_grid.sh with tables/model_table_E10.5473.root, exported with
export_farm.C's export_grid. Their pages go to figures/<name>.pdf (all pages)
and figures/<name>/*.png.
"""

import argparse
import json
import math
import os
import re
import subprocess
import textwrap
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")
DATA = os.path.join(HERE, "data")             # MAID and PKH grids for the drivers
FIG = os.path.join(HERE, "figures")
FIGDATA = os.path.join(FIG, "data")

E = 10.5473            # RG-E beam energy [GeV]
M = 0.938272           # proton mass [GeV]
NAZ = 0.5              # Z/A of deuterium and carbon
BASE = dict(Q2=2.5, y=0.5, z=0.4, pt=0.4)
PHIS = list(range(0, 181, 15))                # 13 points: Simpson on [0, pi]
SCANS = {
    "y": [round(0.30 + 0.05 * i, 2) for i in range(12)],
    "z": [round(0.20 + 0.05 * i, 2) for i in range(13)],
    "pt": [round(0.10 + 0.10 * i, 2) for i in range(10)],
}
Y_PHI = [0.30, 0.50, 0.65, 0.80]              # the curves of fig01

# ---- Reference palette (validated): ordinal blue ramp, categorical slots, chrome ----
RAMP4 = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
RAMP2 = ["#86b6ef", "#1c5cab"]
CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BAND = "#cde2fb"                              # blue step 100: a range, behind the line


# ==== HAPRAD at single points ====

def run_driver(binary, point, table=None):
    """Run a point driver; return {block: {key: float}} from its @key lines."""
    env = dict(os.environ)
    env.pop("HAPRAD_SI_TABLE", None)
    if table:
        env["HAPRAD_SI_TABLE"] = table
    out = subprocess.run([os.path.join(BIN, binary)] + ["%.10g" % v for v in point], capture_output=True, text=True,
                         cwd=DATA, env=env, timeout=1800).stdout
    blocks, cur = {}, None
    for line in out.splitlines():
        if line.startswith("### BEGIN "):
            cur = line[10:].strip()
            blocks[cur] = {}
        elif line.startswith("### END"):
            cur = None
        elif cur and line.startswith("@"):
            k, _, v = line[1:].partition(" ")
            try:
                blocks[cur][k.strip()] = float(v.split()[0].replace("D", "E"))
            except (ValueError, IndexError):
                pass
    return blocks


def x_of(Q2, y):
    return Q2 / (2. * M * E * y)


def scan_points():
    """Every (label, Q2, y, z, pt, phi) the HAPRAD figures need."""
    pts = set()
    for var, vals in SCANS.items():
        for v in vals:
            k = dict(BASE, **{var: v})
            for phi in PHIS:
                pts.add((k["Q2"], k["y"], k["z"], k["pt"], phi))
    for y in Y_PHI:
        for phi in PHIS:
            pts.add((BASE["Q2"], y, BASE["z"], BASE["pt"], phi))
    return sorted(pts)


def load_scan(table, recompute, jobs):
    """{(Q2, y, z, pt, phi): (sib, sig_obs, tai_in, tai_ex)}, cached in scan_haprad.tsv."""
    path = os.path.join(FIGDATA, "scan_haprad.tsv")
    if os.path.exists(path) and not recompute:
        res = {}
        for line in open(path):
            if line.startswith("#"):
                continue
            f = [float(v) for v in line.split()]
            res[tuple(f[:5])] = tuple(f[5:])
        return res
    pts = scan_points()

    def one(p):
        Q2, y, z, pt, phi = p
        r = run_driver("rc_point", (E, x_of(Q2, y), Q2, z, pt, phi), table).get("RES", {})
        ok = r.get("kin_ok", 0) and r.get("sib", 0) > 0
        return p, (r["sib"], r["sig_obs"], r["tai_in"], r["tai_ex"]) if ok else (0., 0., 0., 0.)

    with ThreadPoolExecutor(jobs) as pool:
        res = dict(pool.map(one, pts))
    with open(path, "w") as out:
        out.write("# HAPRAD (C++) with %s at E = %g GeV\n# Q2 y z pt phi sigma_Born sig_obs tai_in tai_ex\n" % (table, E))
        for p in pts:
            out.write("%g %g %g %g %g %.8g %.8g %.8g %.8g\n" % (p + res[p]))
    return res


def simpson(f):
    """Integral over phi in [0, pi] of samples f at PHIS (equal steps, even intervals)."""
    h = math.radians(PHIS[1] - PHIS[0])
    return h / 3. * (f[0] + f[-1] + 4. * sum(f[1:-1:2]) + 2. * sum(f[2:-1:2]))


def setting(scan, Q2, y, z, pt):
    """The phi-integrated view of one kinematic setting, or None if HAPRAD
    rejected any of its points. eta excludes the exclusive tail (as the
    weights do); amplitudes are 2 <cos n phi> of the cross section."""
    rows = [scan[(Q2, y, z, pt, phi)] for phi in PHIS]
    if any(r[0] <= 0. for r in rows):
        return None
    sb = [r[0] for r in rows]
    so = [r[1] + r[2] for r in rows]
    cos = [[math.cos(n * math.radians(p)) for p in PHIS] for n in (1, 2)]
    ib, io = simpson(sb), simpson(so)
    amp = lambda s, i, n: 2. * simpson([a * c for a, c in zip(s, cos[n - 1])]) / i
    eta = [o / b for o, b in zip(so, sb)]
    return dict(eta_avg=io / ib, eta_min=min(eta), eta_max=max(eta), eta=eta,
                vs=(simpson([r[1] for r in rows]) - ib) / ib, inel=simpson([r[2] for r in rows]) / ib,
                excl=NAZ * simpson([r[3] for r in rows]) / ib,
                dA1=amp(so, io, 1) - amp(sb, ib, 1), dA2=amp(so, io, 2) - amp(sb, ib, 2),
                A1=amp(sb, ib, 1), A2=amp(sb, ib, 2))


# ==== Drawing helpers ====

def mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.facecolor": SURFACE,
                         "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE, "legend.frameon": False,
                         "legend.fontsize": 9, "legend.labelcolor": INK2, "lines.linewidth": 2, "lines.markersize": 4})
    return plt


def style(ax, title=None):
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if title:
        ax.set_title(title, color=INK, loc="left", fontsize=12)


def save(fig, name, suptitle, note, table_header, table_rows):
    """Title, footnote, PDF + PNG, and the plotted numbers as TSV."""
    fig.suptitle(suptitle, color=INK, fontsize=12, x=0.01, ha="left")
    lines = textwrap.wrap(note, int(fig.get_figwidth() * 15.5))
    bottom = (0.08 + 0.15 * len(lines)) / fig.get_figheight()
    fig.get_layout_engine().set(rect=(0, bottom, 1, 1 - bottom))
    fig.text(0.01, 0.01, "\n".join(lines), color=MUTED, fontsize=8.5, ha="left", va="bottom")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, "%s.%s" % (name, ext)), dpi=200)
    with open(os.path.join(FIG, name + ".tsv"), "w") as out:
        out.write("\t".join(table_header) + "\n")
        for r in table_rows:
            out.write("\t".join(r if isinstance(r, str) else ("%.6g" % v if isinstance(v, float) else str(v)) for v in r) + "\n")
    print("wrote", name)


KIN_TITLE = r"HAPRAD, $\pi^+$ SIDIS at E = %.4g GeV" % E
MODEL_NOTE = "Structure functions: HAPRAD 2.0's PDF x FF model (GRV94 LO x PKH), tabulated. η without the exclusive tail; Z/A = 0.5."
VAR_LABEL = {"y": "y", "z": "z", "pt": r"$p_t$ [GeV]"}


def base_text(var):
    """The fixed kinematics of a scan, for the panel title."""
    parts = {"Q2": r"Q$^2$ = %g" % BASE["Q2"], "y": "y = %g" % BASE["y"], "z": "z = %g" % BASE["z"],
             "pt": r"$p_t$ = %g" % BASE["pt"]}
    return ", ".join(v for k, v in parts.items() if k != var)


# ==== The figures ====

def fig01(scan):
    plt = mpl()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    rows = []
    for y, c in zip(Y_PHI, RAMP4):
        s = setting(scan, BASE["Q2"], y, BASE["z"], BASE["pt"])
        a1.plot(PHIS, s["eta"], color=c, marker="o", label="y = %.2f  (x = %.3f)" % (y, x_of(BASE["Q2"], y)))
        a1.annotate("y = %.2f" % y, (PHIS[-1], s["eta"][-1]), xytext=(6, 0), textcoords="offset points", va="center",
                    color=INK2, fontsize=10, annotation_clip=False)
        rows += [(y, float(p), e) for p, e in zip(PHIS, s["eta"])]
    a1.axhline(1., color=AXIS, linewidth=1)
    lo = min(r[2] for r in rows)
    a1.set_ylim(lo - 0.05, a1.get_ylim()[1])
    a1.legend(loc="lower right")
    a1.set_ylabel(r"$\eta = \sigma_{\mathrm{obs}} / \sigma_{\mathrm{Born}}$")
    style(a1, "Radiative correction factor")
    pts = [scan[(BASE["Q2"], BASE["y"], BASE["z"], BASE["pt"], p)] for p in PHIS]
    parts = [("virtual + soft photons", [r[1] / r[0] - 1. for r in pts]), ("inelastic tail", [r[2] / r[0] for r in pts]),
             ("exclusive tail × Z/A", [NAZ * r[3] / r[0] for r in pts])]
    for (label, v), c in zip(parts, CAT):
        a2.plot(PHIS, v, color=c, marker="o", label=label)
    a2.axhline(0., color=AXIS, linewidth=1)
    a2.legend(loc="lower right")
    a2.set_ylabel(r"contribution to $\eta - 1$")
    style(a2, "Its parts at y = %.2f" % BASE["y"])
    for ax in (a1, a2):
        ax.set_xlim(0, 180)
        ax.set_xticks(range(0, 181, 30))
        ax.set_xlabel(r"$\phi_h$ [deg]")
    save(fig, "fig01_rc_vs_phi", KIN_TITLE + r":  Q$^2$ = %g GeV$^2$, z = %g, $p_t$ = %g GeV" % (BASE["Q2"], BASE["z"], BASE["pt"]),
         MODEL_NOTE, ["y", "phi_deg", "eta"], rows)


def scan_settings(scan, var):
    out = []
    for v in SCANS[var]:
        k = dict(BASE, **{var: v})
        s = setting(scan, k["Q2"], k["y"], k["z"], k["pt"])
        if s:
            out.append((v, s))
    return out


def fig02(scan):
    plt = mpl()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True, sharey=True)
    rows = []
    for ax, var in zip(axes, ("y", "z", "pt")):
        ss = scan_settings(scan, var)
        xs = [v for v, _ in ss]
        ax.fill_between(xs, [s["eta_min"] for _, s in ss], [s["eta_max"] for _, s in ss], color=BAND, linewidth=0,
                        label=r"range over $\phi_h$")
        ax.plot(xs, [s["eta_avg"] for _, s in ss], color=CAT[0], marker="o", label=r"averaged over $\phi_h$")
        ax.axhline(1., color=AXIS, linewidth=1)
        ax.set_xlabel(VAR_LABEL[var])
        style(ax, base_text(var))
        rows += [(var, v, s["eta_avg"], s["eta_min"], s["eta_max"]) for v, s in ss]
    axes[0].set_ylabel(r"$\eta = \sigma_{\mathrm{obs}} / \sigma_{\mathrm{Born}}$")
    axes[0].legend(loc="lower right")
    save(fig, "fig02_rc_kinematics", KIN_TITLE + r": the correction factor across the kinematics (Q$^2$ = %g GeV$^2$)" % BASE["Q2"],
         MODEL_NOTE + " Average: ratio of the phi-integrated cross sections.", ["variable", "value", "eta_avg", "eta_min", "eta_max"],
         rows)


def fig03(scan):
    plt = mpl()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True, sharey=True)
    rows = []
    labels = ["virtual + soft photons", "inelastic tail", "exclusive tail × Z/A"]
    for ax, var in zip(axes, ("y", "z", "pt")):
        ss = scan_settings(scan, var)
        xs = [v for v, _ in ss]
        for key, label, c in zip(("vs", "inel", "excl"), labels, CAT):
            ax.plot(xs, [s[key] for _, s in ss], color=c, marker="o", label=label)
        ax.axhline(0., color=AXIS, linewidth=1)
        ax.set_xlabel(VAR_LABEL[var])
        style(ax, base_text(var))
        rows += [(var, v, s["vs"], s["inel"], s["excl"]) for v, s in ss]
    axes[0].set_ylabel(r"contribution to $\eta - 1$ (averaged over $\phi_h$)")
    axes[0].legend(loc="upper left")
    save(fig, "fig03_rc_parts", KIN_TITLE + ": the parts of the correction", MODEL_NOTE,
         ["variable", "value", "virtual_soft", "inelastic_tail", "exclusive_tail_x_ZA"], rows)


def fig04(scan):
    plt = mpl()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True, sharey=True)
    rows = []
    for ax, var in zip(axes, ("y", "z", "pt")):
        ss = scan_settings(scan, var)
        xs = [v for v, _ in ss]
        ax.plot(xs, [s["dA1"] for _, s in ss], color=CAT[0], marker="o", label=r"$\cos\phi_h$ amplitude")
        ax.plot(xs, [s["dA2"] for _, s in ss], color=CAT[1], marker="o", label=r"$\cos 2\phi_h$ amplitude")
        ax.axhline(0., color=AXIS, linewidth=1)
        ax.set_xlabel(VAR_LABEL[var])
        style(ax, base_text(var))
        rows += [(var, v, s["A1"], s["dA1"], s["A2"], s["dA2"]) for v, s in ss]
    axes[0].set_ylabel("observed − Born amplitude")
    axes[0].legend(loc="best")
    save(fig, "fig04_fake_modulations", KIN_TITLE + ": azimuthal modulations the radiation adds to the measured cross section",
         MODEL_NOTE + " Amplitude: 2<cos n phi> of the cross section; the .tsv also lists the Born amplitudes.",
         ["variable", "value", "A_cos_Born", "dA_cos", "A_cos2_Born", "dA_cos2"], rows)


# The points compare.py uses: RG-E kinematics, plus EG2-era ones for the
# table-independent parts (kinematics and exclusive tail).
POINTS_RGE = [(10.5473, 0.234, 2.5, 0.34, 0.424, 60.0), (10.5473, 0.234, 2.5, 0.34, 0.424, 180.0),
              (10.5473, 0.150, 2.0, 0.30, 0.350, 90.0), (10.5473, 0.300, 4.0, 0.50, 0.500, 30.0),
              (10.5473, 0.234, 2.5, 0.70, 0.300, 60.0)]
POINTS_EG2 = [(6.0, 0.32, 2.5, 0.50, 0.35, 60.0), (6.0, 0.32, 2.5, 0.50, 0.35, 140.0), (6.0, 0.24, 1.5, 0.40, 0.30, 60.0),
              (5.015, 0.24, 1.5, 0.50, 0.40, 60.0)]


def load_port(table, recompute, jobs):
    """Per point: max relative difference of the kinematics, and the C++ and
    FORTRAN exclusive tail, Born cross section and eta; cached."""
    path = os.path.join(FIGDATA, "port_comparison.tsv")
    cols = ["E", "x", "Q2", "z", "pt", "phi", "kin_maxrel", "tai_ex_cpp", "tai_ex_f77", "sib_cpp", "sib_f77", "eta_cpp", "eta_f77"]
    if os.path.exists(path) and not recompute:
        return [dict(zip(cols, map(float, l.split()))) for l in open(path) if not l.startswith("#")]

    def one(p):
        c, f = run_driver("rc_point", p, table if p[0] == E else None), run_driver("haprad2_point", p)
        ck, fk = c.get("KIN", {}), f.get("KIN", {})
        rel = [abs(ck[k] - fk[k]) / max(abs(ck[k]), abs(fk[k])) for k in set(ck) & set(fk)
               if k != "kin_ok" and max(abs(ck[k]), abs(fk[k])) > 0]
        cr, fr = c.get("RES", {}), f.get("RES", {})
        eta = lambda r: (r["sig_obs"] + r["tai_in"]) / r["sib"] if r.get("sib", 0) > 0 else float("nan")
        return dict(zip(cols, list(p) + [max(rel), cr.get("tai_ex", float("nan")), fr.get("tai_ex", float("nan")),
                                         cr.get("sib", 0.), fr.get("sib", 0.), eta(cr), eta(fr)]))

    with ThreadPoolExecutor(jobs) as pool:
        res = list(pool.map(one, POINTS_RGE + POINTS_EG2))
    with open(path, "w") as out:
        out.write("# C++ (rc_point, table %s at E = %g) vs HAPRAD 2.0 (haprad2_point)\n# %s\n" % (table, E, " ".join(cols)))
        for r in res:
            out.write(" ".join("%.10g" % r[k] for k in cols) + "\n")
    return res


def fig05(port):
    plt = mpl()
    fig, ax = plt.subplots(figsize=(11, 4.6), constrained_layout=True)
    cats = ["kinematics\n(worst quantity per point)", "exclusive tail", "Born cross section\n(model table)",
            "correction factor η\n(model table)"]
    rows = []
    for r in port:
        rge = r["E"] == E
        tag = "x=%.3g Q2=%g z=%g pt=%g phi=%g" % (r["x"], r["Q2"], r["z"], r["pt"], r["phi"])
        vals = [r["kin_maxrel"], abs(r["tai_ex_cpp"] - r["tai_ex_f77"]) / abs(r["tai_ex_f77"])]
        if rge:
            vals += [abs(r["sib_cpp"] - r["sib_f77"]) / r["sib_f77"], abs(r["eta_cpp"] - r["eta_f77"]) / r["eta_f77"]]
        for i, v in enumerate(vals):
            jitter = (POINTS_RGE + POINTS_EG2).index(tuple(r[k] for k in ("E", "x", "Q2", "z", "pt", "phi"))) * 0.05 - 0.2
            ax.scatter(i + jitter, max(v, 1e-17), s=36, color=CAT[0] if rge else CAT[1], edgecolors=SURFACE, linewidths=1.5,
                       zorder=3, label=("RG-E, E = 10.55 GeV" if rge else "EG2-era, E = 5-6 GeV") if i == 0 else None)
            rows.append((cats[i].replace("\n", " "), tag, v))
        if rge and abs(r["eta_cpp"] - r["eta_f77"]) / r["eta_f77"] > 0.015:
            ax.annotate("x = %.2f: FORTRAN\nmodel capped here" % r["x"], (3 + jitter, abs(r["eta_cpp"] - r["eta_f77"]) / r["eta_f77"]),
                        xytext=(0, 10), textcoords="offset points", ha="center", va="bottom", color=INK2, fontsize=9)
    ax.set_yscale("log")
    ax.set_ylim(1e-17, 1)
    ax.axhline(2.2e-16, color=AXIS, linewidth=1)
    ax.annotate("double precision", (3.45, 2.2e-16), xytext=(0, 3), textcoords="offset points", ha="right", color=MUTED, fontsize=9)
    ax.set_xticks(range(4))
    ax.set_xticklabels(cats, color=INK2)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylabel("|C++ − FORTRAN| / FORTRAN")
    handles, labels = ax.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    ax.legend(uniq.values(), uniq.keys(), loc="upper left")
    style(ax)
    save(fig, "fig05_port_validation", "The C++ port against the original HAPRAD 2.0 FORTRAN, point by point",
         "Same inputs to both codes. The Born cross section and η need a structure-function table on the C++ side: HAPRAD 2.0's own "
         "model, tabulated (RG-E points only). Remaining η differences: interpolation in the table, and exp(δ_inf) vs the linear form.",
         ["quantity", "point", "relative_difference"], rows)


# The broad comparison: a grid over the RG-E acceptance, and the RC-grid
# nodes where the data table of the full LD2+C run gave its most extreme
# factors (Q2, nu, z, pt, phi), both codes with HAPRAD 2.0's model.
BROAD_GRID = dict(Q2=[1.5, 2.5, 4.0, 6.0], y=[0.35, 0.5, 0.65, 0.8], z=[0.2, 0.4, 0.6], pt=[0.2, 0.5, 0.8], phi=[0, 60, 120, 180])
EXTREME_NODES = [(1, 8, 0.25, 0.25, 15), (1, 7, 0.15, 0.25, 0), (1, 9, 0.75, 0.5, 15), (1, 6, 0.15, 0.25, 0), (2, 9, 0.35, 0.25, 0),
                 (2, 7, 0.65, 0.25, 15), (1, 6, 0.15, 0.125, 90), (8, 9, 0.25, 0.25, 30), (2, 5, 0.45, 0.375, 0),
                 (2, 8, 0.45, 0.5, 15), (2, 9, 0.85, 1.0, 15), (3, 9, 0.65, 0.375, 30), (7, 9, 0.45, 0.25, 30),
                 (3, 6, 0.85, 0.5, 15), (4, 6, 0.75, 0.75, 15), (2, 8, 0.95, 0.25, 0), (7, 7, 0.75, 0.625, 15)]
BROAD_RES = ["sib", "sig_obs", "tai_in", "tai_ex"]


def broad_points():
    """(set, E, x, Q2, z, pt, phi) of the broad comparison: the grid points
    above the resonance region (W > 2 GeV) with x < 1, then the nodes."""
    import itertools
    pts = []
    g = BROAD_GRID
    for Q2, y, z, pt, phi in itertools.product(g["Q2"], g["y"], g["z"], g["pt"], g["phi"]):
        x = x_of(Q2, y)
        if x < 1 and M * M + Q2 * (1 / x - 1) > 4.0:
            pts.append(("grid", E, x, Q2, z, pt, float(phi)))
    pts += [("node", E, Q2 / (2 * M * nu), Q2, z, pt, float(phi)) for Q2, nu, z, pt, phi in EXTREME_NODES]
    return pts


def load_broad(table, recompute, jobs):
    """Per point of broad_points(): the C++ and FORTRAN sigma_Born, sig_obs,
    tai_in, tai_ex (nan where a code returns nothing: the FORTRAN's model
    stops when its <cos phi> is too large); cached in broad_comparison.tsv."""
    path = os.path.join(FIGDATA, "broad_comparison.tsv")
    cols = ["set", "E", "x", "Q2", "z", "pt", "phi"] + ["%s_%s" % (k, c) for c in ("cpp", "f77") for k in BROAD_RES]
    if os.path.exists(path) and not recompute:
        out = []
        for l in open(path):
            if l.startswith("#"):
                continue
            f = l.split()
            out.append(dict(zip(cols, [f[0]] + [float(v) for v in f[1:]])))
        return out
    if not table:
        raise SystemExit("broad comparison: no cache in %s, and --table not given" % path)

    def res(binary, p):
        r = run_driver(binary, p, table if binary == "rc_point" else None).get("RES", {})
        ok = r.get("kin_ok", 0) and r.get("sib", 0) > 0
        return [r[k] if ok else float("nan") for k in BROAD_RES]

    def one(q):
        p = q[1:]
        return dict(zip(cols, list(q) + res("rc_point", p) + res("haprad2_point", p)))

    with ThreadPoolExecutor(jobs) as pool:
        out = list(pool.map(one, broad_points()))
    with open(path, "w") as f:
        f.write("# C++ (rc_point, table %s) vs HAPRAD 2.0 (haprad2_point), E = %g GeV; nan: no result\n# %s\n"
                % (os.path.basename(table), E, " ".join(cols)))
        for r in out:
            f.write(r["set"] + " " + " ".join("%.10g" % r[k] for k in cols[1:]) + "\n")
    return out


def fig19(broad):
    import numpy as np
    plt = mpl()
    eta = lambda r, c: (r["sig_obs_" + c] + r["tai_in_" + c]) / r["sib_" + c]
    grid = [r for r in broad if r["set"] == "grid"]
    both = [r for r in grid if np.isfinite(r["sib_cpp"]) and np.isfinite(r["sib_f77"])]
    refused = [r for r in grid if np.isfinite(r["sib_cpp"]) and not np.isfinite(r["sib_f77"])]
    born = np.array([abs(r["sib_cpp"] / r["sib_f77"] - 1) for r in both])
    de = np.array([eta(r, "cpp") - eta(r, "f77") for r in both])
    classes = [("table Born within 5%", born < 0.05, CAT[0]), ("5-20%", (born >= 0.05) & (born <= 0.2), CAT[2]),
               ("off by > 20%", born > 0.2, CAT[1])]
    fig, ax = plt.subplot_mosaic([["a", "b"], ["c", "d"]], figsize=(13, 8.6), constrained_layout=True)
    rows = []
    a = ax["a"]
    for label, sel, col in classes:
        ec, ef = np.array([eta(r, "cpp") for r in both])[sel], np.array([eta(r, "f77") for r in both])[sel]
        a.scatter(ef, ec, s=16, color=col, edgecolors="none", alpha=0.85,
                  label="%s: %d points, median |Δη| %.4f" % (label, sel.sum(), np.median(np.abs(de[sel])) if sel.any() else 0))
    lim = (0.6, 1.9)
    a.plot(lim, lim, color=AXIS, linewidth=1)
    a.set_xlim(lim)
    a.set_ylim(lim)
    a.set_xlabel("η, HAPRAD 2.0 FORTRAN")
    a.set_ylabel("η, C++ with the tabulated model")
    a.legend(loc="upper left", title="C++ table Born vs FORTRAN Born:", title_fontsize=9)
    style(a, "η at %d grid points both codes evaluate" % len(both))
    b = ax["b"]
    b.scatter(np.maximum(born, 1e-4), np.maximum(np.abs(de), 1e-5), s=16, color=CAT[0], edgecolors="none", alpha=0.85)
    b.set_xscale("log")
    b.set_yscale("log")
    b.set_xlabel("|σ_Born C++ / σ_Born FORTRAN − 1|: how well the table reproduces the model")
    b.xaxis.label.set_fontsize(10)
    b.set_ylabel("|η C++ − η FORTRAN|")
    b.axhline(np.median(np.abs(de[born < 0.05])), color=AXIS, linewidth=1)
    style(b, "The largest differences: where the table misses the model's Born")
    b.annotate("%d more points: the FORTRAN's model stops\n(|4 <cos φ>| > 0.9), the C++ runs" % len(refused), (0.02, 0.97),
               xycoords="axes fraction", ha="left", va="top", color=INK2, fontsize=9)
    rows += [("grid", r["x"], r["Q2"], r["z"], r["pt"], r["phi"], eta(r, "cpp"), eta(r, "f77"), r["tai_ex_cpp"], r["tai_ex_f77"])
             for r in both]
    c = ax["c"]
    tail = np.array([r["tai_ex_cpp"] / r["tai_ex_f77"] for r in both if r["tai_ex_f77"] > 0])
    order = np.sort(tail)
    c.plot(np.arange(len(order)), order - 1., color=CAT[0], marker="o", markersize=3, linewidth=1)
    c.axhline(0., color=AXIS, linewidth=1)
    c.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    c.set_xlabel("grid points, sorted")
    c.set_ylabel("exclusive tail, C++ / FORTRAN − 1")
    style(c, "Exclusive tail: within %.0e of the FORTRAN at all %d points" % (np.max(np.abs(tail - 1)), len(tail)))
    c.annotate("adaptive cubature with the FORTRAN's D01FCE settings,\nincluding its minimum of 10,000 evaluations per sub-region",
               (0.98, 0.05), xycoords="axes fraction", ha="right", va="bottom", color=INK2, fontsize=9)
    d = ax["d"]
    nodes = [r for r in broad if r["set"] == "node"]
    data = []
    if os.path.exists(os.path.join(FARMTSV, "nodes_0.tsv")):
        for k in (0, 2):
            axes_k, rc, valid, _ = node_grid(k)
            idx = [tuple(int(np.argmin(np.abs(axes_k[i] - v)))
                         for i, v in enumerate((r["Q2"], r["Q2"] / (2 * M * r["x"]), r["z"], r["pt"], r["phi"]))) for r in nodes]
            data.append([rc[i] if valid[i] else np.nan for i in idx])
    xs = np.arange(len(nodes))
    for k, vals, col in zip((0, 2), data, (RAMP3[0], RAMP3[2])):
        v = np.array(vals)
        d.scatter(xs, np.where(v > 0, v, np.nan), s=26, marker="s", color=col, label="data table, iteration %d (≤ 0 not shown)" % k)
    d.scatter(xs, [eta(r, "cpp") for r in nodes], s=30, color=CAT[1], label="model table, C++", zorder=3)
    d.scatter(xs, [eta(r, "f77") for r in nodes], s=60, facecolors="none", edgecolors=INK, linewidths=1.2,
              label="model, FORTRAN (missing: model stops)", zorder=4)
    d.set_yscale("log")
    d.set_xticks(xs)
    d.set_xticklabels(["%g\n%g" % (r["Q2"], round(r["Q2"] / (2 * M * r["x"]))) for r in nodes], fontsize=8)
    d.set_xlabel(r"node: Q$^2$ [GeV$^2$] over $\nu$ [GeV]")
    d.set_ylabel("η at the node")
    d.legend(loc="upper right", fontsize=8)
    style(d, "The full run's extreme nodes, with the model instead")
    rows += [("node", r["x"], r["Q2"], r["z"], r["pt"], r["phi"], eta(r, "cpp"), eta(r, "f77"), r["tai_ex_cpp"], r["tai_ex_f77"])
             for r in nodes]
    g = BROAD_GRID
    save(fig, "fig19_port_broad", "The C++ against the HAPRAD 2.0 FORTRAN across the RG-E kinematics, same physics on both sides",
         "Grid: Q2 in %s GeV^2, y in %s, z in %s, p_t in %s GeV, phi in %s deg, W > 2 GeV (%d points). Both codes with HAPRAD 2.0's "
         "PDF x FF model: the FORTRAN directly, the C++ through the tabulated model (validation/tables, make_model_table). η "
         "without the exclusive tail. The virtual + soft part differs by about 0.4%% by design (exp(δ_inf) in the C++). "
         "Nodes: the 17 RC-grid nodes with the most extreme factors in the full LD2+C run. "
         % (g["Q2"], g["y"], g["z"], g["pt"], g["phi"], len(grid)) + FARM_NOTE,
         ["set", "x", "Q2", "z", "pt", "phi", "eta_cpp", "eta_f77", "tai_ex_cpp", "tai_ex_f77"], rows)


def parse_normalisation(path):
    """{grid: {(comparison, quantity): (median, lo16, hi84)}} from closure_normalisation.py's output."""
    out, grid, comp = {}, None, None
    for line in open(path):
        m = re.match(r"=== (\d+)x\d+x\d+x\d+ grid", line)
        if m:
            grid, comp = int(m.group(1)), None
            out[grid] = {}
            continue
        m = re.match(r"  (producer / model|truth / model|producer / truth) \(", line)
        if m and grid and "median per bin" not in line:
            comp = m.group(1)
            continue
        m = re.match(r"\s+(A|Ac/A|Acc/A)\s+([\d.]+)\s+\[([\d.]+), ([\d.]+)\]", line)
        if m and grid and comp:
            out[grid][(comp, m.group(1))] = tuple(float(m.group(i)) for i in (2, 3, 4))
        if "median per bin" in line:
            comp = None
    return out


def fig06():
    plt = mpl()
    res = parse_normalisation(os.path.join(FIGDATA, "closure_normalisation.txt"))
    grids = sorted(res)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    comps = [("producer / truth", "table / true cell average"), ("truth / model", "true cell average / value at centre"),
             ("producer / model", "table / model (end to end)")]
    rows = []
    for ax, q in zip(axes, ("A", "Ac/A", "Acc/A")):
        for i, ((comp, label), c) in enumerate(zip(comps, CAT)):
            pts = [(g, res[g][(comp, q)]) for g in grids if (comp, q) in res[g]]
            xs = [g + (i - 1) * 0.12 for g, _ in pts]
            med = [v[0] for _, v in pts]
            ax.errorbar(xs, med, yerr=[[v[0] - v[1] for _, v in pts], [v[2] - v[0] for _, v in pts]], color=c, marker="o",
                        capsize=0, elinewidth=1.2, label=label)
            rows += [(q, comp, g, *v) for g, v in pts]
        ax.axhline(1., color=AXIS, linewidth=1)
        ax.set_xticks(grids)
        ax.set_xticklabels(["%d" % g for g in grids])
        ax.set_xlabel("cells per axis of the 4-D table")
        style(ax, q)
    axes[0].set_ylabel("ratio, median over cells [16%, 84%]")
    axes[0].legend(loc="upper right")
    save(fig, "fig06_normalisation", "Closure of the table's absolute normalisation (MakePhiTable)",
         "3M toy pions uniform in (Q2, nu, z, pt2, phi), weighted by HAPRAD 2.0's Born cross section with known luminosity "
         "(closure_normalisation.py). The table matches the true cell average; the rest is bin-centring, which shrinks with cell size.",
         ["quantity", "comparison", "cells_per_axis", "median", "p16", "p84"], rows)


def parse_iteration(path):
    """Toy closure: RC grid vs truth per k (median, p95), and the change between iterations (median, p95)."""
    truth, change = [], []
    for line in open(path):
        m = re.match(r"\s+k=(\d+)\s.*\|\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+cells", line)
        if m:
            truth.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
        m = re.match(r"\s+RC change vs iteration (\d+) over \d+ nodes: 95% (\S+), median (\S+),", line)
        if m:
            change.append((int(m.group(1)), float(m.group(3)), float(m.group(2))))
    return truth, change


FARMTSV = os.path.join(FIG, "farm", "tsv")
FARM_NOTE = ("50 inbending LD2+C runs (dc_fmt pass-1 ntuples), farm workflow rc_inb_C_D2_dc_fmt_20261017_1402, "
             "exported with export_results.sh and export_farm.C.")
RAMP3 = ["#86b6ef", "#3987e5", "#0d366b"]  # ordinal blue steps 250 / 400 / 700: iterations 0, 1, 2
EMPTY = "#e1e0d9"                          # neutral: a cell or region with nothing in it


def farm_tsv(name):
    """A TSV written by export_farm.C, as {column: numpy array}; None if absent."""
    import numpy as np
    path = os.path.join(FARMTSV, name)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    cols = lines[0].rstrip("\n").split("\t")
    data = np.loadtxt(lines[1:], delimiter="\t", ndmin=2)
    return {c: data[:, i] for i, c in enumerate(cols)}


def node_grid(k):
    """Iteration k's RC grid as (axes, rc, valid, nodes): axes the node values
    per dimension (Q2, nu, z, pt, phi), rc = r_vr + r_in and valid as 5-D
    arrays, nodes the raw columns."""
    import numpy as np
    n = farm_tsv("nodes_%d.tsv" % k)
    keys = ["Q2", "nu", "z", "pt", "phi"]
    axes = [np.unique(n[c]) for c in keys]
    idx = tuple(np.searchsorted(a, n[c]) for a, c in zip(axes, keys))
    rc = np.zeros(tuple(len(a) for a in axes))
    valid = np.zeros(rc.shape, bool)
    rc[idx] = n["r_vr"] + n["r_in"]
    valid[idx] = n["valid"] > 0.5
    return axes, rc, valid, n


def grid_interp(axes, rc, valid, v):
    """Multilinear interpolation of rc at v = (Q2, nu, z, pt, phi), as
    RCGrid::Interpolate does: corners without a result are dropped and the
    rest renormalised; nan outside the grid or with no valid corner."""
    import itertools
    import numpy as np
    i0, f = [], []
    for a, x in zip(axes, v):
        t = (x - a[0]) / (a[1] - a[0])
        if t < -1e-9 or t > len(a) - 1 + 1e-9:
            return float("nan")
        i = min(max(int(np.floor(t)), 0), len(a) - 2)
        i0.append(i)
        f.append(min(max(t - i, 0.), 1.))
    s = wsum = 0.
    for corner in itertools.product((0, 1), repeat=5):
        w, ii = 1., []
        for d, u in enumerate(corner):
            w *= f[d] if u else 1. - f[d]
            ii.append(i0[d] + u)
        if w == 0. or not valid[tuple(ii)]:
            continue
        s += w * rc[tuple(ii)]
        wsum += w
    return s / wsum if wsum > 0 else float("nan")


def nearest_rank(v, q):
    """The q-quantile of sorted v as MakeRCGrid compare takes it."""
    return v[min(len(v) - 1, int(q * (len(v) - 1) + 0.5))]


def fig07():
    import numpy as np
    plt = mpl()
    truth, change = parse_iteration(os.path.join(FIGDATA, "closure_iteration.txt"))
    grids = [node_grid(k) for k in range(3)]
    diffs, data = [], []
    for a in range(2):
        both = grids[a][2] & grids[a + 1][2]
        d = np.sort(np.abs(grids[a + 1][1][both] - grids[a][1][both]))
        diffs.append(d)
        data.append((a, nearest_rank(d, 0.5), nearest_rank(d, 0.95)))
    fig, ax = plt.subplot_mosaic([["a", "b", "c"], ["d", "d", "d"]], figsize=(13, 8), constrained_layout=True,
                                 gridspec_kw={"height_ratios": [1, 0.9]})
    a0, a1, a2, a3 = ax["a"], ax["b"], ax["c"], ax["d"]
    a0.plot([t[0] for t in truth], [t[1] for t in truth], color=RAMP2[0], marker="o", label="median node")
    a0.plot([t[0] for t in truth], [t[2] for t in truth], color=RAMP2[1], marker="o", label="95th percentile node")
    a0.set_xlabel("iteration k")
    a0.set_xticks([t[0] for t in truth])
    a0.set_ylabel(r"|$\Delta$ RC factor|")
    a0.legend(loc="upper right")
    style(a0, "Toy: iteration k against the truth")
    for axx, idx, title in ((a1, 1, "Change between iterations: median node"), (a2, 2, "... and 95th percentile node")):
        axx.plot([c[0] + 0.5 for c in change], [c[idx] for c in change], color=CAT[0], marker="o", label="toy (closure)")
        axx.plot([d[0] + 0.5 for d in data], [d[idx] for d in data], color=CAT[1], marker="o", label="RG-E LD2+C data")
        axx.set_xticks([0.5, 1.5, 2.5])
        axx.set_xticklabels(["0→1", "1→2", "2→3"])
        axx.set_xlabel("iterations")
        style(axx, title)
    a1.legend(loc="upper right")
    for axx in (a0, a1, a2):
        axx.set_yscale("log")
        axx.set_ylim(1e-7, 0.5)
    a2.axhline(0.002, color=AXIS, linewidth=1)
    a2.annotate("rc_iterate.py tolerance", (2.6, 0.002), xytext=(0, 3), textcoords="offset points", ha="right", color=MUTED,
                fontsize=9)
    for d, label, c in zip(diffs, ("iteration 0 → 1", "iteration 1 → 2"), RAMP2):
        frac = 1. - np.arange(len(d)) / len(d)
        a3.step(np.maximum(d, 1e-8), frac, where="post", color=c, linewidth=2,
                label="%s: %d nodes, %.1f%% change by > 0.1" % (label, len(d), 100. * np.mean(d > 0.1)))
    a3.set_xscale("log")
    a3.set_yscale("log")
    a3.set_xlim(1e-7, 2e3)
    a3.set_ylim(1e-5, 1.2)
    a3.axvline(0.002, color=AXIS, linewidth=1)
    a3.set_xlabel(r"|$\Delta$ RC factor| between iterations")
    a3.set_ylabel("share of nodes with a larger change")
    a3.legend(loc="lower left")
    style(a3, "Data: every node valid in both iterations")
    rows = [("toy_vs_truth", k, m, p) for k, m, p in truth] + [("toy_change", k, m, p) for k, m, p in change] + \
           [("data_change", k, m, p) for k, m, p in data]
    save(fig, "fig07_convergence", "The RC iteration: the bulk converges like the toy, a tail of nodes does not",
         "Toy: 2M pions weighted by HAPRAD 2.0's Born cross section, observed = Born x RC_true (closure_iteration.py). Data: "
         + FARM_NOTE, ["series", "k", "median", "p95"], rows)


def event_hist(k):
    """Iteration k's selected-event RC factors: (log10 edges, counts, n_le0)."""
    import numpy as np
    path = os.path.join(FARMTSV, "logrc_%d.tsv" % k)
    le0 = float(open(path).readline().split(":")[1].split()[0])
    h = farm_tsv("logrc_%d.tsv" % k)
    return h["log10rc_low"], h["log10rc_high"], h["events"], le0


def share_above(k, t):
    """Share of iteration k's selected events with a factor above t, from the
    histogram (linear within the bin that contains t)."""
    import numpy as np
    lo, hi, n, le0 = event_hist(k)
    L = math.log10(t)
    above = n[lo >= L].sum()
    part = (lo < L) & (hi > L)
    above += (n[part] * (hi[part] - L) / (hi[part] - lo[part])).sum()
    return above / (n.sum() + le0)


def fig08():
    import numpy as np
    plt = mpl()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True, gridspec_kw={"width_ratios": [1.6, 1]})
    a1.axvspan(0.8, 1.2, color=BAND, linewidth=0, label="range of the model's factors (fig02)")
    rows = []
    for k, c in zip(range(3), RAMP3):
        lo, hi, n, le0 = event_hist(k)
        inner = (lo > -100) & (hi < 100)
        lo_i, n_i = lo[inner], n[inner]
        # Rebin 1600 bins of 0.00375 in log10 into 200 of 0.03.
        nb = 8
        m = len(n_i) // nb * nb
        edges = np.append(lo_i[:m:nb], lo_i[m - 1] + (lo_i[1] - lo_i[0]))
        counts = n_i[:m].reshape(-1, nb).sum(axis=1)
        total = n.sum() + le0
        frac = counts / total
        cum = np.cumsum(n_i) / n_i.sum()
        median = 10. ** np.interp(0.5, cum, lo_i + (lo_i[1] - lo_i[0]))
        a1.stairs(np.maximum(frac, 1e-12), 10. ** edges, color=c, linewidth=1.8, label="iteration %d: median %.3f" % (k, median))
        rows += [("iteration_%d" % k, 10. ** a, 10. ** b, f) for a, b, f in zip(edges[:-1], edges[1:], frac)]
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xlim(0.3, 3000)
    a1.set_ylim(1e-9, 0.2)
    a1.set_xlabel("RC factor η (without exclusive tail)")
    a1.set_ylabel("share of selected events per bin")
    a1.legend(loc="upper right")
    style(a1, "RC factor of the 211 M selected events")
    th = [1.5, 2., 10.]
    w = 0.26
    for j, (k, c) in enumerate(zip(range(3), RAMP3)):
        pos = [i + (j - 1) * (w + 0.03) for i in range(len(th))]
        vals = [100. * share_above(k, t) for t in th]
        a2.bar(pos, vals, width=w, color=c, label="iteration %d" % k)
        for p, v in zip(pos, vals):
            a2.annotate("%.1f" % v, (p, v), xytext=(0, 3), textcoords="offset points", ha="center", color=INK2, fontsize=8.5)
        rows += [("above_%g_iteration_%d" % (t, k), t, 0., v) for t, v in zip(th, vals)]
    a2.set_xticks(range(len(th)))
    a2.set_xticklabels(["> %g" % t for t in th])
    a2.set_xlabel("RC factor")
    a2.set_ylabel("share of selected events [%]")
    a2.legend(loc="upper right")
    style(a2, "Events with an implausible factor")
    save(fig, "fig08_data_rc", "RG-E LD2+C: the bulk of the RC factors is sensible, a tail is not, and the iteration grows it",
         "The band marks the range HAPRAD gives with HAPRAD 2.0's own model (fig02), as a scale for plausible factors. Table from data, "
         "not acceptance-corrected; 91% of selected events have full RC-grid support, 0.03% none. " + FARM_NOTE,
         ["series", "rc_low_or_threshold", "rc_high", "share"], rows)


def fig09():
    import numpy as np
    plt = mpl()
    model = [float(l.split("\t")[7]) for l in open(os.path.join(FIGDATA, "shapes_model.tsv")) if not l.startswith(("#", "Q2"))]
    c = farm_tsv("cells_0.tsv")
    fit = c["state"] == 0
    pt = np.sqrt(c["pt2"][fit])
    a, rc, rcc = c["A"][fit], c["Ac_red"][fit] * pt / c["A"][fit], c["Acc_red"][fit] * c["pt2"][fit] / c["A"][fit]
    phis = np.radians(np.arange(0, 361))
    data = np.min(1. + np.outer(rc, np.cos(phis)) + np.outer(rcc, np.cos(2 * phis)), axis=1)
    sets = [("HAPRAD 2.0 model table", np.array(model), CAT[0]), ("RG-E LD2+C data table, 50 runs", data, CAT[1])]
    fig, ax = plt.subplots(figsize=(9.5, 4.6), constrained_layout=True)
    edges = np.linspace(-0.5, 1.2, 69)
    rows = []
    for label, v, col in sets:
        counts = np.histogram(np.clip(v, edges[0] + 1e-9, None), edges)[0] / len(v)
        ax.stairs(counts, edges, color=col, linewidth=2,
                  label="%s: %d cells, %.0f%% at or below 0" % (label, len(v), 100. * np.mean(v <= 0)))
        rows += [(label, a_, b_, k_) for a_, b_, k_ in zip(edges[:-1], edges[1:], counts)]
    ax.axvline(0., color=AXIS, linewidth=1)
    ax.annotate("Born cross section\nnegative at some φ", (-0.03, 0.97), xycoords=("data", "axes fraction"), ha="right",
                va="top", color=INK2, fontsize=9)
    ax.set_xlabel(r"minimum over $\phi_h$ of  1 + (A$_c$/A) cos $\phi_h$ + (A$_{cc}$/A) cos 2$\phi_h$  (at the cell centre)")
    ax.set_ylabel("fraction of fitted cells")
    ax.set_xlim(edges[0], edges[-1])
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.82))
    style(ax)
    save(fig, "fig09_table_shapes", "Why: the data table's phi modulations carry the detector acceptance",
         "Each fitted cell's Born phi shape at its centre. Physical modulations keep it well above 0. Data: the iteration-0 table "
         "of the full run, raw yields (not acceptance-corrected); values below -0.5 are counted in the first bin. " + FARM_NOTE,
         ["table", "bin_low", "bin_high", "fraction"], rows)


def heat_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list("blue_seq", ["#cde2fb", "#6da7ec", "#256abf", "#0d366b"])
    cm.set_bad(EMPTY)
    return cm


def fig11():
    import numpy as np
    plt = mpl()
    projections = [("Q2_y", "y", r"Q$^2$ [GeV$^2$]"), ("z_pt2", r"$p_t^2$ [GeV$^2$]", "z"), ("phi_y", "y", r"$\phi_h$ [deg]")]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6), constrained_layout=True)
    rows, mesh = [], None
    for r, k in enumerate((0, 2)):
        for col, (v, xl, yl) in enumerate(projections):
            m = farm_tsv("map_%s_%d.tsv" % (v, k))
            xe, ye = np.unique(np.append(m["x_low"], m["x_high"])), np.unique(np.append(m["y_low"], m["y_high"]))
            allc = m["all"].reshape(len(xe) - 1, len(ye) - 1).T
            gt2 = m["gt2"].reshape(len(xe) - 1, len(ye) - 1).T
            share = np.ma.masked_where(allc < 1000, gt2 / np.maximum(allc, 1))
            ax = axes[r][col]
            mesh = ax.pcolormesh(xe, ye, share, cmap=heat_cmap(), vmin=0, vmax=1, shading="flat")
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
            style(ax, "Iteration %d" % k if col == 0 else None)
            ax.grid(False)
            tot = 100. * m["gt2"].sum() / m["all"].sum()
            ax.annotate("%.1f%% of these events" % tot, (0.98, 0.97), xycoords="axes fraction", ha="right", va="top", color=INK2,
                        fontsize=9, bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none", alpha=0.85))
            rows += [("%s_iter%d" % (v, k), a, b, c, d, n, g) for a, b, c, d, n, g in
                     zip(m["x_low"], m["x_high"], m["y_low"], m["y_high"], m["all"], m["gt2"])]
    cb = fig.colorbar(mesh, ax=axes, shrink=0.8, pad=0.01)
    cb.set_label("share of the bin's events with RC factor > 2", color=INK2)
    cb.outline.set_visible(False)
    save(fig, "fig11_where_rc_gt2", "Where the implausible factors are: events with RC factor > 2",
         "Selected events with an RC factor, binned in the event's own kinematics; grey: bins with fewer than 1000 events. "
         + FARM_NOTE, ["map", "x_low", "x_high", "y_low", "y_high", "events", "events_rc_gt2"], rows)


def fig12():
    import numpy as np
    plt = mpl()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True, gridspec_kw={"width_ratios": [1.1, 1, 1]})
    a0 = axes[0]
    rows = []
    states = [("fitted from data", CAT[0]), ("filled from fitted neighbours", CAT[2]), ("empty", EMPTY)]
    for k in range(3):
        st = farm_tsv("cells_%d.tsv" % k)["state"]
        left = 0
        for s, (label, col) in enumerate(states):
            n = int(np.sum(st == s))
            a0.barh(k, n, left=left, height=0.55, color=col, label=label if k == 0 else None)
            if s < 2:
                a0.annotate("%d" % n, (left + n / 2., k), ha="center", va="center", color=SURFACE, fontsize=8.5)
            left += n + 400
            rows.append(("cells_iter%d" % k, label, n))
    a0.set_yticks(range(3))
    a0.set_yticklabels(["iteration %d" % k for k in range(3)])
    a0.set_ylim(2.6, -0.9)
    a0.set_xlabel("table cells, thousands (%d in all)" % len(st))
    from matplotlib.ticker import FuncFormatter
    a0.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%d" % (v / 1000)))
    a0.legend(loc="upper right", ncol=1, fontsize=8.5)
    style(a0, "The table: what the tails can draw on")
    mesh = None
    for ax, k in zip(axes[1:], (0, 2)):
        axes_k, rc, valid, _ = node_grid(k)
        q2, nu = axes_k[0], axes_k[1]
        bad = (valid & ((rc < 0.5) | (rc > 2.))).sum(axis=(2, 3, 4))
        nval = valid.sum(axis=(2, 3, 4))
        share = np.ma.masked_where(nval == 0, bad / np.maximum(nval, 1))
        dq, dn = q2[1] - q2[0], nu[1] - nu[0]
        mesh = ax.pcolormesh(np.append(nu - dn / 2, nu[-1] + dn / 2), np.append(q2 - dq / 2, q2[-1] + dq / 2), share,
                             cmap=heat_cmap(), vmin=0, vmax=0.3, shading="flat")
        ax.set_xlabel(r"$\nu$ [GeV]   (y = $\nu$/E)")
        ax.set_ylabel(r"Q$^2$ [GeV$^2$]")
        ax.grid(False)
        style(ax, "RC-grid nodes, iteration %d" % k)
        ax.grid(False)
        ax.annotate("%d of %d valid nodes\noutside [0.5, 2]" % (int(bad.sum()), int(nval.sum())), (0.03, 0.97),
                    xycoords="axes fraction", ha="left", va="top", color=INK2, fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none", alpha=0.85))
        for i, qq in enumerate(q2):
            for j, nn in enumerate(nu):
                rows.append(("nodes_iter%d" % k, "Q2=%g nu=%g" % (qq, nn), "%d/%d" % (bad[i, j], nval[i, j])))
    cb = fig.colorbar(mesh, ax=axes[1:], shrink=0.85, pad=0.01)
    cb.set_label("share of valid nodes outside [0.5, 2]", color=INK2)
    cb.outline.set_visible(False)
    save(fig, "fig12_coverage", "Coverage: a third of the table holds data or a fill, and the extreme grid nodes sit at the edges",
         "Left: 4-D table cells (18 Q2 x 16 nu x 20 z x 24 pt2). Fill: TSemiInclusiveModel's second-neighbour copy (a fitted cell "
         "with A > 0 in the 3^4 block). Maps: valid RC-grid nodes summed over z, p_t and phi; grey: no valid node. " + FARM_NOTE,
         ["series", "category", "value"], rows)


def fig13(scan):
    import numpy as np
    plt = mpl()
    grids = [node_grid(k) for k in range(3)]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True, sharey=True)
    rows = []
    for ax, y in zip(axes, (0.5, 0.8)):
        model = [scan[(BASE["Q2"], y, BASE["z"], BASE["pt"], p)] for p in PHIS]
        ax.plot(PHIS, [(r[1] + r[2]) / r[0] for r in model], color=CAT[1], marker="o", label="HAPRAD 2.0 model table")
        nu = y * E
        for k, c in zip(range(3), RAMP3):
            phis = np.arange(0, 181, 5)
            eta = [grid_interp(grids[k][0], grids[k][1], grids[k][2], (BASE["Q2"], nu, BASE["z"], BASE["pt"], p)) for p in phis]
            ax.plot(phis, eta, color=c, label="data table, iteration %d" % k)
            rows += [("data_iter%d_y%g" % (k, y), float(p), e) for p, e in zip(phis, eta)]
        rows += [("model_y%g" % y, float(p), (r[1] + r[2]) / r[0]) for p, r in zip(PHIS, model)]
        ax.axhline(1., color=AXIS, linewidth=1)
        ax.set_xlim(0, 180)
        ax.set_xticks(range(0, 181, 30))
        ax.set_xlabel(r"$\phi_h$ [deg]")
        style(ax, r"y = %g  ($\nu$ = %.2f GeV, x = %.3f)" % (y, nu, x_of(BASE["Q2"], y)))
    axes[0].set_ylabel(r"$\eta = \sigma_{\mathrm{obs}} / \sigma_{\mathrm{Born}}$")
    axes[0].legend(loc="best")
    save(fig, "fig13_data_vs_model", r"The same kinematics, two tables: Q$^2$ = %g GeV$^2$, z = %g, $p_t$ = %g GeV"
         % (BASE["Q2"], BASE["z"], BASE["pt"]),
         "Model: rc_point with the tabulated HAPRAD 2.0 PDF x FF model. Data: the RC grid of each iteration interpolated as ApplyRC "
         "does. Both without the exclusive tail. " + FARM_NOTE, ["series", "phi_deg", "eta"], rows)


# ==== What the program produces, iteration by iteration ====

ETA_CLIP = 2.0                             # fig14: points above are drawn at the edge, with their value


def eta_cmap():
    """Diverging: blue below 1 (the correction raises the yield), orange above."""
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list("eta_div", ["#0d366b", "#3987e5", "#cde2fb", SURFACE, "#fbd3bd", "#eb6834", "#8a2f0c"])
    cm.set_bad(EMPTY)
    return cm


def fig14():
    import numpy as np
    plt = mpl()
    grids = [node_grid(k) for k in range(3)]
    axes_k = grids[0][0]
    iz, ip = 2, 3                          # z = 0.35, p_t = 0.375 GeV
    q2s, nus = (1, 3, 5), (1, 3, 5, 7)     # Q2 = 2, 4, 6 GeV^2; nu = 3, 5, 7, 9 GeV
    fig, axes = plt.subplots(len(q2s), len(nus), figsize=(13, 8.4), constrained_layout=True, sharex=True, sharey=True)
    rows = []
    phis = axes_k[4]
    for r, iq in enumerate(q2s):
        for c, inu in enumerate(nus):
            ax = axes[r][c]
            Q2, nu = axes_k[0][iq], axes_k[1][inu]
            y = nu / E
            for k, col in enumerate(RAMP3):
                rc, valid = grids[k][1][iq, inu, iz, ip], grids[k][2][iq, inu, iz, ip]
                eta = np.where(valid, rc, np.nan)
                ax.plot(phis, np.minimum(eta, ETA_CLIP), color=col, marker="o", markersize=3.5,
                        linewidth=1.6 if k else 2.4, label="iteration %d" % k)
                rows += [("iter%d" % k, Q2, nu, float(p), e) for p, e in zip(phis, eta)]
                for p, e in zip(phis, eta):
                    if e > ETA_CLIP:
                        ax.annotate("%.1f" % e, (p, ETA_CLIP), xytext=(4, -10 - 9 * k), textcoords="offset points", color=col,
                                    fontsize=8)
            if not np.any(grids[0][2][iq, inu, iz, ip]):
                ax.annotate("no valid node", (0.5, 0.75), xycoords="axes fraction", ha="center", va="center", color=MUTED)
            ax.axhline(1., color=AXIS, linewidth=1)
            style(ax)
            ax.set_title(r"Q$^2$ = %g, y = %.2f, x = %.2f" % (Q2, y, x_of(Q2, y)), color=INK, loc="left", fontsize=10)
            if r == len(q2s) - 1:
                ax.set_xlabel(r"$\phi_h$ [deg]")
            if c == 0:
                ax.set_ylabel(r"$\eta = \sigma_{\mathrm{obs}} / \sigma_{\mathrm{Born}}$")
    ax = axes[0][0]
    ax.set_xlim(-5, 185)
    ax.set_xticks(range(0, 181, 45))
    ax.set_ylim(0.6, ETA_CLIP + 0.05)
    axes[-1][0].legend(*axes[0][0].get_legend_handles_labels(), loc="lower center", fontsize=10)
    save(fig, "fig14_output_phi", r"The program's output: the RC factor against $\phi_h$, z = %g, $p_t$ = %g GeV, three iterations"
         % (axes_k[2][iz], axes_k[3][ip]),
         "Each point is one RC-grid node of the full run: eta = (sigma_F + sigma_in) / sigma_Born from HAPRAD with that iteration's "
         "data table, without the exclusive tail. Events take 1/eta, interpolated between nodes. Values above %g are drawn at the "
         "edge and labelled; missing points: no valid node. Iteration 0 uses the raw-yield table, iterations 1 and 2 the table "
         "refilled with the previous weights. " % ETA_CLIP + FARM_NOTE,
         ["series", "Q2", "nu", "phi_deg", "eta"], rows)


def phi_average(rc, valid):
    """The phi average (Simpson over the 13 phi nodes) of every (Q2, nu, z, pt)
    node line whose phi nodes are all valid; nan otherwise."""
    import numpy as np
    w = np.array([1, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 1], float)
    full = valid.all(axis=4)
    return np.where(full, (rc * w).sum(axis=4) / w.sum(), np.nan)


def fig15():
    import numpy as np
    plt = mpl()
    from matplotlib.colors import TwoSlopeNorm
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), constrained_layout=True, sharey=True)
    rows, med = [], []
    for k in range(3):
        axes_k, rc, valid, _ = node_grid(k)
        avg = phi_average(rc, valid)
        m = np.full(avg.shape[:2], np.nan)
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                v = avg[i, j][np.isfinite(avg[i, j])]
                if len(v):
                    m[i, j] = np.median(v)
        med.append(m)
    q2, nu = axes_k[0], axes_k[1]
    dq, dn = q2[1] - q2[0], nu[1] - nu[0]
    xe, ye = np.append(nu - dn / 2, nu[-1] + dn / 2), np.append(q2 - dq / 2, q2[-1] + dq / 2)
    norm = TwoSlopeNorm(1., 0.7, 1.3)
    for k in range(3):
        ax = axes[k]
        mesh = ax.pcolormesh(xe, ye, np.ma.masked_invalid(med[k]), cmap=eta_cmap(), norm=norm, shading="flat")
        for i in range(len(q2)):
            for j in range(len(nu)):
                if np.isfinite(med[k][i, j]):
                    v = med[k][i, j]
                    ax.annotate("%.2f" % v, (nu[j], q2[i]), ha="center", va="center", fontsize=7.5,
                                color=SURFACE if abs(v - 1.) > 0.18 else INK2)
                    rows.append(("iter%d" % k, q2[i], nu[j], v))
        style(ax, "Iteration %d" % k)
        ax.grid(False)
        ax.set_xlabel(r"$\nu$ [GeV]   (y = $\nu$/E)")
    axes[0].set_ylabel(r"Q$^2$ [GeV$^2$]")
    cb = fig.colorbar(mesh, ax=axes[:3], shrink=0.85, pad=0.01, extend="both")
    cb.set_label(r"median $\langle\eta\rangle_\phi$", color=INK2)
    cb.outline.set_visible(False)
    d = med[2] - med[0]
    ax = axes[3]
    from matplotlib.colors import LinearSegmentedColormap
    dcm = LinearSegmentedColormap.from_list("d", ["#1c5cab", "#86b6ef", SURFACE, "#f3a37b", "#c4480f"])
    dcm.set_bad(EMPTY)
    mesh = ax.pcolormesh(xe, ye, np.ma.masked_invalid(d), cmap=dcm, vmin=-0.1, vmax=0.1, shading="flat")
    for i in range(len(q2)):
        for j in range(len(nu)):
            if np.isfinite(d[i, j]):
                ax.annotate("%+.3f" % d[i, j], (nu[j], q2[i]), ha="center", va="center", fontsize=7, color=INK2)
                rows.append(("change_0_to_2", q2[i], nu[j], d[i, j]))
    style(ax, "Change, iteration 0 → 2")
    ax.grid(False)
    ax.set_xlabel(r"$\nu$ [GeV]   (y = $\nu$/E)")
    cb = fig.colorbar(mesh, ax=ax, shrink=0.85, pad=0.01, extend="both")
    cb.set_label(r"$\Delta$ median $\langle\eta\rangle_\phi$", color=INK2)
    cb.outline.set_visible(False)
    save(fig, "fig15_output_maps", r"The program's output over the kinematics: the $\phi$-averaged RC factor, three iterations",
         "Each cell: the median over z and p_t of the RC-grid nodes' phi average (Simpson over the 13 phi nodes; node lines with "
         "every phi node valid). Below 1 (blue) the correction raises the yield, above 1 (orange) it lowers it. Grey: no complete "
         "node line. Without the exclusive tail. " + FARM_NOTE, ["series", "Q2", "nu", "value"], rows)


def fig16():
    import numpy as np
    plt = mpl()
    nodes = [farm_tsv("nodes_%d.tsv" % k) for k in range(3)]
    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    rows = []
    e0 = np.linspace(0.6, 1.1, 101)
    e1 = np.linspace(-0.5, 2.0, 126)
    for k, (n, col) in enumerate(zip(nodes, RAMP3)):
        v = n["valid"] > 0.5
        h0 = np.histogram(np.clip(n["r_vr"][v], e0[0], e0[-1] - 1e-9), e0)[0] / v.sum()
        h1 = np.histogram(np.clip(n["r_in"][v], e1[0], e1[-1] - 1e-9), e1)[0] / v.sum()
        a0.stairs(h0, e0, color=col, linewidth=2.4 - 0.5 * k, label="iteration %d: median %.4f" % (k, np.median(n["r_vr"][v])))
        a1.stairs(np.maximum(h1, 1e-7), e1, color=col, linewidth=2.4 - 0.5 * k,
                  label="iteration %d: median %.4f" % (k, np.median(n["r_in"][v])))
        rows += [("r_vr_iter%d" % k, a, b, c) for a, b, c in zip(e0[:-1], e0[1:], h0)]
        rows += [("r_in_iter%d" % k, a, b, c) for a, b, c in zip(e1[:-1], e1[1:], h1)]
    a0.set_xlabel(r"$\sigma_F / \sigma_{\mathrm{Born}} = e^{\delta_{\mathrm{inf}}}(1 + \delta_{VR} + \delta_{\mathrm{vac}})$")
    a0.set_ylabel("share of valid nodes per bin")
    a0.legend(loc="upper left")
    style(a0, "Virtual + soft part: kinematics only")
    a1.set_yscale("log")
    a1.set_ylim(1e-6, 1)
    a1.set_xlabel(r"$\sigma_{\mathrm{in}} / \sigma_{\mathrm{Born}}$  (clipped to [-0.5, 2])")
    a1.legend(loc="upper right")
    style(a1, "Inelastic tail: depends on the table")
    both = (nodes[0]["valid"] > 0.5) & (nodes[2]["valid"] > 0.5)
    x, y = nodes[0]["r_in"][both], nodes[2]["r_in"][both]
    lim = (-0.2, 1.0)
    h, xe, ye = np.histogram2d(np.clip(x, *lim), np.clip(y, *lim), bins=120, range=(lim, lim))
    from matplotlib.colors import LogNorm
    cm = heat_cmap()
    cm.set_bad(SURFACE)
    mesh = a2.pcolormesh(xe, ye, np.ma.masked_equal(h.T, 0), cmap=cm, norm=LogNorm(1, h.max()), shading="flat")
    a2.plot(lim, lim, color=AXIS, linewidth=1)
    a2.set_xlim(lim)
    a2.set_ylim(lim)
    a2.set_xlabel(r"$\sigma_{\mathrm{in}} / \sigma_{\mathrm{Born}}$, iteration 0")
    a2.set_ylabel(r"$\sigma_{\mathrm{in}} / \sigma_{\mathrm{Born}}$, iteration 2")
    style(a2, "Node by node, iteration 0 → 2")
    a2.grid(False)
    cb = fig.colorbar(mesh, ax=a2, pad=0.01)
    cb.set_label("nodes", color=INK2)
    cb.outline.set_visible(False)
    save(fig, "fig16_output_parts", "The two parts of the RC factor, iteration by iteration",
         "Every valid RC-grid node of the full run: eta = sigma_F/sigma_Born + sigma_in/sigma_Born (exclusive tail not applied). "
         "sigma_F/sigma_Born depends only on the kinematics, so it barely changes; the inelastic tail integrates the table over "
         "the photon's phase space, so it follows the table. Values outside the axes are counted in the edge bins. " + FARM_NOTE,
         ["series", "bin_low", "bin_high", "fraction"], rows)


def cell_shapes(c):
    """A table's (A, Ac/A, Acc/A) per cell, from the stored reduced amplitudes."""
    a = c["A"]
    with __import__("numpy").errstate(divide="ignore", invalid="ignore"):
        return a, c["Ac_red"] * c["pt2"] ** 0.5 / a, c["Acc_red"] * c["pt2"] / a


def fig17():
    import numpy as np
    plt = mpl()
    cs = [farm_tsv("cells_%d.tsv" % k) for k in range(3)]
    fit = np.all([c["state"] == 0 for c in cs], axis=0)
    sh = [[v[fit] for v in cell_shapes(c)] for c in cs]
    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    rows = []
    e = np.linspace(0.4, 1.4, 101)
    for k, col in ((1, RAMP3[1]), (2, RAMP3[2])):
        r = sh[k][0] / sh[0][0]
        h = np.histogram(np.clip(r, e[0], e[-1] - 1e-9), e)[0] / len(r)
        a0.stairs(h, e, color=col, linewidth=2, label="iteration %d / iteration 0: median %.3f" % (k, np.median(r)))
        rows += [("A%d_over_A0" % k, a, b, c) for a, b, c in zip(e[:-1], e[1:], h)]
    a0.axvline(1., color=AXIS, linewidth=1)
    a0.set_ylim(0, a0.get_ylim()[1] * 1.3)
    a0.set_xlabel("A, refilled with the weights / A from the raw yield")
    a0.set_ylabel("share of cells per bin")
    a0.legend(loc="upper left")
    style(a0, "Corrected yield per cell")
    ed = np.linspace(-0.25, 0.25, 101)
    for ax, j, name in ((a1, 1, r"A$_c$/A"), (a2, 2, r"A$_{cc}$/A")):
        for (k0, k1), col in (((0, 1), RAMP3[1]), ((1, 2), RAMP3[2])):
            d = sh[k1][j] - sh[k0][j]
            h = np.histogram(np.clip(d, ed[0], ed[-1] - 1e-9), ed)[0] / len(d)
            ax.stairs(np.maximum(h, 1e-6), ed, color=col, linewidth=2,
                      label="iteration %d → %d: median %+.4f" % (k0, k1, np.median(d)))
            rows += [("d_%s_%d_%d" % (["", "Ac", "Acc"][j], k0, k1), a, b, c) for a, b, c in zip(ed[:-1], ed[1:], h)]
        ax.axvline(0., color=AXIS, linewidth=1)
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 10.)
        ax.set_xlabel(r"$\Delta$ %s per cell" % name)
        ax.legend(loc="upper left")
        style(ax, r"The cos %s modulation" % (r"$\phi_h$" if j == 1 else r"2$\phi_h$"))
    save(fig, "fig17_output_table", "What the weights do to the table: yields and modulations, iteration by iteration",
         "The %d cells fitted in all three iterations. Iteration k >= 1 refills the table from the same events weighted by "
         "1/eta of iteration k-1, so A_k/A_0 is the cell's mean weight and the changes of Ac/A and Acc/A are the modulations the "
         "correction removes. Values outside the axes are counted in the edge bins. " % fit.sum() + FARM_NOTE,
         ["series", "bin_low", "bin_high", "fraction"], rows)


def fig18():
    import numpy as np
    plt = mpl()
    cs = [farm_tsv("cells_%d.tsv" % k) for k in range(3)]
    fit = np.all([c["state"] == 0 for c in cs], axis=0)
    picks = [(2.0, 3.0), (2.0, 5.0), (2.0, 7.0), (4.0, 5.0), (4.0, 7.0), (4.0, 9.0)]
    z0, pt20 = 0.35, 0.15
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.4), constrained_layout=True, sharex=True)
    phis = np.arange(0, 181, 2)
    rows = []
    c0 = cs[0]
    for ax, (q2, nu) in zip(axes.flat, picks):
        dist = ((c0["Q2"] - q2) / 0.5) ** 2 + ((c0["nu"] - nu) / 0.5) ** 2 + ((c0["z"] - z0) / 0.05) ** 2 + \
               ((c0["pt2"] - pt20) / 0.0625) ** 2
        dist[~fit] = np.inf
        i = int(np.argmin(dist))
        a0 = cs[0]["A"][i]
        for k, col in enumerate(RAMP3):
            a, rc, rcc = (v[i] for v in cell_shapes(cs[k]))
            shape = a * (1. + rc * np.cos(np.radians(phis)) + rcc * np.cos(2. * np.radians(phis))) / a0
            ax.plot(phis, shape, color=col, linewidth=2.4 - 0.5 * k,
                    label="iteration %d: A$_c$/A = %+.3f, A$_{cc}$/A = %+.3f" % (k, rc, rcc))
            rows += [("cell%d_iter%d" % (i, k), float(p), s) for p, s in zip(phis, shape)]
        ax.axhline(0., color=AXIS, linewidth=1)
        y = c0["nu"][i] / E
        style(ax)
        ax.set_title(r"Q$^2$ = %.2f, y = %.2f, z = %.3f, $p_t^2$ = %.3f, %d events" % (c0["Q2"][i], y, c0["z"][i], c0["pt2"][i],
                     c0["events"][i]), color=INK, loc="left", fontsize=9.5)
        ax.legend(loc="lower right", fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel(r"$\phi_h$ [deg]")
        ax.set_xticks(range(0, 181, 45))
    fig.supylabel(r"A (1 + (A$_c$/A) cos $\phi_h$ + (A$_{cc}$/A) cos 2$\phi_h$) / A$_0$  (A$_0$: iteration 0)", color=INK2,
                  fontsize=11)
    save(fig, "fig18_output_cells", r"Single table cells: the fitted $\phi_h$ shape from the raw yield and after each correction",
         "The table cell nearest each kinematic point among those fitted in every iteration, at its centre; normalised to the "
         "iteration-0 A. Iteration 0: raw yields (acceptance included); iterations 1, 2: the same events weighted by 1/eta of "
         "the previous iteration. " + FARM_NOTE, ["series", "phi_deg", "shape"], rows)


# ==== Page sets in the style of V. Klimenko's HAPRAD pi0 study ====
# (CLAS collaboration meeting, May 2025): one page per bin, a grid of small
# panels, the curves of the HAPRAD 2.0 model table's grid and of the full LD2+C
# run's final grid (iteration 2) overlaid, read at the grid nodes; eta without
# the exclusive tail. Triangles: beyond the panel's range, drawn at its edge.
# fig21 divides each table's cross sections by its own phi-averaged sigma_Born,
# since the data table is not absolutely normalised. fig22/23 integrate over
# the 13 phi nodes (Simpson), and use one y range per page, from the 2nd to the
# 98th percentile of its values.

PAGE_TABLES = [("model", "HAPRAD 2.0 model table", CAT[1]), (2, "data table, iteration 2", RAMP3[2])]
SIMPSON13 = [1, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 1]


def node_arrays(k):
    """Grid k (an iteration, or "model") as (axes, eta, sigma_born): 5-D arrays
    over (Q2, nu, z, pt, phi) nodes, nan where the node has no result."""
    import numpy as np
    n = farm_tsv("nodes_%s.tsv" % k)
    keys = ["Q2", "nu", "z", "pt", "phi"]
    axes = [np.unique(n[c]) for c in keys]
    idx = tuple(np.searchsorted(a, n[c]) for a, c in zip(axes, keys))
    shape = tuple(len(a) for a in axes)
    eta, sib = np.full(shape, np.nan), np.full(shape, np.nan)
    ok = n["valid"] > 0.5
    eta[tuple(i[ok] for i in idx)] = (n["r_vr"] + n["r_in"])[ok]
    sib[tuple(i[ok] for i in idx)] = n["sigma_born"][ok]
    return axes, eta, sib


def phi_integrated(eta, sib):
    """The phi-integrated RC factor, int sigma_obs / int sigma_Born (Simpson over
    the 13 phi nodes); nan unless every phi node has a result. The same as the
    ratio of the A terms of A (1 + B cos phi + C cos 2 phi) fits to both."""
    import numpy as np
    w = np.array(SIMPSON13, float)
    return np.sum(eta * sib * w, axis=-1) / np.sum(sib * w, axis=-1)


def page_axes(plt, nrow, ncol):
    plt.rcParams["figure.max_open_warning"] = 0
    # Panels keep one size; small grids get a minimum page, so the title,
    # legend and axis labels have room.
    fig, axes = plt.subplots(nrow, ncol, figsize=(max(1.55 * ncol + 1.6, 10.5), max(1.1 * nrow + 1.6, 4.8)), squeeze=False,
                             sharex=True)
    fig.grid_shape = (nrow, ncol)
    for ax in axes.flat:
        ax.tick_params(labelsize=6, length=2, pad=1)
        style(ax)
        ax.grid(True, color=GRID, linewidth=0.5)
    return fig, axes


def clip_plot(ax, x, v, lo, hi, color, lw, ls="-"):
    """v against x within [lo, hi]; points beyond are drawn on the edge as
    triangles."""
    import numpy as np
    c = np.clip(v, lo, hi)
    ax.plot(x, c, color=color, linewidth=lw, linestyle=ls, marker="o", markersize=1.6)
    for edge, m in ((hi, "^"), (lo, "v")):
        out = np.isfinite(v) & ((v > hi) if m == "^" else (v < lo))
        ax.plot(np.asarray(x)[out], c[out], linestyle="none", marker=m, markersize=3.5, color=color)


def page_limits(curves, model):
    """y range of a panel: what the model spans, widened to the data up to
    0.6 beyond the model on either side."""
    import numpy as np
    m = model[np.isfinite(model)]
    allv = np.concatenate([c[np.isfinite(c)] for c in curves]) if curves else m
    if not len(allv):
        return 0.5, 1.5
    base_lo, base_hi = (m.min(), m.max()) if len(m) else (np.median(allv), np.median(allv))
    lo, hi = max(allv.min(), base_lo - 0.6), min(allv.max(), base_hi + 0.6)
    pad = max(0.02, 0.08 * (hi - lo))
    return lo - pad, hi + pad


def write_pages(name, pages, header, rows):
    """A multi-page PDF, one PNG per page in figures/<name>/, and the TSV."""
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(FIG, name), exist_ok=True)
    with PdfPages(os.path.join(FIG, name + ".pdf")) as pdf:
        for tag, fig in pages:
            pdf.savefig(fig, bbox_inches="tight", pad_inches=0.15)
            fig.savefig(os.path.join(FIG, name, tag + ".png"), dpi=110, bbox_inches="tight", pad_inches=0.15)
            plt.close(fig)
    with open(os.path.join(FIG, name + ".tsv"), "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join("%.6g" % v if isinstance(v, float) else str(v) for v in r) + "\n")
    print("wrote %s (%d pages)" % (name, len(pages)))


def page_legend(fig, entries):
    from matplotlib.lines import Line2D
    fig.legend([Line2D([], [], color=c, linestyle=ls, linewidth=1.6) for _, c, ls in entries], [l for l, _, _ in entries],
               loc="upper right", ncol=len(entries), fontsize=8, frameon=False, bbox_to_anchor=(0.99, 0.995))


def blank_panel(ax):
    """A panel with no node: keep its row/column label, drop ticks, grid and frame."""
    ax.grid(False)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.blank = True


def page_layout(fig, xlabel, ylabel):
    w, h = fig.get_figwidth(), fig.get_figheight()
    nrow, ncol = fig.grid_shape
    left, top = 1.05 / w, 1 - 0.8 / h
    right = min(0.99, left + 1.55 * ncol / w)
    bottom = max(0.75 / h, top - 1.1 * nrow / h)
    for ax in fig.axes:
        if getattr(ax, "blank", False):
            for line in list(ax.lines):
                line.remove()
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom, wspace=0.32, hspace=0.18)
    fig.text(0.12 / w, (top + bottom) / 2, ylabel, rotation=90, ha="left", va="center", fontsize=9, color=INK2)
    fig.text((left + right) / 2, bottom - 0.55 / h, xlabel, ha="center", va="top", fontsize=9, color=INK2)




def fig20_21(which):
    """fig20_rc_phi: eta against phi; fig21_xsec_phi: sigma_obs and sigma_Born
    against phi, each divided by the phi-averaged sigma_Born of its table. One
    page per (Q2, nu), panels z (columns) x p_t (rows)."""
    import numpy as np
    plt = mpl()
    grids = {k: node_arrays(k) for k, _, _ in PAGE_TABLES}
    axes_k = grids["model"][0]
    q2s, nus, zs, pts, phis = axes_k
    tables = PAGE_TABLES
    pages, rows = [], []
    w = np.array(SIMPSON13, float) / sum(SIMPSON13)
    for iq, q2 in enumerate(q2s):
        for inu, nu in enumerate(nus):
            have = np.zeros((len(zs), len(pts)), bool)
            for k, _, _ in tables:
                have |= np.isfinite(grids[k][1][iq, inu]).any(axis=-1)
            if not have.any():
                continue
            zi, pi_ = np.where(have.any(axis=1))[0], np.where(have.any(axis=0))[0][::-1]
            fig, axes = page_axes(plt, len(pi_), len(zi))
            for r, ip in enumerate(pi_):
                for c, iz in enumerate(zi):
                    ax = axes[r][c]
                    curves = []
                    for k, label, col in tables:
                        eta, sib = grids[k][1][iq, inu, iz, ip], grids[k][2][iq, inu, iz, ip]
                        if which == "rc":
                            curves.append((k, col, eta, "-"))
                            rows += [(k, q2, nu, zs[iz], pts[ip], p, e) for p, e in zip(phis, eta)]
                        else:
                            norm = np.sum(sib * w) if np.all(np.isfinite(sib)) else np.nan
                            curves.append((k, col, eta * sib / norm, "-"))
                            curves.append((k, col, sib / norm, "--"))
                            rows += [(k, q2, nu, zs[iz], pts[ip], p, e * s / norm, s / norm) for p, e, s in zip(phis, eta, sib)]
                    model = [v for k, _, v, ls in curves if k == "model"]
                    lo, hi = page_limits([v for _, _, v, _ in curves], np.concatenate(model) if model else np.array([]))
                    for k, col, v, ls in curves:
                        if np.isfinite(v).any():
                            clip_plot(ax, phis, v, lo, hi, col, 1.4 if k == "model" else 1.1, ls)
                    if not any(np.isfinite(v).any() for _, _, v, _ in curves):
                        blank_panel(ax)
                    ax.set_ylim(lo, hi)
                    ax.set_xlim(0, 180)
                    ax.set_xticks([0, 90, 180])
                    if which == "rc":
                        ax.axhline(1., color=AXIS, linewidth=0.6)
                    if c == 0:
                        ax.set_ylabel("%g" % pts[ip], fontsize=8, color=INK2)
                    if r == len(pi_) - 1:
                        ax.set_xlabel("%g" % zs[iz], fontsize=8, color=INK2)
            y = nu / E
            what = (r"RC factor $\eta$ against $\phi_h$" if which == "rc" else
                    r"$\sigma_{\mathrm{obs}}$ (solid) and $\sigma_{\mathrm{Born}}$ (dashed) against $\phi_h$, "
                    r"each / its $\phi$-averaged $\sigma_{\mathrm{Born}}$")
            fig.suptitle(r"%s:  Q$^2$ = %g GeV$^2$, $\nu$ = %g GeV  (y = %.2f, x = %.2f)" % (what, q2, nu, y, x_of(q2, y)),
                         x=0.01, ha="left", fontsize=10, color=INK)
            page_legend(fig, [(l, col, "-") for _, l, col in tables])
            page_layout(fig, r"z  (columns);  each panel: $\phi_h$ = 0-180 deg", r"$p_t$ [GeV]  (rows)")
            pages.append(("Q2_%g_nu_%g" % (q2, nu), fig))
    if which == "rc":
        write_pages("fig20_rc_phi", pages, ["table", "Q2", "nu", "z", "pt", "phi", "eta"], rows)
    else:
        write_pages("fig21_xsec_phi", pages, ["table", "Q2", "nu", "z", "pt", "phi", "sigma_obs_rel", "sigma_born_rel"], rows)


def fig22_23(which):
    """The phi-integrated RC factor, one page per Q2: against p_t^2 with rows z
    (fig22_rc_pt2), or against z with rows p_t (fig23_rc_z); columns nu."""
    import numpy as np
    plt = mpl()
    grids = {k: node_arrays(k) for k, _, _ in PAGE_TABLES}
    q2s, nus, zs, pts, _ = grids["model"][0]
    integ = {k: phi_integrated(g[1], g[2]) for k, g in grids.items()}
    pages, rows = [], []
    for iq, q2 in enumerate(q2s):
        have = np.zeros((len(nus), len(zs), len(pts)), bool)
        for k in integ:
            have |= np.isfinite(integ[k][iq])
        if not have.any():
            continue
        ni = np.where(have.any(axis=(1, 2)))[0]
        if which == "pt2":
            ri, xs, rlab, xlab = np.where(have.any(axis=(0, 2)))[0][::-1], pts ** 2, zs, r"$p_t^2$ [GeV$^2$]"
        else:
            ri, xs, rlab, xlab = np.where(have.any(axis=(0, 1)))[0][::-1], zs, pts, "z"
        fig, axes = page_axes(plt, len(ri), len(ni))
        allv = np.concatenate([integ[k][iq][np.isfinite(integ[k][iq])] for k in integ])
        lo, hi = np.percentile(allv, [2, 98])
        lo, hi = max(lo, 0.5), min(hi, 2.0)
        pad = max(0.02, 0.1 * (hi - lo))
        lo, hi = lo - pad, hi + pad
        for r, i_r in enumerate(ri):
            for c, inu in enumerate(ni):
                ax = axes[r][c]
                drawn = False
                for k, label, col in PAGE_TABLES:
                    v = integ[k][iq, inu, i_r, :] if which == "pt2" else integ[k][iq, inu, :, i_r]
                    if np.isfinite(v).any():
                        clip_plot(ax, xs, v, lo, hi, col, 1.4 if k == "model" else 1.1)
                        drawn = True
                    rows += [(k, q2, nus[inu], rlab[i_r] if which == "pt2" else xx, xx if which == "pt2" else rlab[i_r], vv)
                             for xx, vv in zip(xs, v)]
                ax.axhline(1., color=AXIS, linewidth=0.6)
                ax.set_ylim(lo, hi)
                if not drawn:
                    blank_panel(ax)
                if c == 0:
                    ax.set_ylabel("%g" % rlab[i_r], fontsize=8, color=INK2)
                if r == len(ri) - 1:
                    ax.set_xlabel(r"$\nu$ = %g (y = %.2f)" % (nus[inu], nus[inu] / E), fontsize=7.5, color=INK2)
        fig.suptitle(r"$\phi$-integrated RC factor against %s:  Q$^2$ = %g GeV$^2$" % (xlab, q2), x=0.01, ha="left", fontsize=10,
                     color=INK)
        page_legend(fig, [(l, col, "-") for _, l, col in PAGE_TABLES])
        page_layout(fig, "%s in each panel;  columns: ν [GeV]" % xlab, ("z" if which == "pt2" else r"$p_t$ [GeV]") + "  (rows)")
        pages.append(("Q2_%g" % q2, fig))
    if which == "pt2":
        write_pages("fig22_rc_pt2", pages, ["table", "Q2", "nu", "z", "pt2", "eta_integrated"], rows)
    else:
        write_pages("fig23_rc_z", pages, ["table", "Q2", "nu", "z", "pt", "eta_integrated"], rows)


def fig24():
    """The bulk of the RC grid converges, its tail does not: the change of eta
    = r_vr + r_in between iterations at the nodes valid in all three
    iterations of the full run, by quantile and node by node."""
    import numpy as np
    from matplotlib.colors import LogNorm
    plt = mpl()
    g = [node_grid(k) for k in range(3)]
    both = g[0][2] & g[1][2] & g[2][2]
    eta = [g[k][1][both] for k in range(3)]
    d1, d2 = np.abs(eta[1] - eta[0]), np.abs(eta[2] - eta[1])
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11.5, 4.9), constrained_layout=True, gridspec_kw={"width_ratios": [1, 1.1]})
    rows = []

    # ---- (a) the change between iterations, by quantile of the nodes ----
    quants = [(50, "median"), (90, "90%"), (95, "95%"), (99, "99%"), (99.9, "99.9%"), (100, "largest")]
    ramp = ["#86b6ef", "#5b9be9", "#3987e5", "#1c5cab", "#0d366b", "#06203f"]
    for (q, label), col in zip(quants, ramp):
        v = [np.percentile(d, q) for d in (d1, d2)]
        a0.plot([0, 1], v, color=col, marker="o", markersize=6, linewidth=2)
        a0.annotate("%s  %.2g → %.2g" % (label, v[0], v[1]), (1, v[1]), xytext=(8, 0), textcoords="offset points",
                    va="center", color=INK2, fontsize=9)
        rows.append(("quantile_%s" % label, v[0], v[1]))
    a0.axhline(0.002, color=AXIS, linewidth=1)
    a0.annotate("stopping tolerance 0.002", (2.05, 0.002), xytext=(0, -4), textcoords="offset points", ha="right", va="top",
                color=MUTED, fontsize=8.5)
    a0.set_yscale("log")
    a0.set_xticks([0, 1])
    a0.set_xticklabels(["0 → 1", "1 → 2"])
    a0.set_xlim(-0.25, 2.1)
    a0.set_xlabel("iterations compared")
    a0.set_ylabel(r"|$\Delta\eta$| of the node at that quantile")
    style(a0, "Low quantiles shrink, high ones don't")

    # ---- (b) node by node: second change against first ----
    lo = 1e-6
    x, y = np.maximum(d1, lo), np.maximum(d2, lo)
    edges = np.logspace(-6, 3, 91)
    h, xe, ye = np.histogram2d(x, y, bins=(edges, edges))
    cm = heat_cmap()
    cm.set_bad(SURFACE)
    mesh = a1.pcolormesh(xe, ye, np.ma.masked_equal(h.T, 0), cmap=cm, norm=LogNorm(1, h.max()), shading="flat")
    span = np.array([lo, 1e3])
    a1.plot(span, span, color=CAT[1], linewidth=1.5)
    a1.plot(span, 0.1 * span, color=AXIS, linewidth=1.2, linestyle="--")
    r = (d2 + 1e-7) / (d1 + 1e-7)
    shrink, grow = np.mean(r < 0.3), np.mean(r >= 1)
    from matplotlib.lines import Line2D
    a1.legend([Line2D([], [], color=CAT[1], linewidth=1.5), Line2D([], [], color=AXIS, linewidth=1.2, linestyle="--")],
              ["same change again: not converging", "10× smaller each step"], loc="lower right", fontsize=8.5,
              frameon=True, facecolor=SURFACE, edgecolor="none", framealpha=0.95)
    a1.annotate("%.1f%% of nodes: change shrinks > 3×\n%.1f%% of nodes: change grows" % (100 * shrink, 100 * grow),
                (0.03, 0.97), xycoords="axes fraction", ha="left", va="top", color=INK2, fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc=SURFACE, ec="none", alpha=0.9))
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xlim(lo, 1e3)
    a1.set_ylim(lo, 1e3)
    a1.set_xlabel(r"|$\Delta\eta$|, iteration 0 → 1  (values below 10$^{-6}$ at the edge)")
    a1.set_ylabel(r"|$\Delta\eta$|, iteration 1 → 2")
    style(a1, "Node by node: %d nodes" % len(d1))
    a1.grid(False)
    cb = fig.colorbar(mesh, ax=a1, pad=0.01)
    cb.set_label("nodes", color=INK2)
    cb.outline.set_visible(False)

    save(fig, "fig24_convergence_tail", "Distribution of convergence", "", ["series", "a", "b"], rows)


def fig10(farm):
    plt = mpl()
    t = farm["timing"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    hours = [c / 3600. for c in t["cpu_total_s"]]
    a1.bar(range(len(hours)), hours, width=0.6, color=CAT[0])
    for i, h in enumerate(hours):
        a1.annotate("%.1f h" % h if h >= 0.1 else "< 0.1 h", (i, h), xytext=(0, 3), textcoords="offset points", ha="center",
                    color=INK2, fontsize=9)
    a1.set_xticks(range(len(hours)))
    a1.set_xticklabels(["%s\n%d jobs" % (s, j) for s, j in zip(t["steps"], t["jobs"])])
    a1.set_ylabel("CPU hours, 3 iterations")
    style(a1, "CPU: the RC grid is %.0f%%" % (100. * hours[2] / sum(hours)))
    compute = sum(t["wall_max_s"]) / 3600.
    spans = [s / 3600. for s in t["iteration_span_s"]]
    ys = list(range(len(spans)))
    a2.barh(ys, [compute] * len(spans), height=0.55, color=CAT[0], label="computing (longest job of each step)")
    a2.barh(ys, [s - compute for s in spans], left=[compute + 0.03] * len(spans), height=0.55, color=CAT[1],
            label="waiting in the queue")
    for y, s in zip(ys, spans):
        a2.annotate("%.1f h" % s, (s + 0.03, y), xytext=(4, 0), textcoords="offset points", va="center", color=INK2, fontsize=9)
    a2.set_yticks(ys)
    a2.set_yticklabels(["iteration %d" % i for i in ys])
    a2.set_ylim(len(spans) + 0.4, -0.5)
    a2.set_xlim(0, max(spans) * 1.12)
    a2.set_xlabel("hours")
    a2.legend(loc="lower right")
    style(a2, "Wall time per iteration: mostly queue")
    rows = [(s, j, c / 3600.) for s, j, c in zip(t["steps"], t["jobs"], t["cpu_total_s"])] + \
           [("iteration_%d_span_h" % i, 0, s) for i, s in enumerate(spans)]
    save(fig, "fig10_farm_cost", "The farm run: 456 jobs, %.0f CPU hours, %.1f hours from submission to the last weight"
         % (sum(hours), t["submission_to_end_s"] / 3600.),
         "JLab farm (swif2), 50 inbending LD2+C runs, 3 iterations, 50 grid chunks; timing_report.py. Fill and apply scale with the "
         "number of runs, the grid does not.", ["step", "jobs", "cpu_hours"], rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", help="structure-function table for the C++ (needed unless every scan is cached)")
    ap.add_argument("--only", help="comma-separated figure names, e.g. fig02,fig05")
    ap.add_argument("--recompute", action="store_true", help="re-run HAPRAD instead of using the cached scans")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    args = ap.parse_args()
    os.makedirs(FIGDATA, exist_ok=True)
    table = os.path.realpath(args.table) if args.table else None
    want = set(args.only.split(",")) if args.only else None
    do = lambda name: want is None or name in want

    if any(do(f) for f in ("fig01", "fig02", "fig03", "fig04")):
        scan = load_scan(table, args.recompute, args.jobs)
        for f in (fig01, fig02, fig03, fig04):
            if do(f.__name__):
                f(scan)
    if do("fig05"):
        fig05(load_port(table, args.recompute, args.jobs))
    if do("fig19"):
        fig19(load_broad(table, args.recompute, args.jobs))
    if do("fig06"):
        fig06()
    farm = json.load(open(os.path.join(FIGDATA, "farm_C_D2_20261017.json")))
    if do("fig10"):
        fig10(farm)
    full = [f for f in (fig07, fig08, fig09, fig11, fig12, fig13, fig14, fig15, fig16, fig17, fig18, fig24) if do(f.__name__)]
    if full and not os.path.exists(os.path.join(FARMTSV, "nodes_0.tsv")):
        print("skipping %s: no full-run export in %s (export_results.sh, then export_farm.C)"
              % (", ".join(f.__name__ for f in full), FARMTSV))
        full = []
    for f in full:
        f(load_scan(table, False, args.jobs)) if f is fig13 else f()
    page_sets = [("fig20", fig20_21, "rc"), ("fig21", fig20_21, "xsec"), ("fig22", fig22_23, "pt2"), ("fig23", fig22_23, "z")]
    page_sets = [p for p in page_sets if do(p[0])]
    if page_sets and not os.path.exists(os.path.join(FARMTSV, "nodes_model.tsv")):
        print("skipping %s: no %s (run_rc_grid.sh with the model table, then export_farm.C's export_grid)"
              % (", ".join(p[0] for p in page_sets), os.path.join(FARMTSV, "nodes_model.tsv")))
        page_sets = []
    for _, f, which in page_sets:
        f(which)


if __name__ == "__main__":
    main()
