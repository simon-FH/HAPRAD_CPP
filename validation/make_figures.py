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

Figures:
  fig01_rc_vs_phi        eta against phi_h at four y, and its parts
  fig02_rc_kinematics    phi-averaged eta and its range over phi, against y, z, p_t
  fig03_rc_parts         the parts of eta - 1 against y, z, p_t
  fig04_fake_modulations the change of the cos(phi), cos(2 phi) amplitudes
  fig05_port_validation  C++ against HAPRAD 2.0 FORTRAN, point by point
  fig06_normalisation    the absolute-normalisation closure of the table
  fig07_convergence      the RC iteration: toy closure against real data
  fig08_data_rc          RC factors of the selected LD2+C events
  fig09_table_shapes     why: phi shapes of the data table's cells
  fig10_farm_cost        CPU and wall time of the farm run
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


def fig07(farm):
    plt = mpl()
    truth, change = parse_iteration(os.path.join(FIGDATA, "closure_iteration.txt"))
    data = [(0, farm["compare"]["iter0_to_1"][2], farm["compare"]["iter0_to_1"][3]),
            (1, farm["compare"]["iter1_to_2"][2], farm["compare"]["iter1_to_2"][3])]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True, sharey=True)
    a0, a1, a2 = axes
    a0.plot([t[0] for t in truth], [t[1] for t in truth], color=RAMP2[0], marker="o", label="median node")
    a0.plot([t[0] for t in truth], [t[2] for t in truth], color=RAMP2[1], marker="o", label="95th percentile node")
    a0.set_xlabel("iteration k")
    a0.set_xticks([t[0] for t in truth])
    a0.set_ylabel(r"|$\Delta$ RC factor|")
    a0.legend(loc="upper right")
    style(a0, "Toy: iteration k against the truth")
    for ax, idx, title in ((a1, 1, "Change between iterations: median node"), (a2, 2, "... and 95th percentile node")):
        ax.plot([c[0] + 0.5 for c in change], [c[idx] for c in change], color=CAT[0], marker="o", label="toy (closure)")
        ax.plot([d[0] + 0.5 for d in data], [d[idx] for d in data], color=CAT[1], marker="o", label="RG-E LD2+C data")
        ax.set_xticks([0.5, 1.5, 2.5])
        ax.set_xticklabels(["0→1", "1→2", "2→3"])
        ax.set_xlabel("iterations")
        style(ax, title)
    a1.legend(loc="upper right")
    for ax in axes:
        ax.set_yscale("log")
    a2.axhline(0.002, color=AXIS, linewidth=1)
    a2.annotate("rc_iterate.py tolerance", (2.6, 0.002), xytext=(0, 3), textcoords="offset points", ha="right", color=MUTED,
                fontsize=9)
    rows = [("toy_vs_truth", k, m, p) for k, m, p in truth] + [("toy_change", k, m, p) for k, m, p in change] + \
           [("data_change", k, m, p) for k, m, p in data]
    save(fig, "fig07_convergence", "The RC iteration: the bulk converges like the toy, the extreme nodes do not",
         "Toy: 2M pions weighted by HAPRAD 2.0's Born cross section, observed = Born x RC_true (closure_iteration.py). "
         "Data: 50 inbending LD2+C runs, dc_fmt, farm workflow rc_inb_C_D2_dc_fmt_20261017_1402.",
         ["series", "k", "median", "p95"], rows)


def fig08(farm):
    plt = mpl()
    ev = farm["events"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    pcs = ev["percentiles"]
    xs = list(range(len(pcs)))
    a1.axhspan(0.8, 1.2, color=BAND, linewidth=0, label="plausible range")
    for it, c in zip(("iter0", "iter2"), RAMP2):
        a1.plot(xs, ev[it], color=c, marker="o", label="iteration %s" % it[-1])
        a1.annotate("%.3g" % ev[it][-1], (xs[-1], ev[it][-1]), xytext=(6, 0), textcoords="offset points", va="center",
                    color=INK2, fontsize=9, annotation_clip=False)
    a1.set_yscale("log")
    from matplotlib.ticker import FixedLocator, NullLocator
    a1.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20, 50]))
    a1.yaxis.set_minor_locator(NullLocator())
    a1.set_yticklabels(["0.5", "1", "2", "5", "10", "20", "50"])
    a1.set_xticks(xs)
    a1.set_xticklabels(["%g%%" % p for p in pcs])
    a1.set_xlabel("percentile of the selected events")
    a1.set_ylabel("RC factor (without exclusive tail)")
    a1.legend(loc="upper left")
    style(a1, "RC factor of the %.0f M selected events" % (ev["selected"] / 1e6))
    th = ev["above_pct"]["thresholds"]
    w = 0.36
    for j, (it, c) in enumerate(zip(("iter0", "iter2"), RAMP2)):
        pos = [i + (j - 0.5) * (w + 0.04) for i in range(len(th))]
        a2.bar(pos, ev["above_pct"][it], width=w, color=c, label="iteration %s" % it[-1])
        for p, v in zip(pos, ev["above_pct"][it]):
            a2.annotate("%.1f%%" % v, (p, v), xytext=(0, 3), textcoords="offset points", ha="center", color=INK2, fontsize=9)
    a2.set_xticks(range(len(th)))
    a2.set_xticklabels(["> %g" % t for t in th])
    a2.set_xlabel("RC factor")
    a2.set_ylabel("share of selected events [%]")
    a2.legend(loc="upper right")
    style(a2, "Events with an implausible factor")
    rows = [("percentile_%g" % p, ev["iter0"][i], ev["iter2"][i]) for i, p in enumerate(pcs)] + \
           [("above_%g_pct" % t, ev["above_pct"]["iter0"][i], ev["above_pct"]["iter2"][i]) for i, t in enumerate(th)]
    save(fig, "fig08_data_rc", "RG-E LD2+C, 50 inbending runs: the bulk is sensible, a tail of events is not",
         "weight_summary.C over the RC friends of farm workflow rc_inb_C_D2_dc_fmt_20261017_1402. Table from data, not acceptance-"
         "corrected; 91% of selected events have full RC-grid support, 0.03% none.", ["quantity", "iteration_0", "iteration_2"], rows)


def fig09():
    plt = mpl()
    load = lambda f: [float(l.split("\t")[7]) for l in open(os.path.join(FIGDATA, f)) if not l.startswith(("#", "Q2"))]
    sets = [("HAPRAD 2.0 model table", load("shapes_model.tsv"), CAT[0]), ("data table, run 020026", load("shapes_data_020026.tsv"), CAT[1])]
    fig, ax = plt.subplots(figsize=(9, 4.4), constrained_layout=True)
    edges = [-1.0 + 0.05 * i for i in range(45)]
    rows = []
    for label, v, c in sets:
        n = len(v)
        counts = [sum(1 for x in v if a <= x < b) / n for a, b in zip(edges[:-1], edges[1:])]
        ax.stairs(counts, edges, color=c, linewidth=2, label="%s: %d cells, %.0f%% at or below 0" % (label, n, 100. * sum(x <= 0 for x in v) / n))
        rows += [(label, a, b, k) for a, b, k in zip(edges[:-1], edges[1:], counts)]
    ax.axvline(0., color=AXIS, linewidth=1)
    ax.annotate("Born cross section\nnegative at some φ", (-0.02, ax.get_ylim()[1] * 0.92), ha="right", va="top", color=INK2,
                fontsize=9)
    ax.set_xlabel(r"minimum over $\phi_h$ of  1 + (A$_c$/A) cos $\phi_h$ + (A$_{cc}$/A) cos 2$\phi_h$")
    ax.set_ylabel("fraction of fitted cells")
    ax.set_xlim(-1, 1.2)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.8))
    style(ax)
    save(fig, "fig09_table_shapes", "Why: the data table's phi modulations carry the detector acceptance",
         "Each fitted cell's Born phi shape at its centre (table_shapes.C). Physical modulations keep it well above 0. Data: run 020026 "
         "dc_fmt, raw yields (not acceptance-corrected); 464 of the 492 cells below 0.5 have modulations > 3 sigma from zero.",
         ["table", "bin_low", "bin_high", "fraction"], rows)


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
    if do("fig06"):
        fig06()
    farm = json.load(open(os.path.join(FIGDATA, "farm_C_D2_20261017.json")))
    for f in (fig07, fig08):
        if do(f.__name__):
            f(farm)
    if do("fig09"):
        fig09()
    if do("fig10"):
        fig10(farm)


if __name__ == "__main__":
    main()
