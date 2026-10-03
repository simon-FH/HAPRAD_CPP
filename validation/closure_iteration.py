#!/usr/bin/env python3
"""Closure test of the RC iteration (PLAN.md Phase 4).

Toy pions uniform in (Q2, nu, z, pt2, phi), weighted by HAPRAD 2.0's Born cross
section, give a TRUE Born table T_B (MakePhiTable) and its RC grid RC_true. The
same events weighted by sigma_Born x RC_true are the "observed" data. Started
from them, Utilities/rc_iterate.py must find its way back:

    table_k -> T_B   and   grid_k -> RC_true.

T_B is an exact fixed point of the loop (the producer is linear in the
weights, and observed / RC_true = Born event by event), and the events are the
same throughout, so there is no statistical noise in the comparison: what is
left measures the convergence itself. The exclusive tail is off on both sides.

    ./closure_iteration.py [n_events] [--iterations 4] [--jobs N] [--keep DIR]

About 12 minutes for 2M events and 4 iterations on 14 cores; the RC grids
(5,250 nodes each) dominate.
"""

import argparse
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")
UTIL = os.path.join(HERE, "..", "Utilities")
E = "10.5473"
BOX = "1:9,2:9,0.15:0.95,0:1.5"   # Q2, nu, z, pt2 -- inside the RC grid's nodes

CONFIG = """beam_energy = %s
pid = 211
Q2_min = 0
W2_min = 0
y_max = 1
z_min = 0
z_max = 1
pt2_max = 10
Mx_min = 0
weight_branch = w
bins_Q2  = 8 1 9
bins_nu  = 7 2 9
bins_z   = 8 0.15 0.95
bins_pt2 = 6 0 1.5
bins_phi = 12 -180 180
fit_min_phi_bins = 6
fit_min_events = 100
rc_nodes_Q2  = 5 1 9
rc_nodes_nu  = 5 2 9
rc_nodes_z   = 5 0.15 0.95
rc_nodes_pt  = 6 0 1.25
rc_nodes_phi = 7 0 180
target_naz = 0.5
rc_exclusive_tail = no
""" % E

# Copies DT, replacing the weight w: Born events keep w where the RC grid
# covers them (rc_status <= 1), observed ones get w * rc there; 0 elsewhere.
REWEIGHT = r"""
void reweight(const char* in, const char* friendPath, const char* out, int observed) {
  TFile fi(in); TNtuple* dt = (TNtuple*)fi.Get("DT");
  TFile ff(friendPath); TTree* rc = (TTree*)ff.Get("RC");
  Float_t r; Int_t st; rc->SetBranchAddress("rc", &r); rc->SetBranchAddress("rc_status", &st);
  TString vars; for (auto* b : *dt->GetListOfBranches()) vars += TString(vars.Length() ? ":" : "") + b->GetName();
  const int nv = dt->GetListOfBranches()->GetEntries(), iw = dt->GetListOfBranches()->IndexOf(dt->GetBranch("w"));
  TFile fo(out, "RECREATE"); TNtuple nt("DT", dt->GetTitle(), vars);
  std::vector<Float_t> row(nv); Long64_t kept = 0;
  for (Long64_t i = 0; i < dt->GetEntries(); ++i) {
    dt->GetEntry(i); rc->GetEntry(i); Float_t* v = dt->GetArgs();
    for (int k = 0; k < nv; ++k) row[k] = v[k];
    const bool keep = st <= 1 && v[iw] > 0;
    row[iw] = keep ? (observed ? v[iw] * r : v[iw]) : 0.f;
    kept += keep; nt.Fill(row.data());
  }
  fo.cd(); nt.Write(); fo.Close();
  printf("@reweight %lld of %lld events kept\n", kept, dt->GetEntries());
}
"""

# Table k against the true table, over cells fitted in both: A as a ratio, the
# modulations Ac/A and Acc/A as differences (a ratio blows up where Ac ~ 0).
CMPTABLE = r"""
void cmptable(const char* a, const char* b) {
  TFile fa(a), fb(b);
  THnD *A1 = (THnD*)fa.Get("A"), *A2 = (THnD*)fb.Get("A"), *C1 = (THnD*)fa.Get("Ac"), *C2 = (THnD*)fb.Get("Ac");
  THnD *D1 = (THnD*)fa.Get("Acc"), *D2 = (THnD*)fb.Get("Acc"), *F1 = (THnD*)fa.Get("fitted"), *F2 = (THnD*)fb.Get("fitted");
  std::vector<double> r[3]; long oneOnly = 0;
  for (Long64_t i = 0; i < F1->GetNbins(); ++i) {
    const bool f1 = F1->GetBinContent(i) > 0, f2 = F2->GetBinContent(i) > 0;
    if (f1 != f2) ++oneOnly;
    if (!f1 || !f2) continue;
    const double a1 = A1->GetBinContent(i), a2 = A2->GetBinContent(i);
    r[0].push_back(a1 / a2 - 1);
    r[1].push_back(C1->GetBinContent(i) / a1 - C2->GetBinContent(i) / a2);
    r[2].push_back(D1->GetBinContent(i) / a1 - D2->GetBinContent(i) / a2);
  }
  long nTrue = 0; for (Long64_t i = 0; i < F2->GetNbins(); ++i) nTrue += F2->GetBinContent(i) > 0;
  printf("@table %zu %ld %ld", r[0].size(), oneOnly, nTrue);
  for (int k = 0; k < 3; ++k) {
    std::vector<double> s = r[k]; for (auto& x : s) x = fabs(x); std::sort(s.begin(), s.end());
    double m = 0; for (auto x : r[k]) m += x; m /= r[k].size();
    printf(" %.3e %.3e %+.3e", s[s.size() / 2], s[size_t(0.95 * (s.size() - 1))], m);
  }
  printf("\n");
}
"""


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode:
        sys.exit("failed: %s\n%s" % (" ".join(cmd), (r.stdout + r.stderr)[-3000:]))
    return r.stdout


def at(out, tag):
    return [l.split()[1:] for l in out.splitlines() if l.startswith(tag)][0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("n", nargs="?", default="2000000")
    ap.add_argument("--iterations", type=int, default=4)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--keep", help="work in this directory and keep it")
    args = ap.parse_args()

    tmpdir = None
    if args.keep:
        w = os.path.realpath(args.keep)
        os.makedirs(w, exist_ok=True)
    else:
        tmpdir = tempfile.TemporaryDirectory(prefix="iterclosure_")  # removed when main() returns
        w = tmpdir.name
    p = lambda *x: os.path.join(w, *x)
    open(p("closure.cfg"), "w").write(CONFIG)
    open(p("reweight.C"), "w").write(REWEIGHT)
    open(p("cmptable.C"), "w").write(CMPTABLE)
    cfg, jobs = p("closure.cfg"), str(args.jobs)
    mpt, rcg, apply = os.path.join(UTIL, "bin", "MakePhiTable"), os.path.join(UTIL, "run_rc_grid.sh"), os.path.join(UTIL, "bin", "ApplyRC")
    root = lambda macro: run(["root", "-l", "-b", "-q", macro], cwd=w)

    def table(events, name):
        run([mpt, "fill", cfg, p(name + "_counts.root"), events])
        run([mpt, "fit", cfg, p(name + "_counts.root"), p(name + ".root")])
        return p(name + ".root")

    print("events: %s toy pions weighted by HAPRAD 2.0's Born cross section" % args.n, flush=True)
    with open(p("points.txt"), "w") as f:
        subprocess.run([os.path.join(BIN, "make_weighted_events"), "points", "11", args.n, E, BOX], stdout=f, check=True)
    with open(p("points.txt")) as fi, open(p("harm.txt"), "w") as fo:
        subprocess.run([os.path.join(BIN, "born_harmonics")], stdin=fi, stdout=fo, cwd=os.path.join(HERE, "data"), check=True)
    os.remove(p("points.txt"))
    run([os.path.join(BIN, "make_weighted_events"), "build", "11", args.n, E, BOX, p("harm.txt"), p("events_all.root")])

    # The universe: events the RC grid covers. Defined with a first grid, so
    # that the true table is built from exactly the events the loop can weight.
    print("truth: table and RC grid", flush=True)
    run([rcg, cfg, table(p("events_all.root"), "table_pre"), p("grid_pre.root"), jobs])
    run([apply, cfg, p("grid_pre.root"), p("events_all.root"), p("events_all_rc.root")])
    print("  " + " ".join(at(root('reweight.C("%s","%s","%s",0)' % (p("events_all.root"), p("events_all_rc.root"), p("events_born.root"))), "@reweight")[:3]), flush=True)
    truth = table(p("events_born.root"), "table_true")
    run([rcg, cfg, truth, p("grid_true.root"), jobs])
    run([apply, cfg, p("grid_true.root"), p("events_born.root"), p("events_born_rc.root")])
    kept = at(root('reweight.C("%s","%s","%s",1)' % (p("events_born.root"), p("events_born_rc.root"), p("events_obs.root"))), "@reweight")
    print("  observed events: %s of %s (the rest are outside the true grid's support)" % (kept[0], kept[2]), flush=True)

    print("loop: rc_iterate.py on the observed events, %d iterations" % args.iterations, flush=True)
    r = subprocess.run([os.path.join(UTIL, "rc_iterate.py"), cfg, p("loop"), p("events_obs.root"),
                        "--max-iter", str(args.iterations), "--tol", "0", "--jobs", jobs], capture_output=True, text=True)
    if not os.path.exists(p("loop", "iter_%d" % (args.iterations - 1), "grid.root")):
        sys.exit("rc_iterate failed:\n" + r.stdout + r.stderr)

    print("\nagainst the truth (cells / nodes valid in both; the true table has %s fitted cells):"
          % at(root('cmptable.C("%s","%s")' % (truth, truth)), "@table")[2])
    print("  table: A_k/A_true - 1 as |.| median, |.| 95%, mean; Ac/A and Acc/A differences as |.| median, |.| 95%")
    print("  RC grid: RC_k - RC_true as |.| median, |.| 95%, |.| max, mean")
    for k in range(args.iterations):
        t = at(root('cmptable.C("%s","%s")' % (p("loop", "iter_%d" % k, "table.root"), truth)), "@table")
        g = at(run([os.path.join(UTIL, "bin", "MakeRCGrid"), "compare", p("grid_true.root"), p("loop", "iter_%d" % k, "grid.root")]), "@compare")
        print("  k=%d   A %s %s %s   Ac/A %s %s   Acc/A %s %s   | %s %s %s %s   cells %s (+%s in one only), nodes %s"
              % (k, t[3], t[4], t[5], t[6], t[7], t[9], t[10], g[2], g[3], g[1], g[4], t[0], t[1], g[0]))
    print("\nrc_iterate.py's own convergence measure (change between iterations):")
    print("\n".join("  " + l.strip() for l in r.stdout.splitlines() if "RC change" in l))


if __name__ == "__main__":
    main()
