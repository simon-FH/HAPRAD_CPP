#!/usr/bin/env python3
"""Closure test for MakePhiTable (PLAN.md Phase 2.2).

Toy pions with a known azimuthal modulation go through MakePhiTable fill + fit,
and the fitted modulation in each cell is compared with the truth -- the mean
of the true rc(pt) and rcc(pt) over that cell's events. If filling, binning
and fitting are right, the pulls (fit - truth) / error have mean 0 and width 1.

    ./closure_producer.py [n_events] [seed]

Needs bin/make_toy_events and ../Utilities/bin/MakePhiTable.
"""

import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PRODUCER = os.path.join(HERE, "..", "Utilities", "bin", "MakePhiTable")
TOY = os.path.join(HERE, "bin", "make_toy_events")

# A coarse grid so the cells are well populated with a modest toy sample.
CONFIG = """beam_energy = 10.5473
pid = 211
Q2_min = 1.0
W2_min = 4.0
y_max = 0.95
z_min = 0.15
z_max = 1.0
pt2_max = 1.5
Mx_min = 1.5
bins_Q2  = 4 1 6
bins_nu  = 4 3 9.5
bins_z   = 4 0.15 0.95
bins_pt2 = 4 0 1.2
bins_phi = 12 -180 180
fit_min_phi_bins = 4
fit_min_events = 200
"""

MACRO = r"""
void closure(const char* toy, const char* table) {
  TFile ft(toy), fb(table);
  TNtuple* dt = (TNtuple*)ft.Get("DT");
  THnD* A   = (THnD*)fb.Get("A");     THnD* Ac  = (THnD*)fb.Get("Ac");   THnD* Acc  = (THnD*)fb.Get("Acc");
  THnD* eAc = (THnD*)fb.Get("Ac_err"); THnD* eAcc = (THnD*)fb.Get("Acc_err");
  THnD* fit = (THnD*)fb.Get("fitted");
  // truth: per cell, mean of rc(pt) and rcc(pt) over the events that passed the cuts
  Int_t n[4]; Double_t lo[4], hi[4];
  for (int d=0; d<4; d++) { n[d]=A->GetAxis(d)->GetNbins(); lo[d]=A->GetAxis(d)->GetXmin(); hi[d]=A->GetAxis(d)->GetXmax(); }
  THnD tRc("tRc","",4,n,lo,hi), tRcc("tRcc","",4,n,lo,hi), tN("tN","",4,n,lo,hi);
  Float_t *v = 0;
  const double M = 0.938272, m = 0.1395675;
  for (Long64_t i = 0; i < dt->GetEntries(); ++i) {
    dt->GetEntry(i); v = dt->GetArgs();
    const double Q2=v[2], nu=v[3], W2=v[6], z=v[7], pt2=v[8], th=v[10], y=v[5];
    if (!(Q2 > 1.0 && W2 > 4.0 && y < 0.95 && z > 0.15 && z < 1.0 && pt2 < 1.5)) continue;
    const double Eh=z*nu, ph=sqrt(Eh*Eh-m*m);
    const double mx2 = W2 + m*m - 2*(M+nu)*Eh + 2*sqrt(nu*nu+Q2)*ph*cos(th);
    if (!(mx2 > 1.5*1.5)) continue;
    double x4[4] = {Q2, nu, z, pt2}; bool in = true;
    for (int d=0; d<4; d++) if (x4[d] < lo[d] || x4[d] >= hi[d]) in = false;
    if (!in) continue;
    const double pt = sqrt(pt2);
    tRc.Fill(x4, -0.15*pt/0.5); tRcc.Fill(x4, 0.05*(pt/0.5)*(pt/0.5)); tN.Fill(x4);
  }
  int nc = 0; double s1=0, s2=0, t1=0, t2=0;
  Int_t idx[4];
  for (Long64_t b = 0; b < A->GetNbins(); ++b) {
    if (fit->GetBinContent(b) == 0) continue;
    A->GetBinContent(b, idx);
    const double N = tN.GetBinContent(tN.GetBin(idx));
    if (N == 0) continue;
    const double pt = sqrt(A->GetAxis(3)->GetBinCenter(idx[3])), pt2 = pt*pt;
    const double a = A->GetBinContent(b);
    // measured ratios; the reduced storage and the normalisation cancel
    const double rc = Ac->GetBinContent(b)*pt/a,  erc  = eAc->GetBinContent(b)*pt/a;
    const double rcc = Acc->GetBinContent(b)*pt2/a, ercc = eAcc->GetBinContent(b)*pt2/a;
    const double trc = tRc.GetBinContent(tRc.GetBin(idx))/N, trcc = tRcc.GetBinContent(tRcc.GetBin(idx))/N;
    const double p1 = (rc-trc)/erc, p2 = (rcc-trcc)/ercc;
    s1 += p1; s2 += p2; t1 += p1*p1; t2 += p2*p2; ++nc;
  }
  const double m1 = s1/nc, m2 = s2/nc;
  printf("cells compared: %d\n", nc);
  printf("pull Ac/A   mean %+.3f  rms %.3f\n", m1, sqrt(t1/nc - m1*m1));
  printf("pull Acc/A  mean %+.3f  rms %.3f\n", m2, sqrt(t2/nc - m2*m2));
}
"""


def main():
    nev = sys.argv[1] if len(sys.argv) > 1 else "2000000"
    seed = sys.argv[2] if len(sys.argv) > 2 else "12345"
    tmpdir = tempfile.TemporaryDirectory(prefix="closure_")  # removed when main() returns
    tmp = tmpdir.name
    cfg, toy = os.path.join(tmp, "toy.cfg"), os.path.join(tmp, "toy.root")
    counts, table = os.path.join(tmp, "counts.root"), os.path.join(tmp, "table.root")
    with open(cfg, "w") as f:
        f.write(CONFIG)
    subprocess.run([TOY, toy, nev, seed], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([PRODUCER, "fill", cfg, counts, toy], check=True)
    subprocess.run([PRODUCER, "fit", cfg, counts, table], check=True)
    with open(os.path.join(tmp, "closure.C"), "w") as f:
        f.write(MACRO)
    subprocess.run(["root", "-l", "-b", "-q", "closure.C(\"%s\",\"%s\")" % (toy, table)], cwd=tmp, check=True)


if __name__ == "__main__":
    main()
