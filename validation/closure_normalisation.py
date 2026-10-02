#!/usr/bin/env python3
"""Normalisation closure for MakePhiTable (PLAN.md Phase 1).

Events uniform in the binning variables, weighted by HAPRAD 2.0's Born cross
section times |d(x,y)/d(Q2,nu)|, are a cross-section measurement with a known
answer and a known luminosity. MakePhiTable turns them into a table;
make_model_table builds the table directly from the same model. If A's
definition is right, producer / model = 1 in every cell -- absolutely, not just
up to a constant. That matters: the exclusive tail is added in absolute units,
so the table's scale is part of the RC factor. A wrong Jacobian would also show
as a trend -- a missing factor of nu, for instance, as a factor ~3 across nu.

The modulations are compared as Ac/A and Acc/A, so that they test the fit and
the p_t scaling independently of the normalisation.

Bin-averaging of a steep cross section leaves deviations that shrink as the
cells get smaller; a normalisation error would not. Hence two grids.

    ./closure_normalisation.py [n_events]
"""

import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")
DATA = os.path.join(HERE, "data")
PRODUCER = os.path.join(HERE, "..", "Utilities", "bin", "MakePhiTable")
E = "10.5473"
BOX = (1.5, 6.0, 4.0, 9.5, 0.2, 0.8, 0.0, 0.8)   # Q2, nu, z, pt2 ranges

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
luminosity = %%.10g
bins_Q2  = %%d %g %g
bins_nu  = %%d %g %g
bins_z   = %%d %g %g
bins_pt2 = %%d %g %g
bins_phi = 12 -180 180
fit_min_phi_bins = 12
fit_min_events = 200
""" % ((E,) + BOX)

MACRO = r"""
// Per cell, for A, Ac/A and Acc/A:
//   producer / model    end to end (expected 1 up to bin-centring)
//   truth / model       bin-centring: the cell average vs the value at the centre
//   producer / truth    MakePhiTable alone (expected 1)
// "truth" is the cell average of the events' own harmonics, times nu_c/nu (the
// Jacobian the producer evaluates at the centre), over the draws expected in the
// cell -- so boundary cells are treated as the producer treats them.
THnD* Empty(THnD* h, const char* n) { THnD* c = (THnD*)h->Clone(n); c->Reset(); return c; }
double Num(std::string s) { for (char& c : s) if (c == 'D' || c == 'd') c = 'e'; return atof(s.c_str()); }

void cmp(const char* prod, const char* model, const char* events, const char* mharm) {
  TFile fp(prod), fm(model), fe(events);
  const char* nm[3] = {"A", "Ac", "Acc"};
  THnD *p[3], *m[3];
  for (int i = 0; i < 3; i++) { p[i] = (THnD*)fp.Get(nm[i]); m[i] = (THnD*)fm.Get(nm[i]); }
  THnD* pf = (THnD*)fp.Get("fitted"); THnD* mf = (THnD*)fm.Get("fitted"); THnD* pe = (THnD*)fp.Get("A_err");
  THnD* A = p[0];

  THnD* cen[3] = {Empty(A, "c0"), Empty(A, "c1"), Empty(A, "c2")};   // model harmonics at the centres
  std::ifstream in(mharm); std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line); Int_t j[4]; int ok; std::string s0, s1, s2;
    if (!(ls >> j[0] >> j[1] >> j[2] >> j[3] >> ok >> s0 >> s1 >> s2) || ok != 1) continue;
    cen[0]->SetBinContent(j, Num(s0)); cen[1]->SetBinContent(j, Num(s1)); cen[2]->SetBinContent(j, Num(s2));
  }

  THnD* sum[3] = {Empty(A, "s0"), Empty(A, "s1"), Empty(A, "s2")};   // true cell sums
  TNtuple* dt = (TNtuple*)fe.Get("DT");
  const double M = 0.938272, mpi = 0.1395675;
  Int_t idx[4];
  for (Long64_t i = 0; i < dt->GetEntries(); ++i) {
    dt->GetEntry(i); Float_t* v = dt->GetArgs();
    const double Q2 = v[2], nu = v[3], W2 = v[6], z = v[7], pt2 = v[8], th = v[10];
    const double Eh = z * nu, ph = sqrt(Eh * Eh - mpi * mpi);
    if (!(W2 + mpi * mpi - 2 * (M + nu) * Eh + 2 * sqrt(nu * nu + Q2) * ph * cos(th) > 0)) continue;  // the producer's Mx cut
    const double x4[4] = {Q2, nu, z, pt2}; bool inside = true;
    for (int d = 0; d < 4; d++) { idx[d] = A->GetAxis(d)->FindFixBin(x4[d]); if (idx[d] < 1 || idx[d] > A->GetAxis(d)->GetNbins()) inside = false; }
    if (!inside) continue;
    const double f = A->GetAxis(1)->GetBinCenter(idx[1]) / nu;
    for (int k = 0; k < 3; k++) sum[k]->AddBinContent(idx, v[13 + k] * f);
  }
  const double draws = ((TParameter<double>*)fe.Get("draws"))->GetVal();
  const double box = ((TParameter<double>*)fe.Get("box_volume"))->GetVal();
  double vol = 2 * TMath::Pi();
  for (int d = 0; d < 4; d++) vol *= A->GetAxis(d)->GetBinWidth(1);
  const double nexp = draws * vol / box;

  // r[comparison][quantity]
  std::vector<double> r[3][3]; std::vector<int> ix[4];
  for (Long64_t b = 0; b < pf->GetNbins(); ++b) {
    if (pf->GetBinContent(b) == 0 || mf->GetBinContent(b) == 0) continue;
    if (pe->GetBinContent(b) > 0.05 * p[0]->GetBinContent(b)) continue;   // statistically poor
    pf->GetBinContent(b, idx);
    const double c0 = cen[0]->GetBinContent(idx), s0 = sum[0]->GetBinContent(idx);
    if (c0 == 0 || s0 == 0) continue;
    const double pa = p[0]->GetBinContent(b), ma = m[0]->GetBinContent(b);
    const double pm[3] = {pa / ma, (p[1]->GetBinContent(b) / pa) / (m[1]->GetBinContent(b) / ma),
                          (p[2]->GetBinContent(b) / pa) / (m[2]->GetBinContent(b) / ma)};
    const double tm[3] = {s0 / nexp / c0, (sum[1]->GetBinContent(idx) / s0) / (cen[1]->GetBinContent(idx) / c0),
                          (sum[2]->GetBinContent(idx) / s0) / (cen[2]->GetBinContent(idx) / c0)};
    for (int k = 0; k < 3; k++) { r[0][k].push_back(pm[k]); r[1][k].push_back(tm[k]); r[2][k].push_back(pm[k] / tm[k]); }
    for (int d = 0; d < 4; d++) ix[d].push_back(idx[d]);
  }
  auto quant = [](std::vector<double> v, double q) {
    std::sort(v.begin(), v.end()); return v.empty() ? 0. : v[size_t(q * (v.size() - 1))]; };
  const char* cmpn[3] = {"producer / model (end to end)", "truth / model (bin-centring)", "producer / truth (MakePhiTable alone)"};
  const char* qn[3] = {"A    ", "Ac/A ", "Acc/A"};
  printf("cells compared %zu.  Median over cells, [16%%, 84%%] range:\n", r[0][0].size());
  for (int c = 0; c < 3; c++) {
    printf("  %s\n", cmpn[c]);
    for (int k = 0; k < 3; k++) printf("    %s  %.4f   [%.3f, %.3f]\n", qn[k], quant(r[c][k], .5), quant(r[c][k], .16), quant(r[c][k], .84));
  }
  const char* var[4] = {"Q2", "nu", "z", "pt2"};
  for (int c = 0; c < 3; c += 2) {
    printf("  %s, median per bin:\n", cmpn[c]);
    for (int k = 0; k < 3; k++) for (int d = 0; d < 4; d++) {
      printf("    %s %-4s", qn[k], var[d]);
      for (int b = 1; b <= A->GetAxis(d)->GetNbins(); b++) {
        std::vector<double> sel;
        for (size_t j = 0; j < r[c][k].size(); j++) if (ix[d][j] == b) sel.push_back(r[c][k][j]);
        if (sel.empty()) printf("     --"); else printf(" %6.3f", quant(sel, .5));
      }
      printf("\n");
    }
  }
}
"""


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def main():
    n = sys.argv[1] if len(sys.argv) > 1 else "3000000"
    tmp = tempfile.mkdtemp(prefix="normclosure_")
    box = "%g:%g,%g:%g,%g:%g,%g:%g" % BOX
    pts = run([os.path.join(BIN, "make_weighted_events"), "points", "7", n, E, box]).stdout
    harm = run([os.path.join(BIN, "born_harmonics")], input=pts, cwd=DATA).stdout
    hpath, events = os.path.join(tmp, "harm.txt"), os.path.join(tmp, "events.root")
    open(hpath, "w").write(harm)
    built = run([os.path.join(BIN, "make_weighted_events"), "build", "7", n, E, box, hpath, events]).stdout
    print(built.strip())
    lumi = float(built.split("luminosity")[1])
    open(os.path.join(tmp, "cmp.C"), "w").write(MACRO)

    for nb in (4, 6, 8):
        cfg = os.path.join(tmp, "n%d.cfg" % nb)
        open(cfg, "w").write(CONFIG % (lumi, nb, nb, nb, nb))
        counts, table = os.path.join(tmp, "c%d.root" % nb), os.path.join(tmp, "t%d.root" % nb)
        run([PRODUCER, "fill", cfg, counts, events])
        run([PRODUCER, "fit", cfg, counts, table])
        spec = ",".join("%d:%g:%g" % (nb, BOX[2 * i], BOX[2 * i + 1]) for i in range(4))
        cells = run([os.path.join(BIN, "make_model_table"), "points", spec, E]).stdout
        mh = run([os.path.join(BIN, "born_harmonics")], input=cells, cwd=DATA).stdout
        mhp, model = os.path.join(tmp, "mh%d.txt" % nb), os.path.join(tmp, "m%d.root" % nb)
        open(mhp, "w").write(mh)
        run([os.path.join(BIN, "make_model_table"), "build", spec, E, mhp, model], cwd=DATA)
        print("\n=== %dx%dx%dx%d grid (x 12 phi) ===" % (nb, nb, nb, nb))
        out = run(["root", "-l", "-b", "-q", 'cmp.C("%s","%s","%s","%s")' % (table, model, events, mhp)], cwd=tmp).stdout
        print("\n".join(l for l in out.splitlines() if not l.startswith("Processing")))


if __name__ == "__main__":
    main()
