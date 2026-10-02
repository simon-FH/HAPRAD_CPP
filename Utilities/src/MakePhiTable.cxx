// MakePhiTable -- RG-E ntuples -> structure-function table for HAPRAD
// (PLAN.md Phase 1). Successor to PhiHist/phihist.cpp, deleted in b303487.
//
// Two stages, so the expensive one runs on the cluster:
//
//   MakePhiTable fill <config> <counts.root> <ntuple files or globs ...>
//       Applies the cuts and fills 5-D counts in (Q2, nu, z, pt2, phi).
//       Each job can take any subset of the files; outputs merge with `hadd`.
//
//   MakePhiTable fit <config> <counts.root> <table.root>
//       Fits A + Ac cos(phi) + Acc cos(2 phi) in every 4-D cell and writes the
//       table in the format TSemiInclusiveModel reads.
//
// WHAT A IS. HAPRAD turns the table into structure functions through
// TStructFunctionArray, and sigma_Born = K * (A + Ac cos phi + Acc cos 2 phi)
// with K a kinematic factor of the code. K does not match any standard
// definition of a yield (PLAN.md Phase 1), so the table *defines* A by it:
//
//     A(cell) = [measured cross section in the cell] / K(cell centre)
//
// The measured cross section, in HAPRAD's variables (x, y, z, pt2, phi), is the
// fitted yield times the Jacobian from the binning variables,
// |d(Q2, nu)/d(x, y)| = 2 M E nu, over the cell volume (phi in radians, as in
// HAPRAD's sigma) and the luminosity. Defined this way, whatever K's formula is
// cancels in the calculation: tested in PLAN.md 2.3 by replacing it and seeing
// no change. It requires the conversion to use the real beam energy (PLAN.md
// 2.3, Finding 2). Checked end to end by validation/closure_normalisation.py.
//
// The overall scale does NOT cancel in the RC factor: sigma_Born and the
// inelastic tail scale with the table, the exclusive tail (MAID) does not. With
// `luminosity` unset (1) the table is in arbitrary units and the exclusive
// tail's share of the correction is wrong by that factor.
//
// The table is written with pt_scaling = "reduced" (Ac/pt, Acc/pt^2) and
// interpolation = "log", per PLAN.md 0.2 and 2.3.
//
// ACCEPTANCE is not applied yet; the table records `acceptance_corrected = no`
// so that an uncorrected table cannot be mistaken for a real one.

#include "TRadCor.h"
#include "RGESelection.hxx"
#include "haprad_constants.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "THn.h"
#include "TMath.h"
#include "TNamed.h"
#include "TString.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

using rge::Config;
using rge::LoadConfig;
const char* const* kAxis = rge::kTableAxis;

int Fill(const Config& c, const char* outPath, int nIn, char** in) {
  TChain chain("DT");
  TString sources;
  for (int i = 0; i < nIn; ++i) {
    chain.Add(in[i]);
    sources += TString(in[i]) + "\n";
  }

  rge::Event e;
  Float_t w = 1.f;
  if (!rge::BindEvent(chain, e)) return 1;
  // Optional per-event weight -- e.g. 1/acceptance, or a cross section in a
  // closure test. Any overall scale: the fit uses the cell's mean weight.
  if (!c.weightBranch.empty()) {
    chain.SetBranchStatus(c.weightBranch.c_str(), 1);
    if (chain.SetBranchAddress(c.weightBranch.c_str(), &w) < 0) {
      fprintf(stderr, "input lacks weight branch %s\n", c.weightBranch.c_str());
      return 1;
    }
  }

  THnD counts("counts", "counts", 5, c.n, c.lo, c.hi);
  for (Int_t d = 0; d < 5; ++d) counts.GetAxis(d)->SetName(kAxis[d]);
  counts.Sumw2();

  // Every stage that removes events, in order: the shared selection, then the grid.
  const Int_t nStages = rge::kNStages + 1;
  TH1D flow("cutflow", "events surviving each cut", nStages, 0, nStages);
  for (Int_t i = 0; i < rge::kNStages; ++i) flow.GetXaxis()->SetBinLabel(i + 1, rge::kStage[i]);
  flow.GetXaxis()->SetBinLabel(nStages, "inside grid");

  const Double_t radToDeg = 180. / TMath::Pi();
  const Long64_t n = chain.GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    chain.GetEntry(i);
    const Int_t passed = rge::Stages(c, e);
    for (Int_t s = 0; s < passed; ++s) flow.Fill(s);
    if (passed < rge::kNStages) continue;
    const Double_t v[5] = {e.Q2, e.nu, e.zh, e.pt2, e.phi * radToDeg};  // RG-E stores phi_PQ in radians
    Bool_t inside = true;
    for (Int_t d = 0; d < 5; ++d)
      if (v[d] < c.lo[d] || v[d] >= c.hi[d]) inside = false;
    if (!inside) continue;
    flow.Fill(rge::kNStages);
    counts.Fill(v, w);
  }

  TFile out(outPath, "RECREATE");
  if (out.IsZombie()) return 1;
  counts.Write();
  flow.Write();
  TNamed("config", c.text.c_str()).Write();
  TNamed("source_files", sources.Data()).Write();
  out.Close();

  printf("%s\n", outPath);
  for (Int_t i = 0; i < nStages; ++i) printf("  %-14s %12.0f\n", flow.GetXaxis()->GetBinLabel(i + 1), flow.GetBinContent(i + 1));
  return 0;
}

// Weighted linear least squares for y = p0 + p1 cos(phi) + p2 cos(2 phi), with
// each basis function AVERAGED OVER ITS PHI BIN, not taken at the bin centre:
//
//     <cos(n phi)>_bin = cos(n phi_c) * sin(n w/2) / (n w/2),   w = bin width
//
// A counted bin integrates the distribution, so the centre value attenuates the
// harmonics -- by 1.1% for cos(phi) and 4.5% for cos(2 phi) with 12 bins. The
// closure test saw exactly that, as opposite-sign pull offsets.
// Returns false if the normal matrix is singular.
bool FitHarmonics(const std::vector<Double_t>& phi, Double_t width, const std::vector<Double_t>& y,
                  const std::vector<Double_t>& var, Double_t p[3], Double_t cov[3][3], Double_t& chi2) {
  const Double_t s1 = std::sin(width / 2.) / (width / 2.);
  const Double_t s2 = std::sin(width) / width;
  Double_t a[3][3] = {{0}}, b[3] = {0};
  for (size_t k = 0; k < phi.size(); ++k) {
    const Double_t f[3] = {1., s1 * std::cos(phi[k]), s2 * std::cos(2. * phi[k])};
    for (Int_t i = 0; i < 3; ++i) {
      b[i] += f[i] * y[k] / var[k];
      for (Int_t j = 0; j < 3; ++j) a[i][j] += f[i] * f[j] / var[k];
    }
  }
  // Invert the 3x3 symmetric matrix by cofactors.
  const Double_t det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                       a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
  if (!(std::fabs(det) > 1e-300)) return false;
  cov[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
  cov[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
  cov[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
  cov[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
  cov[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
  cov[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
  cov[1][0] = cov[0][1];
  cov[2][0] = cov[0][2];
  cov[2][1] = cov[1][2];
  for (Int_t i = 0; i < 3; ++i) p[i] = cov[i][0] * b[0] + cov[i][1] * b[1] + cov[i][2] * b[2];
  chi2 = 0.;
  for (size_t k = 0; k < phi.size(); ++k) {
    const Double_t r = y[k] - (p[0] + p[1] * s1 * std::cos(phi[k]) + p[2] * s2 * std::cos(2. * phi[k]));
    chi2 += r * r / var[k];
  }
  return true;
}

// K = sigma_Born / A at a point, from TRadCor with a table that is A = 1
// everywhere. Written to a temporary file because that is what the reader
// takes.
class KFactor {
 public:
  explicit KFactor(const std::string& tmpPath) : fPath(tmpPath) {
    Int_t n[4] = {1, 1, 1, 1};
    Double_t lo[4] = {0, 0, 0, 0}, hi[4] = {1e4, 1e4, 1.0001, 1e4};
    const char* ax[4] = {"Q2", "nu", "z", "pt2"};
    THnD a("A", "A", 4, n, lo, hi), ac("Ac", "Ac", 4, n, lo, hi), acc("Acc", "Acc", 4, n, lo, hi);
    for (THnD* h : {&a, &ac, &acc})
      for (Int_t d = 0; d < 4; ++d) h->GetAxis(d)->SetName(ax[d]);
    Int_t one[4] = {1, 1, 1, 1};
    a.SetBinContent(one, 1.);
    {
      TFile f(fPath.c_str(), "RECREATE");
      a.Write();
      ac.Write();
      acc.Write();
      TNamed("pt_scaling", "raw").Write();
    }
    fOk = fRC.LoadSemiInclusiveTable(fPath.c_str());
  }
  ~KFactor() { std::remove(fPath.c_str()); }
  bool Ok() const { return fOk; }

  Double_t operator()(Double_t E, Double_t Q2, Double_t nu, Double_t z, Double_t pt2) {
    std::ostringstream sink;  // THadronKinematics prints a line per evaluation
    std::streambuf* saved = std::cout.rdbuf(sink.rdbuf());
    const Double_t k = fRC.CalculateBorn(E, Q2 / (2. * kMassProton * nu), Q2, z, std::sqrt(pt2), 90.);
    std::cout.rdbuf(saved);
    return k;
  }

 private:
  std::string fPath;
  TRadCor fRC;
  bool fOk;
};

int Fit(const Config& c, const char* countsPath, const char* outPath) {
  std::unique_ptr<TFile> in(TFile::Open(countsPath, "READ"));
  if (!in || in->IsZombie()) {
    fprintf(stderr, "cannot open %s\n", countsPath);
    return 1;
  }
  std::unique_ptr<THnD> counts(in->Get<THnD>("counts"));
  if (!counts || counts->GetNdimensions() != 5) {
    fprintf(stderr, "%s has no 5-D 'counts'\n", countsPath);
    return 1;
  }
  for (Int_t d = 0; d < 5; ++d)
    if (counts->GetAxis(d)->GetNbins() != c.n[d] || std::fabs(counts->GetAxis(d)->GetXmin() - c.lo[d]) > 1e-12 ||
        std::fabs(counts->GetAxis(d)->GetXmax() - c.hi[d]) > 1e-12) {
      fprintf(stderr, "binning of %s does not match the config (axis %s)\n", countsPath, kAxis[d]);
      return 1;
    }

  Int_t n4[4];
  Double_t lo4[4], hi4[4];
  for (Int_t d = 0; d < 4; ++d) {
    n4[d] = c.n[d];
    lo4[d] = c.lo[d];
    hi4[d] = c.hi[d];
  }
  auto make = [&](const char* name) {
    std::unique_ptr<THnD> h(new THnD(name, name, 4, n4, lo4, hi4));
    for (Int_t d = 0; d < 4; ++d) h->GetAxis(d)->SetName(kAxis[d]);
    return h;
  };
  auto A = make("A"), Ac = make("Ac"), Acc = make("Acc"), fitted = make("fitted");
  auto eA = make("A_err"), eAc = make("Ac_err"), eAcc = make("Acc_err");
  auto rawA = make("yield_a"), rawAc = make("yield_ac"), rawAcc = make("yield_acc");
  auto nev = make("events"), chi2ndf = make("chi2_ndf");

  KFactor K(std::string(outPath) + ".unitK.root");
  if (!K.Ok()) return 1;

  const Int_t nphi = c.n[4];
  std::vector<Double_t> phiC(nphi);
  for (Int_t k = 0; k < nphi; ++k) phiC[k] = counts->GetAxis(4)->GetBinCenter(k + 1) / 180. * TMath::Pi();
  const Double_t phiWidth = (c.hi[4] - c.lo[4]) / nphi / 180. * TMath::Pi();
  const Double_t s1 = std::sin(phiWidth / 2.) / (phiWidth / 2.), s2 = std::sin(phiWidth) / phiWidth;

  // The binning volume of a 5-D cell, phi in radians, and the Jacobian
  // |d(Q2, nu)/d(x, y)| = 2 M E nu that turns a yield density in the binning
  // variables into one in HAPRAD's variables.
  Double_t vol = phiWidth;
  for (Int_t d = 0; d < 4; ++d) vol *= (c.hi[d] - c.lo[d]) / c.n[d];

  Long64_t nCells = 0, nEmpty = 0, nSparse = 0, nNegative = 0, nNoK = 0, nBuilt = 0;
  Int_t idx5[5], idx4[4];
  for (idx4[0] = 1; idx4[0] <= n4[0]; ++idx4[0])
    for (idx4[1] = 1; idx4[1] <= n4[1]; ++idx4[1])
      for (idx4[2] = 1; idx4[2] <= n4[2]; ++idx4[2])
        for (idx4[3] = 1; idx4[3] <= n4[3]; ++idx4[3]) {
          ++nCells;
          std::vector<Double_t> y(nphi), var(nphi);
          Double_t total = 0., total2 = 0.;
          Int_t filled = 0;
          for (Int_t d = 0; d < 4; ++d) idx5[d] = idx4[d];
          for (Int_t k = 0; k < nphi; ++k) {
            idx5[4] = k + 1;
            const Long64_t b = counts->GetBin(idx5);
            y[k] = counts->GetBinContent(b);
            total += y[k];
            total2 += counts->GetBinError2(b);
            if (y[k] > 0.) ++filled;
          }
          if (total == 0.) {
            ++nEmpty;
            continue;
          }
          nev->SetBinContent(idx4, total);
          if (filled < c.fitMinPhiBins || total < c.fitMinEvents) {
            ++nSparse;
            continue;
          }
          // Variances from the fitted prediction, not from each bin's own
          // content: with data-derived variances, bins that fluctuate low get
          // small errors and pull the fit down, by about 1/n for n events per
          // phi bin (validation/closure_normalisation.py). r = sum w^2 / sum w is
          // the cell's mean event weight (1 for unweighted events), so a bin
          // expecting mu has variance r * mu; the floor of one event keeps
          // empty and near-empty bins from dominating. Start flat, iterate.
          const Double_t r = total2 / total;
          for (Int_t k = 0; k < nphi; ++k) var[k] = r * std::max(total / nphi, r);
          Double_t p[3], cov[3][3], chi2;
          bool ok = true;
          for (Int_t it = 0; ok && it < 4; ++it) {
            if (it > 0)
              for (Int_t k = 0; k < nphi; ++k) {
                const Double_t mu = p[0] + p[1] * s1 * std::cos(phiC[k]) + p[2] * s2 * std::cos(2. * phiC[k]);
                var[k] = r * std::max(mu, r);
              }
            ok = FitHarmonics(phiC, phiWidth, y, var, p, cov, chi2);
          }
          if (!ok) {
            ++nSparse;
            continue;
          }
          if (!(p[0] > 0.)) {  // log interpolation needs A > 0
            ++nNegative;
            continue;
          }
          const Double_t Q2 = A->GetAxis(0)->GetBinCenter(idx4[0]);
          const Double_t nu = A->GetAxis(1)->GetBinCenter(idx4[1]);
          const Double_t z = A->GetAxis(2)->GetBinCenter(idx4[2]);
          const Double_t pt2 = A->GetAxis(3)->GetBinCenter(idx4[3]);
          const Double_t k = K(c.E, Q2, nu, z, pt2);
          if (!(k > 0.) || !std::isfinite(k)) {
            ++nNoK;
            continue;
          }
          // yield -> cross section in (x, y, z, pt2, phi); then divide by K.
          const Double_t s = 2. * kMassProton * c.E * nu / vol / c.lumi / k;
          const Double_t pt = std::sqrt(pt2);
          A->SetBinContent(idx4, p[0] * s);
          Ac->SetBinContent(idx4, p[1] * s / pt);
          Acc->SetBinContent(idx4, p[2] * s / pt2);
          eA->SetBinContent(idx4, std::sqrt(cov[0][0]) * s);
          eAc->SetBinContent(idx4, std::sqrt(cov[1][1]) * s / pt);
          eAcc->SetBinContent(idx4, std::sqrt(cov[2][2]) * s / pt2);
          rawA->SetBinContent(idx4, p[0]);
          rawAc->SetBinContent(idx4, p[1]);
          rawAcc->SetBinContent(idx4, p[2]);
          chi2ndf->SetBinContent(idx4, nphi > 3 ? chi2 / (nphi - 3) : 0.);
          fitted->SetBinContent(idx4, 1.);
          ++nBuilt;
        }

  TFile out(outPath, "RECREATE");
  if (out.IsZombie()) return 1;
  for (THnD* h : {A.get(), Ac.get(), Acc.get(), fitted.get(), eA.get(), eAc.get(), eAcc.get(), rawA.get(), rawAc.get(),
                  rawAcc.get(), nev.get(), chi2ndf.get()})
    h->Write();
  TNamed("pt_scaling", "reduced").Write();
  TNamed("interpolation", "log").Write();
  TNamed("acceptance_corrected", "no").Write();
  TNamed("normalisation",
         "A = (fitted yield harmonic) x 2 M E nu / (5-D cell volume, phi in rad) / luminosity / "
         "K(cell centre), K = sigma_Born / A of TRadCor at the beam energy. Valid only if "
         "TStructFunctionArray converts with the real beam energy (PLAN.md 2.3, Finding 2).")
      .Write();
  TNamed("luminosity", c.lumi == 1. ? "1 (arbitrary units: the exclusive tail is NOT on the same scale)"
                                    : Form("%.10g events/nb", c.lumi))
      .Write();
  TNamed("beam_energy", Form("%.6g", c.E)).Write();
  TNamed("config", c.text.c_str()).Write();
  if (TNamed* src = in->Get<TNamed>("source_files")) src->Write();
  if (TH1D* flow = in->Get<TH1D>("cutflow")) flow->Write();
  out.Close();

  printf("%s\n", outPath);
  printf("  4-D cells                     %8lld\n", nCells);
  printf("    no events                   %8lld\n", nEmpty);
  printf("    too sparse to fit           %8lld   (< %d filled phi bins or < %g events)\n", nSparse, c.fitMinPhiBins,
         c.fitMinEvents);
  printf("    fitted A <= 0               %8lld\n", nNegative);
  printf("    HAPRAD rejects the centre   %8lld\n", nNoK);
  printf("    written                     %8lld\n", nBuilt);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc >= 5 && std::string(argv[1]) == "fill") {
    Config c;
    if (!LoadConfig(argv[2], c)) return 1;
    return Fill(c, argv[3], argc - 4, argv + 4);
  }
  if (argc == 5 && std::string(argv[1]) == "fit") {
    Config c;
    if (!LoadConfig(argv[2], c)) return 1;
    return Fit(c, argv[3], argv[4]);
  }
  fprintf(stderr,
          "usage: MakePhiTable fill <config> <counts.root> <ntuple files or globs ...>\n"
          "       MakePhiTable fit  <config> <counts.root> <table.root>\n"
          "Counts from several fill jobs merge with `hadd`.\n");
  return 2;
}
