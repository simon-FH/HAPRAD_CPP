// make_model_table -- build a structure-function table from HAPRAD 2.0's own
// PDF x FF model, for PLAN.md Phase 2.3.
//
// Tier 3 compares the C++ against the FORTRAN, but the two normally use
// different physics: the C++ a table of measured amplitudes, the FORTRAN its
// built-in model. A table made FROM that model puts the same physics under both,
// so whatever Tier 3 then reports is a defect or an approximation in the C++
// chain -- table discretisation, interpolation, or the inversion at shifted
// kinematics.
//
// How the table is made. The C++ Born cross section is linear in the amplitudes,
//
//     sigma_Born(phi) = K(kinematics) * (A + Ac cos(phi) + Acc cos(2 phi)),
//
// (the A/Ac/Acc -> H1..H4 round trip already showed this), and the FORTRAN Born
// is exactly B0 + Bc cos(phi) + Bcc cos(2 phi). So at each cell centre:
//
//     A = B0 / K,   Ac / p_t = Bc / (K p_t),   Acc / p_t^2 = Bcc / (K p_t^2)
//
// K comes from the C++ itself, as sigma_Born with a unit table (A = 1,
// Ac = Acc = 0). Both codes then agree on the Born cross section at every cell
// centre by construction -- which `build` checks before writing anything.
//
// Two steps, because the FORTRAN and the C++ library cannot share a process
// (both carry pkhff):
//
//   make_model_table points <bins> <E>                    > cells.txt
//   born_harmonics < cells.txt                            > harm.txt
//   make_model_table build  <bins> <E> harm.txt table.root
//
// <bins> = "nQ2:lo:hi,nnu:lo:hi,nz:lo:hi,npt2:lo:hi"; axes are (Q2, nu, z, pt2).

#include "TRadCor.h"
#include "TSemiInclusiveModel.h"
#include "haprad_constants.h"

#include "TFile.h"
#include "THn.h"
#include "TNamed.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Binning {
  Int_t n[4];
  Double_t lo[4], hi[4];
};

const char* kAxisNames[4] = {"Q2", "nu", "z", "pt2"};

bool ParseBins(const char* spec, Binning& b) {
  return std::sscanf(spec, "%d:%lf:%lf,%d:%lf:%lf,%d:%lf:%lf,%d:%lf:%lf", &b.n[0], &b.lo[0], &b.hi[0], &b.n[1], &b.lo[1],
                     &b.hi[1], &b.n[2], &b.lo[2], &b.hi[2], &b.n[3], &b.lo[3], &b.hi[3]) == 12;
}

Double_t Centre(const Binning& b, Int_t d, Int_t i) { return b.lo[d] + (i - 0.5) * (b.hi[d] - b.lo[d]) / b.n[d]; }

// THadronKinematics prints a line per evaluation; keep it off the terminal.
struct MuteCout {
  std::ostringstream sink;
  std::streambuf* saved;
  MuteCout() : saved(std::cout.rdbuf(sink.rdbuf())) {}
  ~MuteCout() { std::cout.rdbuf(saved); }
};

std::unique_ptr<THnD> NewTable(const char* name, const Binning& b) {
  std::unique_ptr<THnD> h(new THnD(name, name, 4, b.n, b.lo, b.hi));
  for (Int_t d = 0; d < 4; ++d) h->GetAxis(d)->SetName(kAxisNames[d]);
  return h;
}

int Points(const Binning& b, Double_t E) {
  for (Int_t i0 = 1; i0 <= b.n[0]; ++i0)
    for (Int_t i1 = 1; i1 <= b.n[1]; ++i1)
      for (Int_t i2 = 1; i2 <= b.n[2]; ++i2)
        for (Int_t i3 = 1; i3 <= b.n[3]; ++i3) {
          const Double_t Q2 = Centre(b, 0, i0), nu = Centre(b, 1, i1);
          const Double_t z = Centre(b, 2, i2), pt2 = Centre(b, 3, i3);
          const Double_t x = Q2 / (2. * kMassProton * nu);
          if (x >= 1. || nu >= E) continue;  // no such kinematics
          printf("%d %d %d %d %.10g %.12g %.12g %.12g %.12g\n", i0, i1, i2, i3, E, x, Q2, z, std::sqrt(pt2));
        }
  return 0;
}

int Build(const Binning& b, Double_t E, const char* harmPath, const char* outPath) {
  // FORTRAN harmonics, keyed by cell. PDFLIB prints a banner to stdout, so keep
  // only lines that parse as data.
  std::map<std::vector<Int_t>, std::vector<Double_t>> harm;
  {
    std::ifstream in(harmPath);
    std::string line;
    while (std::getline(in, line)) {
      std::istringstream ls(line);
      Int_t j[4], ok;
      std::string s0, s1, s2;
      if (!(ls >> j[0] >> j[1] >> j[2] >> j[3] >> ok >> s0 >> s1 >> s2)) continue;
      for (std::string* s : {&s0, &s1, &s2})
        for (char& c : *s)
          if (c == 'D' || c == 'd') c = 'e';  // FORTRAN exponents
      if (ok != 1) continue;
      harm[{j[0], j[1], j[2], j[3]}] = {std::atof(s0.c_str()), std::atof(s1.c_str()), std::atof(s2.c_str())};
    }
  }
  if (harm.empty()) {
    fprintf(stderr, "no usable cells in %s\n", harmPath);
    return 1;
  }

  // Unit table, for K. One cell per axis, wide enough to cover everything.
  const TString unitPath = TString(outPath) + ".unit.root";
  {
    Binning u = {{1, 1, 1, 1}, {0., 0., 0., 0.}, {1000., 1000., 1.0001, 1000.}};
    auto a = NewTable("A", u), ac = NewTable("Ac", u), acc = NewTable("Acc", u);
    Int_t one[4] = {1, 1, 1, 1};
    a->SetBinContent(one, 1.);
    TFile f(unitPath, "RECREATE");
    a->Write();
    ac->Write();
    acc->Write();
    TNamed("pt_scaling", "raw").Write();
  }

  auto A = NewTable("A", b), Ac = NewTable("Ac", b), Acc = NewTable("Acc", b), fitted = NewTable("fitted", b);
  Long64_t nCells = 0, nBuilt = 0, nNoK = 0;
  {
    TRadCor rc;
    if (!rc.LoadSemiInclusiveTable(unitPath)) return 1;
    MuteCout mute;
    for (const auto& kv : harm) {
      ++nCells;
      const std::vector<Int_t>& j = kv.first;
      const Double_t Q2 = Centre(b, 0, j[0]), nu = Centre(b, 1, j[1]);
      const Double_t z = Centre(b, 2, j[2]), pt2 = Centre(b, 3, j[3]);
      const Double_t x = Q2 / (2. * kMassProton * nu), pt = std::sqrt(pt2);
      const Double_t K = rc.CalculateBorn(E, x, Q2, z, pt, 90.);
      if (!(K > 0.) || !std::isfinite(K)) {
        ++nNoK;
        continue;
      }
      Int_t idx[4] = {j[0], j[1], j[2], j[3]};
      A->SetBinContent(idx, kv.second[0] / K);
      Ac->SetBinContent(idx, kv.second[1] / (K * pt));
      Acc->SetBinContent(idx, kv.second[2] / (K * pt2));
      fitted->SetBinContent(idx, 1.);
      ++nBuilt;
    }
  }
  std::remove(unitPath.Data());

  {
    TFile f(outPath, "RECREATE");
    A->Write();
    Ac->Write();
    Acc->Write();
    fitted->Write();
    TNamed("pt_scaling", "reduced").Write();
    const char* interp = getenv("MODEL_TABLE_INTERP");
    TNamed("interpolation", interp ? interp : "linear").Write();
    TNamed("source", Form("HAPRAD 2.0 PDF x FF model (GRV94 LO x PKH), Born harmonics, E = %g GeV", E)).Write();
  }
  printf("cells with FORTRAN harmonics : %lld\n", nCells);
  printf("  built                      : %lld\n", nBuilt);
  printf("  C++ rejects the kinematics : %lld\n", nNoK);

  // Self-check: through the real reader, the C++ Born at each cell centre must
  // reproduce the FORTRAN harmonics at any phi.
  TRadCor rc;
  if (!rc.LoadSemiInclusiveTable(outPath)) return 1;
  Double_t worst = 0.;
  Long64_t nChecked = 0;
  struct Res { Double_t d, Q2, nu, z, pt2, phi; };
  std::vector<Res> res;
  {
    MuteCout mute;
    for (const auto& kv : harm) {
      const std::vector<Int_t>& j = kv.first;
      Int_t idx[4] = {j[0], j[1], j[2], j[3]};
      if (fitted->GetBinContent(idx) == 0.) continue;
      const Double_t Q2 = Centre(b, 0, j[0]), nu = Centre(b, 1, j[1]);
      const Double_t z = Centre(b, 2, j[2]), pt2 = Centre(b, 3, j[3]);
      const Double_t x = Q2 / (2. * kMassProton * nu);
      for (Double_t phi : {0., 70., 145., 180.}) {
        const Double_t f = kv.second[0] + kv.second[1] * std::cos(phi / kRadianDeg) + kv.second[2] * std::cos(2. * phi / kRadianDeg);
        const Double_t c = rc.CalculateBorn(E, x, Q2, z, std::sqrt(pt2), phi);
        if (f == 0.) continue;
        const Double_t d = std::fabs(c - f) / std::fabs(f);
        worst = std::max(worst, d);
        res.push_back({d, Q2, nu, z, pt2, phi});
        ++nChecked;
      }
    }
  }
  printf("self-check: C++ Born vs FORTRAN Born at %lld (cell, phi) pairs, worst reldiff %.2e\n", nChecked, worst);
  if (getenv("MODEL_TABLE_DIAG")) {
    std::sort(res.begin(), res.end(), [](const Res& a, const Res& b) { return a.d > b.d; });
    printf("  %-10s %6s %6s %5s %6s %6s %6s\n", "reldiff", "Q2", "nu", "z", "pt2", "y", "phi");
    for (size_t i = 0; i < res.size() && i < 8; ++i)
      printf("  %-10.2e %6.2f %6.2f %5.2f %6.3f %6.3f %6.0f\n", res[i].d, res[i].Q2, res[i].nu, res[i].z, res[i].pt2,
             res[i].nu / E, res[i].phi);
    // Residual against y, binned.
    Double_t ysum[5] = {0}, ymax[5] = {0};
    Long64_t yn[5] = {0};
    for (const Res& r : res) {
      const Int_t k = std::min(4, Int_t(5. * r.nu / E));
      ysum[k] += r.d; ymax[k] = std::max(ymax[k], r.d); ++yn[k];
    }
    printf("  y bin       mean reldiff   max reldiff\n");
    for (Int_t k = 0; k < 5; ++k)
      if (yn[k]) printf("  %.1f-%.1f    %10.2e   %10.2e\n", 0.2 * k, 0.2 * (k + 1), ysum[k] / yn[k], ymax[k]);
  }
  printf("wrote %s\n", outPath);
  // Not exactly zero: the C++ inversion builds lambda = S X Q2 - M^2 Q2^2 and
  // drops the -m_e^2 lambda_q term (TStructFunctionArray, "electron part is
  // ignored"), while the Born uses the full lambda in b^2. The neglected
  // fraction m_e^2 lambda_q / lambda grows with y and falls with Q2: measured
  // 1e-11 below y = 0.2, 1.6e-6 at (Q2, y) = (1.03, 0.94), 1.0e-5 at
  // (0.43, 0.98). A real disagreement would sit far above 1e-4.
  return worst < 1e-4 ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  Binning b;
  if (argc >= 4 && std::string(argv[1]) == "points" && ParseBins(argv[2], b)) return Points(b, std::atof(argv[3]));
  if (argc >= 6 && std::string(argv[1]) == "build" && ParseBins(argv[2], b))
    return Build(b, std::atof(argv[3]), argv[4], argv[5]);
  fprintf(stderr,
          "usage: make_model_table points <bins> <E>                > cells.txt\n"
          "       born_harmonics < cells.txt                        > harm.txt\n"
          "       make_model_table build  <bins> <E> harm.txt table.root\n"
          "  <bins> = nQ2:lo:hi,nnu:lo:hi,nz:lo:hi,npt2:lo:hi\n");
  return 2;
}
