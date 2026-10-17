// MakeRCGrid -- HAPRAD evaluated on a regular grid of nodes in (Q2, nu, z, pt,
// phi), from which ApplyRC interpolates an RC factor for every event; plus the
// check of that interpolation and the comparison of two grids.
//
//   MakeRCGrid run   <config> <table.root> <chunk> <nchunks> <grid_chunk.root>
//   hadd grid.root grid_chunk*.root
//   MakeRCGrid check <config> <grid.root> <table.root> <n> <chunk> <nchunks> <ntuple files ...>
//   MakeRCGrid compare <gridA.root> <gridB.root> [n]
//   MakeRCGrid inspect <grid.root> [n]
//
// `run` evaluates every nchunks-th node, starting at `chunk` -- interleaved, so
// chunks cost about the same. One HAPRAD call takes ~0.5 s; run_rc_grid.sh
// runs the chunks in parallel and merges them. The node grid (rc_nodes_* keys),
// beam energy and particle come from the same config as the table.
//
// `check` is the interpolation error budget: it samples n selected events
// (fixed seed, so every chunk sees the same sample), calls HAPRAD directly at
// each one's own kinematics, and prints that next to the interpolated grid.
// check_rc_grid.py runs the chunks and summarises.
//
// `compare` is the convergence measure of the RC iteration (rc_iterate.py):
// the change of the RC factor without the exclusive tail, r_vr + r_in, over the
// nodes valid in both grids. It prints one machine-readable line,
//   @compare <nodes> <max |d|> <median |d|> <95% |d|> <mean d> <valid in one only>
// and, given n, the n nodes that changed most.
//
// `inspect` shows the distribution of the RC factor over one grid's nodes and
// the n nodes furthest from 1, with what HAPRAD saw there.
//
// The grid belongs to ONE table: every iteration of the RC loop needs a new one.

#include "RCGrid.hxx"
#include "RGESelection.hxx"
#include "TRadCor.h"
#include "TSemiInclusiveModel.h"
#include "haprad_constants.h"

#include "TChain.h"
#include "TError.h"
#include "TFile.h"
#include "THn.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "TStopwatch.h"

#include <algorithm>
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

// HAPRAD's lower bound on M_x^2 for the semi-inclusive part, as GetRC uses it:
// the pi+ n threshold.
const Double_t kMx2Threshold = (kMassNeutron + kMassPion) * (kMassNeutron + kMassPion);

// p_t = 0 makes the C++ return sigma_Born = 0; at 1 MeV the RC factor has
// already lost its phi dependence (checked: 0.94573 at phi = 0, 90, 180).
const Double_t kPtMin = 1e-3;

// GSL's "failed to reach tolerance" from the inner integrals arrives through
// ROOT's error handler, hundreds of times per HAPRAD call, also at points
// whose result agrees with the FORTRAN. Printed, it would run to gigabytes over
// a grid; it is counted instead, per node. Every other message goes to ROOT's
// default handler.
Long64_t gIntegratorWarnings = 0;
void CountGSLWarnings(int level, Bool_t abort, const char* location, const char* msg) {
  if (location && std::string(location).find("GSLError") != std::string::npos) {
    ++gIntegratorWarnings;
    return;
  }
  DefaultErrorHandler(level, abort, location, msg);
}

// What one HAPRAD call gives: whether the node is kinematically allowed,
// whether HAPRAD returned a finite result, sigma_Born and the three ratios
// r_vr, r_in, r_ex (RCGrid.hxx), and the table-lookup diagnostics.
struct Result {
  bool physical = false, valid = false;
  Double_t sib = 0., r[3] = {0., 0., 0.};
  Double_t outside = 0., empty = 0., second = 0.;  // fraction of table lookups
  Long64_t warnings = 0;              // integrator warnings during the call
};

// One HAPRAD call at (Q2, nu, z, pt, phi[deg]), with the library's per-call
// printout swallowed. Kinematically forbidden points (x >= 1, nu >= E, W below
// the pion threshold, E_h <= m_pi, pt beyond the hadron momentum) are not sent
// to HAPRAD. pt is raised to kPtMin. The ratios are formed only for a positive,
// finite sigma_Born.
Result Evaluate(TRadCor& rc, const Config& c, Double_t Q2, Double_t nu, Double_t z, Double_t pt, Double_t phi) {
  Result res;
  const Double_t M = kMassProton, m = kMassPion;
  const Double_t x = Q2 / (2. * M * nu), W2 = M * M + 2. * M * nu - Q2, Eh = z * nu;
  // ---- Kinematic limits: forbidden nodes are not sent to HAPRAD ----
  pt = std::max(pt, kPtMin);
  res.physical = x < 1. && nu < c.E && W2 > (M + m) * (M + m) && Eh > m && pt * pt < Eh * Eh - m * m;
  if (!res.physical) return res;

  // ---- The HAPRAD call, its printout swallowed, integrator warnings counted ----
  std::ostringstream sink;
  std::streambuf* out = std::cout.rdbuf(sink.rdbuf());
  std::streambuf* err = std::cerr.rdbuf(sink.rdbuf());
  const Long64_t before = gIntegratorWarnings;
  rc.CalculateRCFactor(c.E, x, Q2, z, pt, phi, kMx2Threshold, c.targetNAZ);
  res.warnings = gIntegratorWarnings - before;
  std::cout.rdbuf(out);
  std::cerr.rdbuf(err);

  // ---- The three ratios to sigma_Born ----
  res.sib = rc.GetSigBorn();
  if (res.sib > 0. && std::isfinite(res.sib)) {
    res.r[0] = rc.GetSigObs() / res.sib;
    res.r[1] = rc.GetTail(0) / res.sib;
    res.r[2] = rc.GetTail(1) / res.sib;
    res.valid = std::isfinite(res.r[0]) && std::isfinite(res.r[1]) && std::isfinite(res.r[2]);
  }
  // ---- Table-lookup diagnostics of this call ----
  if (const TSemiInclusiveModel* sf = rc.GetSemiInclusiveModel())
    if (sf->NLookups() > 0) {
      res.outside = double(sf->NOutOfRange()) / sf->NLookups();
      res.empty = double(sf->NEmptyCell()) / sf->NLookups();
      res.second = double(sf->NSecondNeighbour()) / sf->NLookups();
    }
  return res;
}

// The table's "luminosity" key, copied into the grid so that ApplyRC can tell
// whether the exclusive tail is on the table's scale. Tables without the key
// (not made by MakePhiTable, e.g. model tables) are absolute.
std::string TableLuminosity(const char* tablePath) {
  std::unique_ptr<TFile> f(TFile::Open(tablePath, "READ"));
  if (!f || f->IsZombie()) return "";
  TNamed* t = f->Get<TNamed>("luminosity");
  return t ? t->GetTitle() : "absent (not a MakePhiTable table; taken as absolute)";
}

// The `run` mode: evaluate nodes chunk, chunk + nChunks, ... and write their
// results into one-bin-per-node histograms (RCGrid.hxx). Every node of the
// chunk is marked in "computed", valid or not, so that the merge can be
// checked for completeness. Per node it also stores the fractions of table
// lookups that were out of range, empty or second-neighbour, and the number of
// integrator warnings.
int Run(const Config& c, const char* tablePath, Long64_t chunk, Long64_t nChunks, const char* outPath) {
  // ==== Load the table and book one histogram per stored quantity ====
  const rge::RCNodes g(c);
  TRadCor rc;
  if (!rc.LoadSemiInclusiveTable(tablePath)) return 1;

  std::unique_ptr<THnD> sib(g.MakeHist("sigma_born")), comp[3], valid(g.MakeHist("valid")),
      computed(g.MakeHist("computed")), outside(g.MakeHist("table_outside")), empty(g.MakeHist("table_empty")),
      second(g.MakeHist("table_second")),
      warn(g.MakeHist("integrator_warnings"));
  for (Int_t k = 0; k < 3; ++k) comp[k].reset(g.MakeHist(rge::kRCComp[k]));

  Long64_t nNodes = 0, nUnphysical = 0, nRejected = 0, nValid = 0;
  TStopwatch clock;
  const Long64_t size = g.Size();
  Int_t i[5], bin[5];
  // ==== Evaluate this chunk's nodes ====
  // Every nChunks-th node; each is marked computed, then filled if valid.
  for (Long64_t k = chunk; k < size; k += nChunks) {
    g.Decode(k, i);
    for (Int_t d = 0; d < 5; ++d) bin[d] = i[d] + 1;
    const Result r = Evaluate(rc, c, g.Node(0, i[0]), g.Node(1, i[1]), g.Node(2, i[2]), g.Node(3, i[3]), g.Node(4, i[4]));
    ++nNodes;
    computed->SetBinContent(bin, 1.);
    if (!r.physical) {
      ++nUnphysical;
      continue;
    }
    outside->SetBinContent(bin, r.outside);
    empty->SetBinContent(bin, r.empty);
    second->SetBinContent(bin, r.second);
    warn->SetBinContent(bin, r.warnings);
    if (!r.valid) {
      ++nRejected;
      continue;
    }
    ++nValid;
    valid->SetBinContent(bin, 1.);
    sib->SetBinContent(bin, r.sib);
    for (Int_t m = 0; m < 3; ++m) comp[m]->SetBinContent(bin, r.r[m]);
    if (nValid % 500 == 0) {
      fprintf(stderr, "chunk %lld: %lld nodes, %.0f s\n", chunk, nNodes, clock.RealTime());
      clock.Continue();
    }
  }

  // ==== Write the chunk ====
  // Before creating the output: opening and closing the table resets ROOT's
  // current directory, and the writes below would go nowhere.
  const std::string tableLumi = TableLuminosity(tablePath);
  TFile out(outPath, "RECREATE");
  if (out.IsZombie()) return 1;
  for (THnD* h : {sib.get(), comp[0].get(), comp[1].get(), comp[2].get(), valid.get(), computed.get(), outside.get(),
                  empty.get(), second.get(), warn.get()})
    h->Write();
  TNamed("table", tablePath).Write();
  TNamed("table_luminosity", tableLumi.c_str()).Write();
  TNamed("beam_energy", Form("%.6g", c.E)).Write();
  TNamed("mx2_threshold", Form("%.6g", kMx2Threshold)).Write();
  TNamed("config", c.text.c_str()).Write();
  out.Close();

  printf("%s: chunk %lld of %lld, %lld nodes -- %lld unphysical, %lld rejected by HAPRAD, %lld valid; "
         "%lld integrator warnings; %.0f s\n",
         outPath, chunk, nChunks, nNodes, nUnphysical, nRejected, nValid, gIntegratorWarnings, clock.RealTime());
  return 0;
}

// The `check` mode: draw a reservoir sample of n selected events inside the
// grid, then for every nChunks-th one (starting at `chunk`) print the
// interpolated ratios next to a direct HAPRAD call at the event's kinematics,
// one "@" line per event.
int Check(const Config& c, const char* gridPath, const char* tablePath, Long64_t n, Long64_t chunk, Long64_t nChunks,
          int nIn, char** in) {
  // ==== Load the grid and the table it was built from ====
  rge::RCGrid grid;
  std::string why;
  if (!grid.Load(gridPath, why)) {
    fprintf(stderr, "%s\n", why.c_str());
    return 1;
  }
  if (grid.Info("table") != tablePath)
    fprintf(stderr, "WARNING: the grid was built from %s, not %s\n", grid.Info("table").c_str(), tablePath);
  TRadCor rc;
  if (!rc.LoadSemiInclusiveTable(tablePath)) return 1;

  // ==== Bind the input ntuples ====
  TChain chain("DT");
  for (int k = 0; k < nIn; ++k) chain.Add(in[k]);
  rge::Event e;
  if (!rge::BindEvent(chain, e)) return 1;

  // ==== Draw the event sample ====
  // Reservoir sample of selected events inside the grid (each event replaces
  // a random member with probability n / seen), with a fixed seed: the same
  // sample in every chunk.
  std::vector<rge::Event> sample;
  TRandom3 rng(4357);
  Long64_t seen = 0;
  const Long64_t entries = chain.GetEntries();
  for (Long64_t k = 0; k < entries; ++k) {
    chain.GetEntry(k);
    if (rge::Stages(c, e) < rge::kNStages) continue;
    Double_t v[5], r[3], cover;
    rge::GridPoint(e, v);
    if (grid.Interpolate(v, r, cover) == rge::RCGrid::kOutside) continue;
    ++seen;
    if (Long64_t(sample.size()) < n)
      sample.push_back(e);
    else {
      const Long64_t j = Long64_t(rng.Rndm() * seen);
      if (j < n) sample[j] = e;
    }
  }

  // ==== Interpolated vs direct, for this chunk's share of the sample ====
  printf("# id Q2 nu z pt phi status cover  interp: r_vr r_in r_ex  direct: r_vr r_in r_ex\n");
  for (Long64_t j = chunk; j < Long64_t(sample.size()); j += nChunks) {
    Double_t v[5], r[3], cover;
    rge::GridPoint(sample[j], v);
    const int status = grid.Interpolate(v, r, cover);
    const Result d = Evaluate(rc, c, v[0], v[1], v[2], v[3], v[4]);
    if (!d.valid) {
      printf("@ %lld %.5g %.5g %.5g %.5g %.5g %d %.4f  no direct result\n", j, v[0], v[1], v[2], v[3], v[4], status,
             cover);
      continue;
    }
    printf("@ %lld %.5g %.5g %.5g %.5g %.5g %d %.4f  %.8g %.8g %.8g  %.8g %.8g %.8g\n", j, v[0], v[1], v[2], v[3], v[4],
           status, cover, r[0], r[1], r[2], d.r[0], d.r[1], d.r[2]);
    fflush(stdout);
  }
  return 0;
}

// The kinematics of a grid node from its global bin, as "Q2 nu z pt phi x".
std::string NodeKinematics(const THnD* h, Long64_t bin) {
  Int_t idx[5];
  h->GetBinContent(bin, idx);
  Double_t v[5];
  for (Int_t d = 0; d < 5; ++d) v[d] = h->GetAxis(d)->GetBinCenter(idx[d]);
  return Form("Q2 %5.2f nu %5.2f z %4.2f pt %5.3f phi %5.1f  x %5.3f", v[0], v[1], v[2], v[3], v[4],
              v[0] / (2. * kMassProton * v[1]));
}

// The `compare` mode: over the nodes valid in both grids, the change d of
// r_vr + r_in from grid A to grid B; prints the count, max, median and 95th
// percentile of |d|, the mean d, and the number of nodes valid in only one.
// With nWorst > 0, also lists the nWorst nodes that changed most.
int Compare(const char* pathA, const char* pathB, Int_t nWorst) {
  std::unique_ptr<TFile> fa(TFile::Open(pathA, "READ")), fb(TFile::Open(pathB, "READ"));
  if (!fa || fa->IsZombie() || !fb || fb->IsZombie()) {
    fprintf(stderr, "cannot open %s or %s\n", pathA, pathB);
    return 1;
  }
  // ==== Read both grids; they must have the same nodes ====
  const char* names[3] = {"valid", "r_vr", "r_in"};
  THnD *a[3], *b[3];
  for (Int_t k = 0; k < 3; ++k) {
    a[k] = fa->Get<THnD>(names[k]);
    b[k] = fb->Get<THnD>(names[k]);
    if (!a[k] || !b[k] || a[k]->GetNbins() != b[k]->GetNbins()) {
      fprintf(stderr, "%s and %s are not RC grids on the same nodes\n", pathA, pathB);
      return 1;
    }
  }
  for (Int_t d = 0; d < 5; ++d)
    if (a[0]->GetAxis(d)->GetNbins() != b[0]->GetAxis(d)->GetNbins() ||
        std::fabs(a[0]->GetAxis(d)->GetXmin() - b[0]->GetAxis(d)->GetXmin()) > 1e-9 ||
        std::fabs(a[0]->GetAxis(d)->GetXmax() - b[0]->GetAxis(d)->GetXmax()) > 1e-9) {
      fprintf(stderr, "%s and %s have different nodes\n", pathA, pathB);
      return 1;
    }
  // ==== Change of r_vr + r_in over the nodes valid in both ====
  std::vector<Double_t> diff;
  std::vector<std::pair<Double_t, Long64_t>> byNode;  // (|d|, bin), for the listing
  Long64_t oneOnly = 0;
  Double_t sum = 0.;
  for (Long64_t i = 0; i < a[0]->GetNbins(); ++i) {
    const bool va = a[0]->GetBinContent(i) > 0.5, vb = b[0]->GetBinContent(i) > 0.5;
    if (va != vb) ++oneOnly;
    if (!va || !vb) continue;
    const Double_t d = (b[1]->GetBinContent(i) + b[2]->GetBinContent(i)) - (a[1]->GetBinContent(i) + a[2]->GetBinContent(i));
    diff.push_back(std::fabs(d));
    byNode.emplace_back(std::fabs(d), i);
    sum += d;
  }
  if (diff.empty()) {
    fprintf(stderr, "no node is valid in both grids\n");
    return 1;
  }
  // ==== Summary statistics, one machine-readable line ====
  const Double_t mean = sum / diff.size();
  std::sort(diff.begin(), diff.end());
  auto q = [&](Double_t f) { return diff[std::min(diff.size() - 1, size_t(f * (diff.size() - 1) + 0.5))]; };
  printf("@compare %zu %.6g %.6g %.6g %+.6g %lld\n", diff.size(), diff.back(), q(0.5), q(0.95), mean, oneOnly);

  // ==== Optional: the nodes that changed most ====
  if (nWorst > 0) {
    std::sort(byNode.begin(), byNode.end(), [](const auto& p, const auto& r) { return p.first > r.first; });
    printf("largest changes of r_vr + r_in (A -> B):\n");
    for (size_t k = 0; k < byNode.size() && k < size_t(nWorst); ++k) {
      const Long64_t i = byNode[k].second;
      printf("  %s   %10.4g -> %10.4g\n", NodeKinematics(a[0], i).c_str(), a[1]->GetBinContent(i) + a[2]->GetBinContent(i),
             b[1]->GetBinContent(i) + b[2]->GetBinContent(i));
    }
  }
  return 0;
}

// The `inspect` mode: how the RC factor r_vr + r_in is distributed over the
// valid nodes of one grid, and the nWorst nodes furthest from 1 (any factor
// <= 0 first, then by |log rc|), each with its kinematics, sigma_Born, the
// three ratios and the fractions of its table lookups that were outside the
// table, empty, or filled from second neighbours. Tells noisy table cells
// (small sigma_Born) from table edges (many outside / empty lookups).
int Inspect(const char* path, Int_t nWorst) {
  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  if (!f || f->IsZombie()) {
    fprintf(stderr, "cannot open %s\n", path);
    return 1;
  }
  // ==== Read the stored quantities ====
  const char* names[8] = {"valid", "sigma_born", "r_vr", "r_in", "r_ex", "table_outside", "table_empty", "table_second"};
  THnD* h[8];
  for (Int_t k = 0; k < 8; ++k) {
    h[k] = f->Get<THnD>(names[k]);
    if (!h[k]) {
      fprintf(stderr, "%s has no %s: not an RC grid\n", path, names[k]);
      return 1;
    }
  }
  // ==== The RC factor at every valid node ====
  std::vector<std::pair<Double_t, Long64_t>> rc;  // (r_vr + r_in, bin)
  for (Long64_t i = 0; i < h[0]->GetNbins(); ++i)
    if (h[0]->GetBinContent(i) > 0.5) rc.emplace_back(h[2]->GetBinContent(i) + h[3]->GetBinContent(i), i);
  if (rc.empty()) {
    fprintf(stderr, "%s has no valid node\n", path);
    return 1;
  }
  std::vector<Double_t> v;
  for (const auto& p : rc) v.push_back(p.first);
  std::sort(v.begin(), v.end());
  auto q = [&](Double_t x) { return v[std::min(v.size() - 1, size_t(x * (v.size() - 1) + 0.5))]; };
  auto count = [&](auto pred) { return std::count_if(v.begin(), v.end(), pred); };
  printf("%zu valid nodes; RC factor r_vr + r_in:\n", v.size());
  printf("  min %.4g  1%% %.4g  5%% %.4g  median %.4g  95%% %.4g  99%% %.4g  max %.4g\n", v.front(), q(0.01), q(0.05),
         q(0.5), q(0.95), q(0.99), v.back());
  printf("  <= 0: %ld   < 0.5: %ld   > 1.5: %ld   > 2: %ld   > 10: %ld\n", count([](Double_t x) { return x <= 0.; }),
         count([](Double_t x) { return x < 0.5; }), count([](Double_t x) { return x > 1.5; }),
         count([](Double_t x) { return x > 2.; }), count([](Double_t x) { return x > 10.; }));

  // ==== Table lookups: all nodes against those outside [0.5, 2] ====
  Double_t sAll[3] = {0., 0., 0.}, sOdd[3] = {0., 0., 0.};
  Long64_t nOdd = 0;
  for (const auto& p : rc) {
    const bool odd = !(p.first >= 0.5 && p.first <= 2.);
    for (Int_t k = 0; k < 3; ++k) {
      const Double_t x = h[5 + k]->GetBinContent(p.second);
      sAll[k] += x;
      if (odd) sOdd[k] += x;
    }
    nOdd += odd;
  }
  printf("mean fraction of table lookups   outside   empty   2nd nbr\n");
  printf("  all valid nodes               %7.3f %7.3f %7.3f\n", sAll[0] / rc.size(), sAll[1] / rc.size(), sAll[2] / rc.size());
  if (nOdd)
    printf("  nodes with rc outside [0.5,2] %7.3f %7.3f %7.3f   (%lld nodes)\n", sOdd[0] / nOdd, sOdd[1] / nOdd,
           sOdd[2] / nOdd, nOdd);

  // ==== The nodes furthest from 1 ====
  auto badness = [](Double_t x) { return x <= 0. ? 1e300 : std::fabs(std::log(x)); };
  std::sort(rc.begin(), rc.end(), [&](const auto& p, const auto& r) { return badness(p.first) > badness(r.first); });
  printf("nodes furthest from rc = 1:\n");
  printf("  %-55s %10s %9s %9s %9s %9s  %6s %6s %6s\n", "node", "sigma_Born", "rc", "r_vr", "r_in", "r_ex", "out", "empty",
         "2nd");
  for (size_t k = 0; k < rc.size() && k < size_t(nWorst); ++k) {
    const Long64_t i = rc[k].second;
    printf("  %-55s %10.4g %9.4g %9.4g %9.4g %9.4g  %6.3f %6.3f %6.3f\n", NodeKinematics(h[0], i).c_str(),
           h[1]->GetBinContent(i), rc[k].first, h[2]->GetBinContent(i), h[3]->GetBinContent(i), h[4]->GetBinContent(i),
           h[5]->GetBinContent(i), h[6]->GetBinContent(i), h[7]->GetBinContent(i));
  }
  return 0;
}

}  // namespace

// Parse the mode and its arguments, install the GSL-warning counter, check
// that the config has an RC grid and that the MAID file is reachable, then
// dispatch.
int main(int argc, char** argv) {
  const std::string mode = argc > 1 ? argv[1] : "";
  SetErrorHandler(CountGSLWarnings);
  if (mode == "compare" && (argc == 4 || argc == 5)) return Compare(argv[2], argv[3], argc == 5 ? atoi(argv[4]) : 0);
  if (mode == "inspect" && (argc == 3 || argc == 4)) return Inspect(argv[2], argc == 4 ? atoi(argv[3]) : 20);
  Config c;
  if ((mode == "run" && argc == 7) || (mode == "check" && argc >= 9)) {
    if (!rge::LoadConfig(argv[2], c)) return 1;
    // HAPRAD reads its MAID grid from the working directory (run_rc_grid.sh
    // and check_rc_grid.py change into one that has it).
    if (!std::ifstream("pi_n_maid.dat")) {
      fprintf(stderr,
              "pi_n_maid.dat is not in the working directory; the exclusive tail needs it. Run from a "
              "directory that has it (haprad2/), or use run_rc_grid.sh / check_rc_grid.py.\n");
      return 1;
    }
    if (!c.hasRCGrid) {
      fprintf(stderr, "%s has no rc_nodes_* keys\n", argv[2]);
      return 1;
    }
    if (mode == "run") {
      const Long64_t chunk = atoll(argv[4]), nChunks = atoll(argv[5]);
      if (nChunks < 1 || chunk < 0 || chunk >= nChunks) {
        fprintf(stderr, "need 0 <= chunk < nchunks\n");
        return 2;
      }
      return Run(c, argv[3], chunk, nChunks, argv[6]);
    }
    const Long64_t chunk = atoll(argv[6]), nChunks = atoll(argv[7]);
    if (nChunks < 1 || chunk < 0 || chunk >= nChunks) {
      fprintf(stderr, "need 0 <= chunk < nchunks\n");
      return 2;
    }
    return Check(c, argv[3], argv[4], atoll(argv[5]), chunk, nChunks, argc - 8, argv + 8);
  }
  fprintf(stderr,
          "usage: MakeRCGrid run   <config> <table.root> <chunk> <nchunks> <grid_chunk.root>\n"
          "       MakeRCGrid check <config> <grid.root> <table.root> <n> <chunk> <nchunks> <ntuple files ...>\n"
          "       MakeRCGrid compare <gridA.root> <gridB.root> [n: list the n largest changes]\n"
          "       MakeRCGrid inspect <grid.root> [n: list the n nodes furthest from 1, default 20]\n"
          "Chunks of `run` merge with `hadd`; run_rc_grid.sh does both.\n");
  return 2;
}
