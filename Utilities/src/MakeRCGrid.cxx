// MakeRCGrid -- HAPRAD on a grid of nodes, for weighting events (PLAN.md Phase 3).
//
//   MakeRCGrid run   <config> <table.root> <chunk> <nchunks> <grid_chunk.root>
//   hadd grid.root grid_chunk*.root
//   MakeRCGrid check <config> <grid.root> <table.root> <n> <chunk> <nchunks> <ntuple files ...>
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
// ROOT's error handler, hundreds of times per HAPRAD call (also at the test
// points that converge in Tier 3). Printed, it would run to gigabytes over a
// grid; it is counted instead, per node.
Long64_t gIntegratorWarnings = 0;
void CountGSLWarnings(int level, Bool_t abort, const char* location, const char* msg) {
  if (location && std::string(location).find("GSLError") != std::string::npos) {
    ++gIntegratorWarnings;
    return;
  }
  DefaultErrorHandler(level, abort, location, msg);
}

struct Result {
  bool physical = false, valid = false;
  Double_t sib = 0., r[3] = {0., 0., 0.};
  Double_t outside = 0., empty = 0.;  // fraction of table lookups
  Long64_t warnings = 0;              // integrator warnings during the call
};

// One HAPRAD call, with the library's per-call printout swallowed.
Result Evaluate(TRadCor& rc, const Config& c, Double_t Q2, Double_t nu, Double_t z, Double_t pt, Double_t phi) {
  Result res;
  const Double_t M = kMassProton, m = kMassPion;
  const Double_t x = Q2 / (2. * M * nu), W2 = M * M + 2. * M * nu - Q2, Eh = z * nu;
  pt = std::max(pt, kPtMin);
  res.physical = x < 1. && nu < c.E && W2 > (M + m) * (M + m) && Eh > m && pt * pt < Eh * Eh - m * m;
  if (!res.physical) return res;

  std::ostringstream sink;
  std::streambuf* out = std::cout.rdbuf(sink.rdbuf());
  std::streambuf* err = std::cerr.rdbuf(sink.rdbuf());
  const Long64_t before = gIntegratorWarnings;
  rc.CalculateRCFactor(c.E, x, Q2, z, pt, phi, kMx2Threshold, c.targetNAZ);
  res.warnings = gIntegratorWarnings - before;
  std::cout.rdbuf(out);
  std::cerr.rdbuf(err);

  res.sib = rc.GetSigBorn();
  if (res.sib > 0. && std::isfinite(res.sib)) {
    res.r[0] = rc.GetSigObs() / res.sib;
    res.r[1] = rc.GetTail(0) / res.sib;
    res.r[2] = rc.GetTail(1) / res.sib;
    res.valid = std::isfinite(res.r[0]) && std::isfinite(res.r[1]) && std::isfinite(res.r[2]);
  }
  if (const TSemiInclusiveModel* sf = rc.GetSemiInclusiveModel())
    if (sf->NLookups() > 0) {
      res.outside = double(sf->NOutOfRange()) / sf->NLookups();
      res.empty = double(sf->NEmptyCell()) / sf->NLookups();
    }
  return res;
}

std::string TableLuminosity(const char* tablePath) {
  std::unique_ptr<TFile> f(TFile::Open(tablePath, "READ"));
  if (!f || f->IsZombie()) return "";
  TNamed* t = f->Get<TNamed>("luminosity");
  return t ? t->GetTitle() : "absent (not a MakePhiTable table; taken as absolute)";
}

int Run(const Config& c, const char* tablePath, Long64_t chunk, Long64_t nChunks, const char* outPath) {
  const rge::RCNodes g(c);
  TRadCor rc;
  if (!rc.LoadSemiInclusiveTable(tablePath)) return 1;

  std::unique_ptr<THnD> sib(g.MakeHist("sigma_born")), comp[3], valid(g.MakeHist("valid")),
      computed(g.MakeHist("computed")), outside(g.MakeHist("table_outside")), empty(g.MakeHist("table_empty")),
      warn(g.MakeHist("integrator_warnings"));
  for (Int_t k = 0; k < 3; ++k) comp[k].reset(g.MakeHist(rge::kRCComp[k]));

  Long64_t nNodes = 0, nUnphysical = 0, nRejected = 0, nValid = 0;
  TStopwatch clock;
  const Long64_t size = g.Size();
  Int_t i[5], bin[5];
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

  // Before creating the output: opening and closing the table resets ROOT's
  // current directory, and the writes below would go nowhere.
  const std::string tableLumi = TableLuminosity(tablePath);
  TFile out(outPath, "RECREATE");
  if (out.IsZombie()) return 1;
  for (THnD* h : {sib.get(), comp[0].get(), comp[1].get(), comp[2].get(), valid.get(), computed.get(), outside.get(),
                  empty.get(), warn.get()})
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

int Check(const Config& c, const char* gridPath, const char* tablePath, Long64_t n, Long64_t chunk, Long64_t nChunks,
          int nIn, char** in) {
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

  TChain chain("DT");
  for (int k = 0; k < nIn; ++k) chain.Add(in[k]);
  rge::Event e;
  if (!rge::BindEvent(chain, e)) return 1;

  // Reservoir sample of selected events inside the grid: the same sample in
  // every chunk.
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

}  // namespace

int main(int argc, char** argv) {
  const std::string mode = argc > 1 ? argv[1] : "";
  SetErrorHandler(CountGSLWarnings);
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
          "Chunks of `run` merge with `hadd`; run_rc_grid.sh does both.\n");
  return 2;
}
