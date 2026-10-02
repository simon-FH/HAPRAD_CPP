// RGESelection.hxx -- the analysis configuration and event selection shared by
// MakePhiTable, MakeRCGrid and ApplyRC.
//
// One config file (config/rge_pip.cfg) drives all three, so the table, the RC
// grid and the weights cannot disagree on beam energy, particle or cuts.

#ifndef RGESELECTION_HXX
#define RGESELECTION_HXX

#include "ConfigFile.h"
#include "haprad_constants.h"

#include "TChain.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

namespace rge {

// The structure-function table: cells in (Q2, nu, z, pt2), phi binned for the fit.
const char* const kTableAxis[5] = {"Q2", "nu", "z", "pt2", "phi"};
const char* const kTableBinKey[5] = {"bins_Q2", "bins_nu", "bins_z", "bins_pt2", "bins_phi"};

// The RC grid: NODES in (Q2, nu, z, pt, phi). p_t, not pt2, because the RC
// factor's phi modulation grows like p_t; phi in degrees, 0..180 (the
// unpolarised cross section is symmetric under phi -> 360 - phi).
const char* const kRCAxis[5] = {"Q2", "nu", "z", "pt", "phi"};
const char* const kRCNodeKey[5] = {"rc_nodes_Q2", "rc_nodes_nu", "rc_nodes_z", "rc_nodes_pt", "rc_nodes_phi"};

struct Config {
  std::string text;  // verbatim, stored in every output
  Double_t E, Etol;
  Int_t pid;
  Double_t Q2min, W2min, ymax, zmin, zmax, pt2max, Mxmin, vzmin, vzmax;
  Int_t n[5];
  Double_t lo[5], hi[5];
  Int_t fitMinPhiBins;
  Double_t fitMinEvents;
  std::string weightBranch;  // empty: every event counts 1
  Double_t lumi;             // events per nb of cross section; 1 = arbitrary units
  // RC grid and weights (MakeRCGrid, ApplyRC). hasRCGrid is false when the
  // rc_nodes_* keys are absent.
  bool hasRCGrid;
  Int_t rcN[5];
  Double_t rcLo[5], rcHi[5];
  Double_t targetNAZ;    // Z/A: the share of the exclusive tail (MAID is a proton model)
  bool rcExclusiveTail;  // include the exclusive tail in the weight
};

inline bool LoadConfig(const char* path, Config& c) {
  std::ifstream in(path);
  if (!in) {
    fprintf(stderr, "cannot read config %s\n", path);
    return false;
  }
  std::stringstream ss;
  ss << in.rdbuf();
  c.text = ss.str();

  try {
    ConfigFile cf(path);
    c.E = cf.read<Double_t>("beam_energy");
    c.Etol = cf.read<Double_t>("beam_energy_tolerance", 0.005);
    c.pid = cf.read<Int_t>("pid", 211);
    c.Q2min = cf.read<Double_t>("Q2_min", 1.);
    c.W2min = cf.read<Double_t>("W2_min", 4.);
    c.ymax = cf.read<Double_t>("y_max", 0.85);
    c.zmin = cf.read<Double_t>("z_min", 0.);
    c.zmax = cf.read<Double_t>("z_max", 1.);
    c.pt2max = cf.read<Double_t>("pt2_max", 1e9);
    c.Mxmin = cf.read<Double_t>("Mx_min");
    c.vzmin = cf.read<Double_t>("vz_min", -1e9);
    c.vzmax = cf.read<Double_t>("vz_max", 1e9);
    c.fitMinPhiBins = cf.read<Int_t>("fit_min_phi_bins", 4);
    c.fitMinEvents = cf.read<Double_t>("fit_min_events", 50.);
    c.weightBranch = cf.read<std::string>("weight_branch", "");
    c.lumi = cf.read<Double_t>("luminosity", 1.);
    if (!(c.lumi > 0.)) {
      fprintf(stderr, "luminosity must be positive\n");
      return false;
    }
    for (Int_t d = 0; d < 5; ++d) {
      std::istringstream b(cf.read<std::string>(kTableBinKey[d]));
      if (!(b >> c.n[d] >> c.lo[d] >> c.hi[d]) || c.n[d] < 1 || !(c.hi[d] > c.lo[d])) {
        fprintf(stderr, "%s must be '<nbins> <low> <high>'\n", kTableBinKey[d]);
        return false;
      }
    }
    c.hasRCGrid = cf.keyExists(kRCNodeKey[0]);
    for (Int_t d = 0; c.hasRCGrid && d < 5; ++d) {
      std::istringstream b(cf.read<std::string>(kRCNodeKey[d]));
      if (!(b >> c.rcN[d] >> c.rcLo[d] >> c.rcHi[d]) || c.rcN[d] < 2 || !(c.rcHi[d] > c.rcLo[d])) {
        fprintf(stderr, "%s must be '<nodes, at least 2> <first> <last>'\n", kRCNodeKey[d]);
        return false;
      }
    }
    c.targetNAZ = cf.read<Double_t>("target_naz", 0.5);
    const std::string ex = cf.read<std::string>("rc_exclusive_tail", "no");
    if (ex != "yes" && ex != "no") {
      fprintf(stderr, "rc_exclusive_tail must be yes or no\n");
      return false;
    }
    c.rcExclusiveTail = (ex == "yes");
  } catch (ConfigFile::key_not_found& e) {
    fprintf(stderr, "config %s: required key '%s' is missing\n", path, e.key.c_str());
    return false;
  }
  return true;
}

// M_x^2 = (P + q - p_h)^2 from quantities in the ntuple. The hadron is taken to
// be a pion; the result matched HAPRAD's own missing-mass rejections exactly
// on run 020026 (PLAN.md 0.1).
inline Double_t Mx2(Double_t W2, Double_t nu, Double_t Q2, Double_t zh, Double_t thetaPQ) {
  const Double_t m = kMassPion, M = kMassProton;
  const Double_t Eh = zh * nu;
  const Double_t ph = std::sqrt(std::max(0., Eh * Eh - m * m));
  return W2 + m * m - 2. * (M + nu) * Eh + 2. * std::sqrt(nu * nu + Q2) * ph * std::cos(thetaPQ);
}

// One entry of the RG-E `DT` ntuple: the branches the selection reads.
struct Event {
  Float_t pid, Eb, Q2, nu, yb, W2, zh, pt2, phi, theta, vz;
};

// Enables and binds the branches of Event; everything else is switched off.
inline bool BindEvent(TChain& chain, Event& e) {
  struct B {
    const char* name;
    Float_t* v;
  } br[] = {{"pid", &e.pid}, {"E_beam", &e.Eb}, {"Q2", &e.Q2},       {"nu", &e.nu},
            {"y_bjorken", &e.yb}, {"W2", &e.W2}, {"z_h", &e.zh},    {"p_T2", &e.pt2},
            {"phi_PQ", &e.phi},  {"theta_PQ", &e.theta}, {"vz", &e.vz}};
  chain.SetBranchStatus("*", 0);
  for (const B& b : br) {
    chain.SetBranchStatus(b.name, 1);
    if (chain.SetBranchAddress(b.name, b.v) < 0) {
      fprintf(stderr, "input lacks branch %s\n", b.name);
      return false;
    }
  }
  return true;
}

// The selection, stage by stage, in order.
const char* const kStage[] = {"read", "pid", "beam energy", "Q2", "W2", "y", "z", "pt2", "vertex", "missing mass"};
const Int_t kNStages = sizeof(kStage) / sizeof(kStage[0]);

// How many stages the event survives: 1 (read only) .. kNStages (selected).
inline Int_t Stages(const Config& c, const Event& e) {
  Int_t s = 1;
  if (Int_t(e.pid) != c.pid) return s;
  ++s;
  if (std::fabs(e.Eb - c.E) > c.Etol) return s;  // one beam energy per table
  ++s;
  if (!(e.Q2 > c.Q2min)) return s;
  ++s;
  if (!(e.W2 > c.W2min)) return s;
  ++s;
  if (!(e.yb < c.ymax)) return s;
  ++s;
  if (!(e.zh > c.zmin && e.zh < c.zmax)) return s;
  ++s;
  if (!(e.pt2 >= 0. && e.pt2 < c.pt2max)) return s;
  ++s;
  if (!(e.vz > c.vzmin && e.vz < c.vzmax)) return s;
  ++s;
  if (!(Mx2(e.W2, e.nu, e.Q2, e.zh, e.theta) > c.Mxmin * c.Mxmin)) return s;
  return ++s;
}

}  // namespace rge

#endif
