// phi_convention -- check that RG-E's phi_PQ is the phi_h HAPRAD expects.
//
// A convention mismatch here would be silent: shifting phi by pi flips the sign
// of cos(phi), and therefore of A_c, and nothing downstream complains.
//
// The test uses the invariant V1 = 2 k1.p_h. It needs only the beam and the
// hadron, both available in every ntuple row:
//
//   from four-vectors:  V1 = 2 E (E_h - p_z)       beam along +z, m_e neglected
//   from HAPRAD:        THadronKinematics::V1() evaluated at
//                       (E, x_B, Q2, z_h, sqrt(p_T2), phi_PQ)
//
// HAPRAD writes V1 = 2 (a1 + b cos phi_h). Because V1 depends on phi only through
// cos(phi), it cannot tell phi from -phi -- which does not matter for an
// unpolarised analysis -- but it does separate phi from phi + pi.
//
// HAPRAD's p_t input assumes the hadron goes forward in the virtual-photon
// frame (p_l = +sqrt(p_h^2 - p_t^2)), so backward hadrons are counted but left
// out of the comparison.
//
// usage: phi_convention "<file glob>" [max_entries]

#include "THapradConfig.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "THapradException.h"
#include "haprad_constants.h"

#include "TChain.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// HAPRAD's V1 for one kinematic point, or NaN if the point is rejected.
double HapradV1(double E, double x, double Q2, double z, double pt, double phi_deg) {
  THapradConfig cfg;
  TKinematicalVariables kin(x, -Q2, z, pt, phi_deg / kRadianDeg, E);
  TLorentzInvariants inv(&cfg, &kin);
  THadronKinematics had(&cfg, &kin, &inv);
  try {
    inv.Evaluate();
    had.Evaluate();
  } catch (TKinematicException&) {
    return NAN;
  }
  return had.V1();
}

double Quantile(std::vector<double> v, double q) {
  if (v.empty()) return NAN;
  std::sort(v.begin(), v.end());
  return v[std::min(v.size() - 1, size_t(q * v.size()))];
}

double FractionBelow(const std::vector<double>& v, double cut) {
  if (v.empty()) return NAN;
  return double(std::count_if(v.begin(), v.end(), [cut](double d) { return d < cut; })) / v.size();
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: phi_convention \"<file glob>\" [max_entries]\n");
    return 2;
  }
  const Long64_t maxEntries = (argc > 2) ? atoll(argv[2]) : -1;

  TChain chain("DT");
  chain.Add(argv[1]);

  Float_t pid, Ebeam, px, py, pz, mass, Q2, xb, yb, W2, zh, pt2, phiPQ, thetaPQ;
  chain.SetBranchAddress("pid", &pid);
  chain.SetBranchAddress("E_beam", &Ebeam);
  chain.SetBranchAddress("px", &px);
  chain.SetBranchAddress("py", &py);
  chain.SetBranchAddress("pz", &pz);
  chain.SetBranchAddress("mass", &mass);
  chain.SetBranchAddress("Q2", &Q2);
  chain.SetBranchAddress("x_bjorken", &xb);
  chain.SetBranchAddress("y_bjorken", &yb);
  chain.SetBranchAddress("W2", &W2);
  chain.SetBranchAddress("z_h", &zh);
  chain.SetBranchAddress("p_T2", &pt2);
  chain.SetBranchAddress("phi_PQ", &phiPQ);
  chain.SetBranchAddress("theta_PQ", &thetaPQ);

  // THadronKinematics prints a diagnostic line for every evaluation; mute it.
  std::ostringstream sink;
  std::streambuf* coutBuf = std::cout.rdbuf();

  std::vector<double> diffAsIs, diffShifted, phiLever;
  Long64_t nSel = 0, nBackward = 0, nRejected = 0;
  Long64_t nRejPx2 = 0, nRejY = 0, nRejT = 0, nRejEh = 0, nRejOther = 0;

  const Long64_t nEntries = chain.GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    chain.GetEntry(i);
    if (pid != 211) continue;
    if (!(Q2 > 1 && W2 > 4 && yb < 0.85 && zh > 0.15 && zh < 1.0 && pt2 > 0 && pt2 < 1.5)) continue;
    ++nSel;
    if (thetaPQ > TMath::PiOver2()) {
      ++nBackward;
      continue;
    }

    const double Eh = std::sqrt(double(mass) * mass + double(px) * px + double(py) * py + double(pz) * pz);
    const double V1lab = 2. * Ebeam * (Eh - pz);
    const double phiDeg = phiPQ * kRadianDeg;  // RG-E stores radians

    std::cout.rdbuf(sink.rdbuf());
    const double v0 = HapradV1(Ebeam, xb, Q2, zh, std::sqrt(pt2), phiDeg);
    const double v1 = HapradV1(Ebeam, xb, Q2, zh, std::sqrt(pt2), phiDeg + 180.);
    std::cout.rdbuf(coutBuf);

    if (!std::isfinite(v0) || !std::isfinite(v1)) {
      // HAPRAD prints which check failed just before throwing; classify by it.
      const std::string why = sink.str();
      if (why.find("px2:") != std::string::npos)
        ++nRejPx2;
      else if (why.find("y:") != std::string::npos)
        ++nRejY;
      else if (why.find("t:") != std::string::npos)
        ++nRejT;
      else if (why.find("E_h:") != std::string::npos || why.find("p_h:") != std::string::npos)
        ++nRejEh;
      else
        ++nRejOther;
      ++nRejected;
      sink.str("");
      continue;
    }
    sink.str("");
    diffAsIs.push_back(std::fabs(v0 - V1lab) / V1lab);
    diffShifted.push_back(std::fabs(v1 - V1lab) / V1lab);
    phiLever.push_back(std::fabs(v0 - v1) / V1lab);

    if (maxEntries > 0 && Long64_t(diffAsIs.size()) >= maxEntries) break;
  }

  printf("pi+ passing DIS cuts        : %lld\n", nSel);
  printf("  backward (theta_PQ > 90)  : %lld  (%.2f%%)  -- excluded, see header\n", nBackward,
         nSel ? 100. * nBackward / nSel : 0.);
  printf("  rejected by HAPRAD        : %lld\n", nRejected);
  printf("      missing mass p_x^2 out of range : %lld\n", nRejPx2);
  printf("      y out of range                  : %lld\n", nRejY);
  printf("      t out of range                  : %lld\n", nRejT);
  printf("      E_h < m_h or p_h < p_t          : %lld\n", nRejEh);
  printf("      other                           : %lld\n", nRejOther);
  printf("  compared                  : %zu\n\n", diffAsIs.size());

  printf("Sensitivity: |V1(phi) - V1(phi+pi)| / V1  -- how hard the test can bite\n");
  printf("  median %.3e    5%% quantile %.3e\n\n", Quantile(phiLever, 0.5), Quantile(phiLever, 0.05));

  printf("%-22s %12s %12s %12s %14s\n", "hypothesis", "median", "95%", "max", "frac < 1e-4");
  printf("%-22s %12.3e %12.3e %12.3e %14.4f\n", "phi_h = phi_PQ", Quantile(diffAsIs, 0.5), Quantile(diffAsIs, 0.95),
         Quantile(diffAsIs, 1.0), FractionBelow(diffAsIs, 1e-4));
  printf("%-22s %12.3e %12.3e %12.3e %14.4f\n", "phi_h = phi_PQ + pi", Quantile(diffShifted, 0.5),
         Quantile(diffShifted, 0.95), Quantile(diffShifted, 1.0), FractionBelow(diffShifted, 1e-4));
  return 0;
}
