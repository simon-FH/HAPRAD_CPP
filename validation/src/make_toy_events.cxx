// make_toy_events -- toy pions in the RG-E ntuple format (tree "DT"), with a
// KNOWN azimuthal modulation, for closure tests of MakePhiTable (PLAN.md 2.2).
//
// Kinematics are uniform in (Q2, nu, z, pt2) inside physical limits; phi is
// drawn from
//
//     1 + rc cos(phi) + rcc cos(2 phi),   rc = -0.15 (pt/0.5),  rcc = 0.05 (pt/0.5)^2
//
// so every cell has its own true modulation, with the p_t scaling the Born cross
// section requires. Only the branches MakePhiTable reads are written, plus the
// derived ones (x_bjorken, ...) for completeness.
//
// usage: make_toy_events <out.root> <n_events> [seed]

#include "haprad_constants.h"

#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TRandom3.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

double TrueRc(double pt) { return -0.15 * pt / 0.5; }
double TrueRcc(double pt) { return 0.05 * (pt / 0.5) * (pt / 0.5); }

int main(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: make_toy_events <out.root> <n_events> [seed]\n");
    return 2;
  }
  const Long64_t nWanted = atoll(argv[2]);
  TRandom3 rng(argc > 3 ? atoi(argv[3]) : 12345);

  const double E = 10.5473, M = kMassProton, m = kMassPion;
  TFile out(argv[1], "RECREATE");
  TNtuple dt("DT", "toy pions with known azimuthal modulation",
             "pid:E_beam:Q2:nu:x_bjorken:y_bjorken:W2:z_h:p_T2:phi_PQ:theta_PQ:vz");

  Long64_t n = 0;
  while (n < nWanted) {
    const double Q2 = rng.Uniform(1.0, 6.0);
    const double nu = rng.Uniform(3.0, 9.5);
    const double z = rng.Uniform(0.15, 0.95);
    const double pt2 = rng.Uniform(0.0, 1.2);
    const double Eh = z * nu;
    if (Eh <= m) continue;
    const double ph2 = Eh * Eh - m * m;
    if (pt2 >= ph2) continue;  // forward hadron: p_l = +sqrt(ph^2 - pt^2)
    const double pt = std::sqrt(pt2);
    const double theta = std::acos(std::sqrt(ph2 - pt2) / std::sqrt(ph2));

    // phi by accept-reject on the known modulation
    const double rc = TrueRc(pt), rcc = TrueRcc(pt);
    const double fmax = 1. + std::fabs(rc) + std::fabs(rcc);
    double phi;
    do {
      phi = rng.Uniform(-TMath::Pi(), TMath::Pi());
    } while (rng.Uniform(0., fmax) > 1. + rc * std::cos(phi) + rcc * std::cos(2. * phi));

    const double x = Q2 / (2. * M * nu), y = nu / E;
    const double W2 = M * M + 2. * M * nu - Q2;
    const Float_t row[12] = {211.f, Float_t(E), Float_t(Q2), Float_t(nu), Float_t(x), Float_t(y), Float_t(W2),
                             Float_t(z), Float_t(pt2), Float_t(phi), Float_t(theta), 0.f};
    dt.Fill(row);
    ++n;
  }
  dt.Write();
  out.Close();
  printf("wrote %lld toy pions to %s\n", n, argv[1]);
  return 0;
}
