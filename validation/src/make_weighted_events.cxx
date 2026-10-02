// make_weighted_events -- toy pions weighted by HAPRAD 2.0's Born cross section,
// for the normalisation closure of MakePhiTable (PLAN.md Phase 1).
//
// Events are uniform in the BINNING variables (Q2, nu, z, pt2, phi). Weighting
// each by
//
//     w = sigma_Born(x, y, z, pt2, phi) * |d(x, y)/d(Q2, nu)|,   |...| = 1/(2 M E nu)
//
// makes the weighted sample a measurement of that cross section, binned the way
// real data will be. MakePhiTable must then turn it into a table that agrees
// with the one make_model_table builds directly from the same model -- up to
// one overall constant, and the usual bin-averaging of a steep function.
//
// Two steps, because the FORTRAN runs as a separate program:
//
//   make_weighted_events points <seed> <n> <E> <box>              > points.txt
//   born_harmonics < points.txt                                   > harm.txt
//   make_weighted_events build  <seed> <n> <E> <box> harm.txt out.root
//
// <box> = "Q2lo:Q2hi,nulo:nuhi,zlo:zhi,pt2lo:pt2hi". The same seed regenerates
// the same events in `build`. Weights are scaled to a mean of one; `build`
// prints the luminosity that scaling amounts to,
//
//     L = (draws, physical or not) x (weight scale) / (box volume, phi in rad),
//
// in events per unit of HAPRAD's sigma (nb). Given it, MakePhiTable must
// reproduce the model's A exactly, not just up to a constant.
//
// Each event also carries its own harmonics (B0, Bc, Bcc; zero where HAPRAD
// rejects the point), and the file the draw count and box volume, so that the
// closure can compute the true cell AVERAGES and separate bin-centring from
// errors in the producer. MakePhiTable reads none of these.

#include "haprad_constants.h"

#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TParameter.h"
#include "TRandom3.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Box {
  double lo[4], hi[4];
};

bool ParseBox(const char* s, Box& b) {
  return std::sscanf(s, "%lf:%lf,%lf:%lf,%lf:%lf,%lf:%lf", &b.lo[0], &b.hi[0], &b.lo[1], &b.hi[1], &b.lo[2], &b.hi[2],
                     &b.lo[3], &b.hi[3]) == 8;
}

struct Event {
  double Q2, nu, z, pt2, phi;
};

// The n-th accepted event of the sequence. Rejected draws (unphysical) are
// skipped identically in both steps, so ids line up. Returns the number of
// draws, which fixes the event density in the box.
template <class F>
long Generate(unsigned seed, long n, double E, const Box& b, F f) {
  TRandom3 rng(seed);
  const double M = kMassProton, m = kMassPion;
  long draws = 0;
  for (long id = 0; id < n;) {
    ++draws;
    Event e;
    e.Q2 = rng.Uniform(b.lo[0], b.hi[0]);
    e.nu = rng.Uniform(b.lo[1], b.hi[1]);
    e.z = rng.Uniform(b.lo[2], b.hi[2]);
    e.pt2 = rng.Uniform(b.lo[3], b.hi[3]);
    e.phi = rng.Uniform(-TMath::Pi(), TMath::Pi());
    const double x = e.Q2 / (2. * M * e.nu), Eh = e.z * e.nu;
    if (x >= 1. || e.nu >= E || Eh <= m || e.pt2 >= Eh * Eh - m * m) continue;
    f(id, e);
    ++id;
  }
  return draws;
}

}  // namespace

int main(int argc, char** argv) {
  Box b;
  if (argc < 6 || !ParseBox(argv[5], b)) {
    fprintf(stderr,
            "usage: make_weighted_events points <seed> <n> <E> <box>                  > points.txt\n"
            "       make_weighted_events build  <seed> <n> <E> <box> harm.txt out.root\n"
            "  <box> = Q2lo:Q2hi,nulo:nuhi,zlo:zhi,pt2lo:pt2hi\n");
    return 2;
  }
  const std::string mode = argv[1];
  const unsigned seed = std::atoi(argv[2]);
  const long n = std::atol(argv[3]);
  const double E = std::atof(argv[4]);
  const double M = kMassProton;

  if (mode == "points") {
    Generate(seed, n, E, b, [&](long id, const Event& e) {
      printf("%ld 0 0 0 %.10g %.12g %.12g %.12g %.12g\n", id, E, e.Q2 / (2. * M * e.nu), e.Q2, e.z, std::sqrt(e.pt2));
    });
    return 0;
  }

  if (mode == "build" && argc >= 8) {
    // Harmonics by event id; PDFLIB's banner and rejected points are skipped.
    std::vector<double> B0(n, 0.), Bc(n, 0.), Bcc(n, 0.);
    std::vector<char> ok(n, 0);
    std::ifstream in(argv[6]);
    std::string line;
    while (std::getline(in, line)) {
      std::istringstream ls(line);
      long id, j1, j2, j3;
      int k;
      std::string s0, s1, s2;
      if (!(ls >> id >> j1 >> j2 >> j3 >> k >> s0 >> s1 >> s2) || id < 0 || id >= n || k != 1) continue;
      for (std::string* s : {&s0, &s1, &s2})
        for (char& c : *s)
          if (c == 'D' || c == 'd') c = 'e';
      B0[id] = std::atof(s0.c_str());
      Bc[id] = std::atof(s1.c_str());
      Bcc[id] = std::atof(s2.c_str());
      ok[id] = 1;
    }

    std::vector<Event> ev(n);
    std::vector<double> w(n, 0.);
    double sumW = 0.;
    long nOk = 0;
    const long draws = Generate(seed, n, E, b, [&](long id, const Event& e) {
      ev[id] = e;
      if (!ok[id]) return;  // HAPRAD rejects the point: no cross section
      const double sigma = B0[id] + Bc[id] * std::cos(e.phi) + Bcc[id] * std::cos(2. * e.phi);
      w[id] = sigma / (2. * M * E * e.nu);
      sumW += w[id];
      ++nOk;
    });
    if (nOk == 0) {
      fprintf(stderr, "no event has a cross section -- wrong harmonics file?\n");
      return 1;
    }
    const double scale = double(n) / sumW;
    double volume = 2. * TMath::Pi();
    for (int d = 0; d < 4; ++d) volume *= b.hi[d] - b.lo[d];
    const double lumi = draws * scale / volume;

    TFile out(argv[7], "RECREATE");
    TNtuple dt("DT", "toy pions weighted by HAPRAD 2.0 Born",
               "pid:E_beam:Q2:nu:x_bjorken:y_bjorken:W2:z_h:p_T2:phi_PQ:theta_PQ:vz:w:B0:Bc:Bcc");
    const double m = kMassPion;
    for (long id = 0; id < n; ++id) {
      const Event& e = ev[id];
      const double Eh = e.z * e.nu, ph2 = Eh * Eh - m * m;
      const double theta = std::acos(std::sqrt(ph2 - e.pt2) / std::sqrt(ph2));
      const Float_t row[16] = {211.f, Float_t(E), Float_t(e.Q2), Float_t(e.nu), Float_t(e.Q2 / (2. * M * e.nu)),
                               Float_t(e.nu / E), Float_t(M * M + 2. * M * e.nu - e.Q2), Float_t(e.z), Float_t(e.pt2),
                               Float_t(e.phi), Float_t(theta), 0.f, Float_t(w[id] * scale),
                               Float_t(B0[id]), Float_t(Bc[id]), Float_t(Bcc[id])};
      dt.Fill(row);
    }
    dt.Write();
    TParameter<double>("draws", draws).Write();
    TParameter<double>("box_volume", volume).Write();  // phi in rad
    TParameter<double>("luminosity", lumi).Write();
    out.Close();
    printf("wrote %ld events (%ld with a cross section) to %s\n", n, nOk, argv[7]);
    printf("luminosity %.10g\n", lumi);
    return 0;
  }
  fprintf(stderr, "unknown mode '%s'\n", mode.c_str());
  return 2;
}
