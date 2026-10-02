// make_toy_table -- write a structure-function table with KNOWN amplitudes, in
// the format TSemiInclusiveModel reads (see TSemiInclusiveModel.h).
//
// For testing the consumer before any real table exists: the amplitudes are
// smooth analytic functions, so what HAPRAD gets back can be checked against
// what went in.
//
// Variables follow PLAN.md 0.3 -- (Q2, nu, z, pt2), the ones RG-E's acc_corr
// bins in -- so the same acceptance map can serve both.
//
// usage: make_toy_table <out.root> [raw|reduced]

#include "TFile.h"
#include "THn.h"
#include "TMath.h"
#include "TNamed.h"

#include <cmath>
#include <cstdio>
#include <cstring>

namespace {

// Known amplitudes. A falls with z and has a Gaussian p_t^2 dependence whose
// width grows with z; the cos(phi) and cos(2 phi) modulations are -15% and +5%
// at p_t = 0.5 GeV and scale like p_t and p_t^2, as the Born cross section
// requires.
double ToyA(double /*Q2*/, double /*nu*/, double z, double pt2) {
  const double width = 0.2 + 0.25 * z * z;
  return 1000. * std::pow(1. - 0.9 * z, 1.5) * std::exp(-pt2 / width) / width;
}
double ToyAcOverPt(double Q2, double nu, double z, double pt2) { return -0.15 / 0.5 * ToyA(Q2, nu, z, pt2); }
double ToyAccOverPt2(double Q2, double nu, double z, double pt2) { return 0.05 / 0.25 * ToyA(Q2, nu, z, pt2); }

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: make_toy_table <out.root> [raw|reduced]\n");
    return 2;
  }
  const bool reduced = (argc < 3) || std::strcmp(argv[2], "raw") != 0;

  const char* names[4] = {"Q2", "nu", "z", "pt2"};
  Int_t nbins[4] = {9, 8, 10, 6};
  Double_t lo[4] = {1., 2., 0., 0.};
  Double_t hi[4] = {10., 10., 1., 1.5};

  THnD A("A", "A", 4, nbins, lo, hi);
  THnD Ac("Ac", "Ac", 4, nbins, lo, hi);
  THnD Acc("Acc", "Acc", 4, nbins, lo, hi);
  THnD fitted("fitted", "fitted", 4, nbins, lo, hi);
  for (THnD* h : {&A, &Ac, &Acc, &fitted})
    for (Int_t d = 0; d < 4; ++d) h->GetAxis(d)->SetName(names[d]);

  Int_t idx[4];
  for (idx[0] = 1; idx[0] <= nbins[0]; ++idx[0])
    for (idx[1] = 1; idx[1] <= nbins[1]; ++idx[1])
      for (idx[2] = 1; idx[2] <= nbins[2]; ++idx[2])
        for (idx[3] = 1; idx[3] <= nbins[3]; ++idx[3]) {
          double v[4];
          for (Int_t d = 0; d < 4; ++d) v[d] = A.GetAxis(d)->GetBinCenter(idx[d]);
          const double pt2 = v[3];
          A.SetBinContent(idx, ToyA(v[0], v[1], v[2], pt2));
          if (reduced) {
            Ac.SetBinContent(idx, ToyAcOverPt(v[0], v[1], v[2], pt2));
            Acc.SetBinContent(idx, ToyAccOverPt2(v[0], v[1], v[2], pt2));
          } else {
            Ac.SetBinContent(idx, ToyAcOverPt(v[0], v[1], v[2], pt2) * std::sqrt(pt2));
            Acc.SetBinContent(idx, ToyAccOverPt2(v[0], v[1], v[2], pt2) * pt2);
          }
          fitted.SetBinContent(idx, 1.);
        }

  TFile out(argv[1], "RECREATE");
  A.Write();
  Ac.Write();
  Acc.Write();
  fitted.Write();
  TNamed("pt_scaling", reduced ? "reduced" : "raw").Write();
  out.Close();
  printf("wrote %s: %d cells, axes Q2 nu z pt2, pt_scaling = %s\n", argv[1],
         nbins[0] * nbins[1] * nbins[2] * nbins[3], reduced ? "reduced" : "raw");
  return 0;
}
