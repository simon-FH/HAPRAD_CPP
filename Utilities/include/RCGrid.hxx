// RCGrid.hxx -- the radiative-correction grid: HAPRAD evaluated at the nodes of
// a regular grid in (Q2, nu, z, pt, phi), and interpolated to any event.
//
// Written by MakeRCGrid, read by ApplyRC and MakeRCGrid's error check.
//
// Each node holds sigma_Born and three ratios to it, as TRadCor computes them:
//
//     r_vr = sig_obs / sigma_Born   (Born x virtual and soft-photon factor)
//     r_in = tai[0]  / sigma_Born   (inelastic radiative tail)
//     r_ex = tai[1]  / sigma_Born   (exclusive radiative tail, MAID)
//
// The RC factor sigma_obs / sigma_Born is r_vr + r_in + c r_ex, with c = 0
// (TRadCor's GetFactor1) or Z/A (GetFactor3). Keeping the parts separate lets
// that choice be made when the weights are applied -- it has to be, because the
// exclusive tail is meaningful only for a table in absolute units (PLAN.md
// Phase 1).
//
// Storage: one THnD per quantity with ONE BIN PER NODE, centred on it (bin i+1
// holds node i), so that grids computed in chunks merge with `hadd`. The
// "computed" histogram counts how often each node was evaluated; after a
// complete merge it is 1 everywhere, which Load() checks.

#ifndef RCGRID_HXX
#define RCGRID_HXX

#include "RGESelection.hxx"

#include "TFile.h"
#include "THn.h"
#include "TNamed.h"

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace rge {

const char* const kRCComp[3] = {"r_vr", "r_in", "r_ex"};

struct RCNodes {
  Int_t n[5];
  Double_t lo[5], hi[5];

  explicit RCNodes(const Config& c) {
    for (Int_t d = 0; d < 5; ++d) {
      n[d] = c.rcN[d];
      lo[d] = c.rcLo[d];
      hi[d] = c.rcHi[d];
    }
  }
  RCNodes() {}
  Double_t Step(Int_t d) const { return (hi[d] - lo[d]) / (n[d] - 1); }
  Double_t Node(Int_t d, Int_t i) const { return lo[d] + i * Step(d); }
  Long64_t Size() const {
    Long64_t s = 1;
    for (Int_t d = 0; d < 5; ++d) s *= n[d];
    return s;
  }
  // Linear index <-> per-axis node indices (0-based), last axis fastest.
  Long64_t Index(const Int_t i[5]) const {
    Long64_t k = 0;
    for (Int_t d = 0; d < 5; ++d) k = k * n[d] + i[d];
    return k;
  }
  void Decode(Long64_t k, Int_t i[5]) const {
    for (Int_t d = 4; d >= 0; --d) {
      i[d] = k % n[d];
      k /= n[d];
    }
  }
  THnD* MakeHist(const char* name) const {
    Double_t a[5], b[5];
    for (Int_t d = 0; d < 5; ++d) {
      a[d] = lo[d] - Step(d) / 2.;
      b[d] = hi[d] + Step(d) / 2.;
    }
    THnD* h = new THnD(name, name, 5, n, a, b);
    for (Int_t d = 0; d < 5; ++d) h->GetAxis(d)->SetName(kRCAxis[d]);
    return h;
  }
};

class RCGrid {
 public:
  enum Status { kInside = 0, kPartial = 1, kOutside = 2, kNoValidCorner = 3 };

  // Reads a merged grid. Returns false, with the reason in `why`, if the file
  // is not a complete grid.
  bool Load(const char* path, std::string& why) {
    std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
    if (!f || f->IsZombie()) {
      why = std::string("cannot open ") + path;
      return false;
    }
    std::unique_ptr<THnD> computed(f->Get<THnD>("computed")), valid(f->Get<THnD>("valid"));
    std::unique_ptr<THnD> comp[3];
    for (Int_t k = 0; k < 3; ++k) comp[k].reset(f->Get<THnD>(kRCComp[k]));
    if (!computed || !valid || !comp[0] || !comp[1] || !comp[2] || computed->GetNdimensions() != 5) {
      why = std::string(path) + " is not an RC grid";
      return false;
    }
    for (Int_t d = 0; d < 5; ++d) {
      const TAxis* ax = computed->GetAxis(d);
      fG.n[d] = ax->GetNbins();
      const Double_t step = ax->GetBinWidth(1);
      fG.lo[d] = ax->GetXmin() + step / 2.;
      fG.hi[d] = ax->GetXmax() - step / 2.;
    }
    const Long64_t size = fG.Size();
    for (Int_t k = 0; k < 3; ++k) fR[k].assign(size, 0.f);
    fValid.assign(size, 0);
    Long64_t missing = 0, twice = 0;
    Int_t i[5], bin[5];
    for (Long64_t k = 0; k < size; ++k) {
      fG.Decode(k, i);
      for (Int_t d = 0; d < 5; ++d) bin[d] = i[d] + 1;
      const Double_t done = computed->GetBinContent(bin);
      if (done < 0.5) ++missing;
      if (done > 1.5) ++twice;
      if (valid->GetBinContent(bin) > 0.5) {
        fValid[k] = 1;
        for (Int_t c = 0; c < 3; ++c) fR[c][k] = comp[c]->GetBinContent(bin);
      }
    }
    if (missing || twice) {
      why = Form("%s is incomplete: %lld nodes never computed, %lld computed more than once "
                 "(a chunk missing from the merge, or merged twice)",
                 path, missing, twice);
      return false;
    }
    for (const char* key : {"table", "table_luminosity", "beam_energy"})
      if (TNamed* t = f->Get<TNamed>(key)) fInfo[key] = t->GetTitle();
    return true;
  }

  // Multilinear interpolation of the three ratios at v = (Q2, nu, z, pt,
  // phi[deg, 0..180]). Corners where HAPRAD gave no result are left out and the
  // rest renormalised; `cover` is the interpolation weight that had a result
  // (1 inside the valid region). Outside the grid nothing is extrapolated.
  Status Interpolate(const Double_t v[5], Double_t r[3], Double_t& cover) const {
    r[0] = r[1] = r[2] = 0.;
    cover = 0.;
    Int_t i0[5];
    Double_t f[5];
    for (Int_t d = 0; d < 5; ++d) {
      const Double_t t = (v[d] - fG.lo[d]) / fG.Step(d);
      if (!(t >= -1e-9 && t <= fG.n[d] - 1 + 1e-9)) return kOutside;
      i0[d] = std::min(std::max(Int_t(std::floor(t)), 0), fG.n[d] - 2);
      f[d] = std::min(std::max(t - i0[d], 0.), 1.);
    }
    Double_t sum[3] = {0., 0., 0.};
    Int_t i[5];
    for (Int_t corner = 0; corner < 32; ++corner) {
      Double_t w = 1.;
      for (Int_t d = 0; d < 5; ++d) {
        const Int_t up = (corner >> d) & 1;
        i[d] = i0[d] + up;
        w *= up ? f[d] : 1. - f[d];
      }
      if (w == 0.) continue;
      const Long64_t k = fG.Index(i);
      if (!fValid[k]) continue;
      cover += w;
      for (Int_t c = 0; c < 3; ++c) sum[c] += w * fR[c][k];
    }
    if (cover <= 0.) return kNoValidCorner;
    for (Int_t c = 0; c < 3; ++c) r[c] = sum[c] / cover;
    return cover > 1. - 1e-9 ? kInside : kPartial;
  }

  const RCNodes& Nodes() const { return fG; }
  std::string Info(const std::string& key) const {
    auto it = fInfo.find(key);
    return it == fInfo.end() ? "" : it->second;
  }

 private:
  RCNodes fG;
  std::vector<float> fR[3];
  std::vector<char> fValid;
  std::map<std::string, std::string> fInfo;
};

// The event's position in the grid's variables. phi_PQ is in radians in RG-E;
// the grid holds phi in degrees on 0..180, using phi -> 360 - phi symmetry.
inline void GridPoint(const Event& e, Double_t v[5]) {
  v[0] = e.Q2;
  v[1] = e.nu;
  v[2] = e.zh;
  v[3] = std::sqrt(std::max(0.f, e.pt2));
  v[4] = std::fabs(std::remainder(Double_t(e.phi), 2. * M_PI)) * 180. / M_PI;
}

// Whether the exclusive tail may be used with a table: only if the table is in
// absolute units. Tables from MakePhiTable record "luminosity"; one without it
// set says "1 (arbitrary units ...". Model tables carry no such key and are
// absolute (they come from HAPRAD 2.0's cross section).
inline bool TableIsAbsolute(const std::string& tableLuminosity) {
  return tableLuminosity.rfind("1 (arbitrary", 0) != 0;
}

}  // namespace rge

#endif
