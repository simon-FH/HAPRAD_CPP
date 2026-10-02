#include "TSemiInclusiveModel.h"
#include "haprad_constants.h"
#include "square_power.h"

#include "TAxis.h"
#include "TFile.h"
#include "THn.h"
#include "TMath.h"
#include "TNamed.h"

#include <cmath>
#include <iostream>

namespace {

Bool_t SameAxis(const TAxis* a, const TAxis* b) {
  if (a->GetNbins() != b->GetNbins()) return false;
  for (Int_t i = 1; i <= a->GetNbins() + 1; ++i)
    if (a->GetBinLowEdge(i) != b->GetBinLowEdge(i)) return false;
  return true;
}

Bool_t SameBinning(const THnD* a, const THnD* b) {
  if (a->GetNdimensions() != b->GetNdimensions()) return false;
  for (Int_t d = 0; d < a->GetNdimensions(); ++d)
    if (!SameAxis(a->GetAxis(d), b->GetAxis(d))) return false;
  return true;
}

}  // namespace

TSemiInclusiveModel::TSemiInclusiveModel()
    : fPtReduced(false), fInterpolate(true), fLogA(false), fNLookups(0), fNOutOfRange(0), fNEmptyCell(0), fNUnphysical(0) {
  for (Int_t i = 0; i < 4; ++i) fVar[i] = kQ2;
}

TSemiInclusiveModel::~TSemiInclusiveModel() {}

void TSemiInclusiveModel::Clear() {
  fPath.clear();
  fA.reset();
  fAc.reset();
  fAcc.reset();
  fFitted.reset();
  fPtReduced = false;
  fLogA = false;
}

void TSemiInclusiveModel::ResetCounters() const {
  fNLookups = 0;
  fNOutOfRange = 0;
  fNEmptyCell = 0;
  fNUnphysical = 0;
}

Bool_t TSemiInclusiveModel::Load(const char* path) {
  Clear();
  ResetCounters();

  const char* where = "TSemiInclusiveModel::Load";
  std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
  if (!file || file->IsZombie()) {
    std::cerr << where << ": cannot open '" << path << "'" << std::endl;
    return false;
  }

  std::unique_ptr<THnD> a(file->Get<THnD>("A"));
  std::unique_ptr<THnD> ac(file->Get<THnD>("Ac"));
  std::unique_ptr<THnD> acc(file->Get<THnD>("Acc"));
  if (!a || !ac || !acc) {
    std::cerr << where << ": '" << path << "' must hold THnD objects named A, Ac and Acc" << std::endl;
    return false;
  }
  if (a->GetNdimensions() != 4) {
    std::cerr << where << ": tables must be 4-dimensional, got " << a->GetNdimensions() << std::endl;
    return false;
  }
  if (!SameBinning(a.get(), ac.get()) || !SameBinning(a.get(), acc.get())) {
    std::cerr << where << ": A, Ac and Acc do not share the same binning" << std::endl;
    return false;
  }

  std::unique_ptr<THnD> fitted(file->Get<THnD>("fitted"));
  if (fitted && !SameBinning(a.get(), fitted.get())) {
    std::cerr << where << ": 'fitted' does not share the binning of A" << std::endl;
    return false;
  }

  // Each axis name says which kinematic variable it is.
  for (Int_t d = 0; d < 4; ++d) {
    const TString n = a->GetAxis(d)->GetName();
    if (n == "Q2")
      fVar[d] = kQ2;
    else if (n == "x" || n == "xB" || n == "Xb")
      fVar[d] = kX;
    else if (n == "nu")
      fVar[d] = kNu;
    else if (n == "z" || n == "zh" || n == "Zh")
      fVar[d] = kZ;
    else if (n == "pt" || n == "Pt")
      fVar[d] = kPt;
    else if (n == "pt2" || n == "Pt2" || n == "p_T2")
      fVar[d] = kPt2;
    else {
      std::cerr << where << ": axis " << d << " is named '" << n << "', which is not a known variable"
                << " (Q2, x, nu, z, pt, pt2)" << std::endl;
      return false;
    }
  }

  Bool_t reduced = false;
  if (TNamed* s = file->Get<TNamed>("pt_scaling")) {
    const TString v = s->GetTitle();
    delete s;
    if (v == "reduced")
      reduced = true;
    else if (v != "raw") {
      std::cerr << where << ": pt_scaling must be 'raw' or 'reduced', got '" << v << "'" << std::endl;
      return false;
    }
  }

  Bool_t logA = false;
  if (TNamed* s = file->Get<TNamed>("interpolation")) {
    const TString v = s->GetTitle();
    delete s;
    if (v == "log")
      logA = true;
    else if (v != "linear") {
      std::cerr << where << ": interpolation must be 'linear' or 'log', got '" << v << "'" << std::endl;
      return false;
    }
  }

  fA = std::move(a);
  fAc = std::move(ac);
  fAcc = std::move(acc);
  fFitted = std::move(fitted);
  fPtReduced = reduced;
  fLogA = logA;
  fPath = path;
  return true;
}

void TSemiInclusiveModel::Evaluate(Double_t q2, Double_t X, Double_t /*Y*/, Double_t Z, Double_t pt2, Double_t mx2,
                                   Double_t /*pl*/, Double_t& A, Double_t& Ac, Double_t& Acc) const {
  // Every early return leaves the amplitudes at zero, never at whatever the
  // caller's stack held -- the same convention as strf() in ihaprad.f.
  A = 0.;
  Ac = 0.;
  Acc = 0.;
  if (!fA) return;

  ++fNLookups;

  if (X <= 0 || X > 1 || Z < 0 || Z > 1 || mx2 < SQ(kMassProton + kMassPion)) {
    ++fNUnphysical;
    return;
  }

  const Double_t pt2pos = TMath::Max(0., pt2);

  // For each axis: the two neighbouring cells to interpolate between, and the
  // weight of the upper one.
  Int_t lo[4], hi[4];
  Double_t w[4];
  Bool_t outside = false;
  for (Int_t d = 0; d < 4; ++d) {
    Double_t v = 0.;
    switch (fVar[d]) {
      case kQ2: v = q2; break;
      case kX: v = X; break;
      case kNu: v = q2 / (2. * kMassProton * X); break;
      case kZ: v = Z; break;
      case kPt: v = std::sqrt(pt2pos); break;
      case kPt2: v = pt2pos; break;
    }
    const TAxis* ax = fA->GetAxis(d);
    const Int_t n = ax->GetNbins();
    if (v < ax->GetXmin() || v >= ax->GetXmax()) outside = true;

    Int_t b = ax->FindFixBin(v);
    if (b < 1) b = 1;
    if (b > n) b = n;

    if (!fInterpolate) {
      lo[d] = hi[d] = b;
      w[d] = 0.;
      continue;
    }
    // Neighbours by cell centre; constant beyond the outermost centres.
    Int_t b1 = (v < ax->GetBinCenter(b)) ? b - 1 : b;
    Int_t b2 = b1 + 1;
    if (b1 < 1) {
      lo[d] = hi[d] = 1;
      w[d] = 0.;
    } else if (b2 > n) {
      lo[d] = hi[d] = n;
      w[d] = 0.;
    } else {
      const Double_t c1 = ax->GetBinCenter(b1), c2 = ax->GetBinCenter(b2);
      lo[d] = b1;
      hi[d] = b2;
      w[d] = (v - c1) / (c2 - c1);
    }
  }
  if (outside) ++fNOutOfRange;

  // Sum over the 2^4 corners, dropping unfitted ones and renormalising.
  Double_t sumW = 0., sA = 0., sAc = 0., sAcc = 0.;
  Int_t idx[4];
  for (Int_t corner = 0; corner < 16; ++corner) {
    Double_t wc = 1.;
    for (Int_t d = 0; d < 4; ++d) {
      const Bool_t up = (corner >> d) & 1;
      if (up && hi[d] == lo[d]) {
        wc = 0.;
        break;
      }
      idx[d] = up ? hi[d] : lo[d];
      wc *= up ? w[d] : 1. - w[d];
    }
    if (wc == 0.) continue;
    const Long64_t bin = fA->GetBin(idx);
    if (fFitted && fFitted->GetBinContent(bin) == 0.) continue;
    const Double_t a = fA->GetBinContent(bin);
    if (fLogA) {
      if (!(a > 0.)) continue;
      sumW += wc;
      sA += wc * std::log(a);
      sAc += wc * fAc->GetBinContent(bin) / a;
      sAcc += wc * fAcc->GetBinContent(bin) / a;
    } else {
      sumW += wc;
      sA += wc * a;
      sAc += wc * fAc->GetBinContent(bin);
      sAcc += wc * fAcc->GetBinContent(bin);
    }
  }
  if (sumW <= 0.) {
    ++fNEmptyCell;
    return;
  }

  if (fLogA) {
    A = std::exp(sA / sumW);
    Ac = A * sAc / sumW;
    Acc = A * sAcc / sumW;
  } else {
    A = sA / sumW;
    Ac = sAc / sumW;
    Acc = sAcc / sumW;
  }
  if (fPtReduced) {
    Ac *= std::sqrt(pt2pos);
    Acc *= pt2pos;
  }
}
