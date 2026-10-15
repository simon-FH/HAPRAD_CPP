#ifndef TSEMIINCLUSIVEMODEL_H
#define TSEMIINCLUSIVEMODEL_H

#include "TROOT.h"
#include "THn.h"

#include <memory>
#include <string>
#include <vector>

//______________________________________________________________________________
//
// Table of measured azimuthal amplitudes, A + Ac cos(phi) + Acc cos(2 phi), in
// 4-D kinematic cells. TStructFunctionArray inverts these into the structure
// functions H1..H4. Nothing is fitted here: the table is built outside HAPRAD
// (one fit of the phi distribution per cell) and only looked up.
//
// FILE FORMAT
//
//   THnD "A", "Ac", "Acc"   the amplitudes, one value per cell. All three must
//                           have identical binning. Each axis NAME says which
//                           variable it is:
//                             Q2
//                             x    (aliases: xB, Xb)
//                             nu
//                             z    (aliases: zh, Zh)
//                             pt   (alias: Pt)
//                             pt2  (aliases: Pt2, p_T2)
//                           Any four, in any order. The binning lives in the
//                           file, so the table and its reader cannot disagree.
//
//   THnD "fitted"           optional. Non-zero where the cell holds a fit.
//                           Absent means every cell is fitted.
//
//   TNamed "interpolation"  optional, "linear" (default) or "log".
//                             log: interpolate log(A), and the ratios Ac/A and
//                                  Acc/A, linearly. A carries a steep,
//                                  cross-section-like kinematic dependence
//                                  (roughly 1/Q^4) that linear interpolation
//                                  between cell centres follows badly; its
//                                  logarithm and the modulation ratios are
//                                  smooth. Requires A > 0 in every fitted cell
//                                  used; cells with A <= 0 are treated as
//                                  unfitted.
//
//   TNamed "pt_scaling"     optional, "raw" (default) or "reduced".
//                             raw:     the THnDs hold Ac and Acc themselves.
//                             reduced: they hold Ac/p_t and Acc/p_t^2, and the
//                                      lookup multiplies back by the p_t of the
//                                      point being evaluated.
//                           In the Born cross section Ac ~ p_t and Acc ~ p_t^2.
//                           A raw table returns a non-zero constant across the
//                           lowest p_t cell, and TStructFunctionArray divides
//                           by p_t and p_t^2, which diverges where the shifted
//                           p_t reaches zero (at phi = 180 deg). A reduced
//                           table vanishes with the right power and avoids it.
//
// INTERPOLATION
//
//   Values are interpolated multilinearly in all four variables between cell
//   centres -- the stored quantity, i.e. Ac/p_t and Acc/p_t^2 for a reduced
//   table, before the p_t factors are applied. A nearest-cell lookup makes the
//   integrand piecewise constant, and the non-adaptive Gauss-Kronrod rule used
//   for the inner R integral (TRV2LN) cannot converge across those steps (with
//   a toy table, thousands of tolerance failures per call, against none for a
//   smooth one). SetInterpolation(false) selects the nearest-cell lookup, for
//   comparison only.
//
// OUT-OF-RANGE AND EMPTY CELLS
//
//   Beyond the outermost cell centres the value is held constant (clamped). A
//   point outside the table's edges is counted as out of range. Corners that fall
//   in unfitted cells are dropped and the remaining weights renormalised. If no
//   corner is fitted, the second neighbours are used: the same interpolation in
//   a filled copy of the table, made at Load(), in which every unfitted cell
//   holds the average of the fitted cells next to it (its 3x3x3x3 block,
//   weighted by 1 / squared distance in cells). Seen from the point, that
//   reaches the second ring of cells; it is continuous and costs no more than
//   an ordinary lookup (counted as NSecondNeighbour). Only if the filled copy
//   has no usable corner either are the amplitudes zero and the point counted
//   as empty.
//   The counters exist because the tail integrals sample shifted kinematics
//   well outside the measured region; they show how much of an RC factor rests
//   on clamped, filled or missing table values.

class TSemiInclusiveModel {
 public:
  TSemiInclusiveModel();
  ~TSemiInclusiveModel();

  // Read a table. Returns false, with a message on stderr, if the file is
  // unreadable or malformed; the previous table (if any) is then discarded.
  Bool_t Load(const char* path);
  Bool_t IsLoaded() const { return fA != nullptr; }
  const std::string& GetPath() const { return fPath; }
  Bool_t IsPtReduced() const { return fPtReduced; }
  Bool_t IsLogInterpolation() const { return fLogA; }

  void SetInterpolation(Bool_t on) { fInterpolate = on; }
  Bool_t GetInterpolation() const { return fInterpolate; }
  // Second-neighbour fallback for points with no fitted corner (see above). On
  // by default; off leaves such points empty, for comparison.
  void SetSecondNeighbours(Bool_t on) { fSecondNeighbours = on; }
  Bool_t GetSecondNeighbours() const { return fSecondNeighbours; }

  // Interpolate A, Ac and Acc at one kinematic point (the variables named by
  // the table's axes are taken from q2, X, Z, pt2 and mx2, with nu derived from
  // q2 and X). Same argument list as the original free function; Y and pl are
  // unused. Updates the diagnostic counters below.
  void Evaluate(Double_t q2, Double_t X, Double_t Y, Double_t Z, Double_t pt2, Double_t mx2, Double_t pl, Double_t& A,
                Double_t& Ac, Double_t& Acc) const;

  // Diagnostics, accumulated over Evaluate() calls since the last reset.
  void ResetCounters() const;
  Long64_t NLookups() const { return fNLookups; }
  Long64_t NOutOfRange() const { return fNOutOfRange; }
  Long64_t NEmptyCell() const { return fNEmptyCell; }
  Long64_t NSecondNeighbour() const { return fNSecondNeighbour; }
  Long64_t NUnphysical() const { return fNUnphysical; }

 private:
  enum EVariable { kQ2, kX, kNu, kZ, kPt, kPt2 };

  void Clear();
  // Build fFillA/fFillAc/fFillAcc/fFillOk: the filled copy of the table used
  // for second neighbours. Called at the end of Load().
  void BuildSecondNeighbours();

  std::string fPath;
  std::unique_ptr<THnD> fA;
  std::unique_ptr<THnD> fAc;
  std::unique_ptr<THnD> fAcc;
  std::unique_ptr<THnD> fFitted;
  std::unique_ptr<THnD> fFillA;  // filled copies for second neighbours
  std::unique_ptr<THnD> fFillAc;
  std::unique_ptr<THnD> fFillAcc;
  std::vector<char> fFillOk;  // per global bin: usable in the filled copy
  EVariable fVar[4];
  Bool_t fPtReduced;
  Bool_t fInterpolate;
  Bool_t fLogA;
  Bool_t fSecondNeighbours;

  mutable Long64_t fNLookups;
  mutable Long64_t fNOutOfRange;
  mutable Long64_t fNEmptyCell;
  mutable Long64_t fNSecondNeighbour;
  mutable Long64_t fNUnphysical;
};

#endif
