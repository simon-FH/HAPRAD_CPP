#include "TStructFunctionArray.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "TSemiInclusiveModel.h"
#include "haprad_constants.h"
#include "square_power.h"

//______________________________________________________________________________
//
// The function calculates deep inelastic (ita = 1), elastic (ita = 2),
// quasielastic (ita = 3) structure functions in kinematical point (tau,R).
//
//     R = S_x - tt
//     tau = (t - Q2) / R
//
// where
//
//     tt = t + mf2 - M^2
//
// mf2 is invarint mass of final hadrons

TStructFunctionArray::TStructFunctionArray(const TRadCor *rc) {
  fKin = rc->GetKinematicalVariables();
  fInv = rc->GetLorentzInvariants();
  fHadKin = rc->GetHadronKinematics();
  fModel = rc->GetSemiInclusiveModel();
}

TStructFunctionArray::~TStructFunctionArray() {
  // Do nothing
}

void TStructFunctionArray::Evaluate(Double_t tau, Double_t mu, Double_t R) {
  using namespace TMath;

  // strf() in ihaprad.f opens by zeroing sfm(1..8); do the same. Degenerate
  // shifted kinematics must contribute nothing, and TPODINL constructs a fresh
  // TStructFunctionArray for every single integrand evaluation, so anything
  // left unassigned here is read back as uninitialised stack.
  for (Int_t i = 0; i < 4; ++i) fArray[i] = 0.;

  const Double_t &M = kMassProton;
  const Double_t &m_h = kMassDetectedHadron;

  Double_t tldQ2 = fInv->Q2() + R * tau;
  Double_t tldSx = fInv->Sx() - R;

  Double_t tldLq = SQ(tldSx) + 4. * SQ(M) * tldQ2;
  Double_t tldNu = tldSx / (2 * M);

  Double_t tldT = fKin->T() - R * (tau - mu);

  Double_t tldPx2 = fHadKin->Px2() - R * (1. + tau - mu);
  Double_t PhQ = (SQ(m_h) - tldQ2 - tldT) / 2.;

  Double_t tldX = tldQ2 / (2 * M * tldNu);
  Double_t tldZ = fHadKin->Eh() / tldNu;

  Double_t tld_sq = Sqrt(tldQ2 + SQ(tldNu));
  Double_t tldPl = (fHadKin->Eh() * tldNu - PhQ) / tld_sq;
  Double_t tldPt2 = SQ(fHadKin->Ph()) - SQ(tldPl);

  Double_t epsNu = kEpsMachine * Sqrt(SQ(fHadKin->Nu()) + SQ(fHadKin->Nu() - tldNu) + SQ(tldNu));

  Double_t epsT = kEpsMachine * Sqrt(SQ(fKin->T()) + SQ((tau - mu) * R) + SQ(R * tau) + SQ(R * mu));

  Double_t epsPhq = Sqrt(SQ(2 * SQ(m_h) * kEpsMachine) + SQ(tldQ2 * kEpsMachine) + SQ(epsT)) / 2;

  Double_t epsQ2 = Sqrt(SQ(tldQ2 * kEpsMachine) + SQ(2 * SQ(tldNu) * epsNu)) / (2. * Sqrt(tldQ2 + SQ(tldNu)));

  Double_t epsPl = Sqrt(SQ(fHadKin->Eh() * kEpsMachine * tldNu / tld_sq) + SQ(fHadKin->Eh() * epsNu / tld_sq) + SQ(epsPhq / tld_sq) +
                        SQ((fHadKin->Eh() * tldNu - PhQ) / SQ(tld_sq) * epsQ2));

  Double_t epsPt2 = 2. * Sqrt(SQ(SQ(fHadKin->Ph())) * SQ(kEpsMachine) + SQ(tldPl) * SQ(epsPl));

  if (tldPt2 >= 0 || (SQ(tldPt2) - SQ(epsPt2)) <= 0) {
    Double_t a = fInv->S() / (2 * M) * (fInv->S() / (2 * M) - fHadKin->Nu()) * tldQ2 / fInv->Q2();

    Double_t tldE = 0.5 * (tldNu + Sqrt(SQ(tldNu) + 4 * a));
    Double_t tldY = tldNu / tldE;

    // Measured azimuthal amplitudes at the (shifted) kinematics of this
    // evaluation. These were three uninitialised doubles, filled only by an
    // external script that sed-substituted constants into this file.
    Double_t A = 0., Ac = 0., Acc = 0.;
    if (fModel) fModel->Evaluate(tldQ2, tldX, tldY, tldZ, tldPt2, tldPx2, tldPl, A, Ac, Acc);

    // Relation between data yield and Hizs
    Double_t rlt = 0.14;  // Constant separating the long. & trans. cross-sections
    Double_t coef = (M * tldQ2 / 4. / tldZ);
    Double_t tldXg = 2. * M * (fKin->E() - tldNu);
    Double_t tldSp = fInv->S() + tldXg;
    // HAY QUE TENER CUIDADO CON LOS ZEROS
    Double_t tldLambda = fInv->S() * tldXg * tldQ2 - SQ(M * tldQ2);  // electron part is ignored: m_e*lamndaQ
    Double_t LambdaRatio = tldLq / tldLambda;
    Double_t RelH1H2 = (1. + 4. * SQ(M * tldX) / tldQ2) / 2. / tldX / (1 + rlt);
    Double_t N = SQ(tldQ2) * Sqrt(tldQ2 + SQ(tldY)) / SQ(tldY);
    Double_t H2z = N * A * coef / (RelH1H2 * SQ(tldQ2) + (fInv->S() * tldXg - SQ(M) * tldQ2) * tldX);
    Double_t H1z = H2z * RelH1H2;
    Double_t H3z = N * Ac * Sqrt(LambdaRatio) * coef / tldSp / Sqrt(tldPt2);
    Double_t H4z = N * Acc * LambdaRatio * coef * tldX / tldPt2;

    // Including kinematic shift due to photon emission
    Double_t aa = tldSx * (tldZ - 2 * M * tldPl / Sqrt(tldLq)) / (2. * SQ(M));

    Double_t h1 = tldZ * (2. * tldQ2 * tldX * H1z - tldPt2 * H4z) / (M * tldQ2 * tldX);

    Double_t h2 =
        2 * tldZ * (2 * tldX * (tldX * H2z + aa * H3z) * tldLq + H4z * (SQ(aa) * tldLq - 2 * tldPt2 * tldQ2)) / M / tldX / tldQ2 / tldLq;

    Double_t h3 = 2 * H4z * tldZ / M / tldQ2 / tldX;

    Double_t h4 = -2 * tldZ * (tldX * H3z + aa * H4z) / M / tldQ2 / tldX;

    // H3z and H4z carry 1/sqrt(tldPt2) and 1/tldPt2 because this model inverts
    // the measured azimuthal amplitudes A_c and A_cc, which enter the Born
    // cross section multiplied by b ~ p_t and b^2 ~ p_t^2. The FORTRAN model
    // has no such inversion (semi_inclusive_model.f builds H3/H4 directly from
    // PDFs x FFs), so it stays regular here.
    //
    // The guard above admits tldPt2 == 0 (it tests >= 0), and there H4z is
    // infinite while the h_i combine it as tldPt2 * H4z, i.e. 0 * inf = NaN.
    // That NaN propagated all the way out through the tau and R integrals and
    // poisoned the whole inelastic tail. It is reachable in practice: at
    // phi_h = 180 deg the shifted p_t passes exactly through zero somewhere in
    // the integration domain.
    //
    // Treat any degenerate evaluation as no contribution, which is what strf()
    // does for tldPt2 < 0 and is the only defensible reading -- these are
    // isolated points of the integration domain, not a region.
    //
    // Note this is a robustness fix, not a physics one: the 1/p_t^2 blow-up is
    // an artifact of feeding this inversion azimuthal amplitudes that do not
    // carry their proper p_t dependence (see hurdle H1). A model in which
    // A_c ~ p_t and A_cc ~ p_t^2 keeps H3z and H4z finite.
    if (Finite(h1) && Finite(h2) && Finite(h3) && Finite(h4)) {
      fArray[0] = h1;
      fArray[1] = h2;
      fArray[2] = h3;
      fArray[3] = h4;
    }
  }
}
