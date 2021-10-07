#include "TDelta.h"
#include "TRadCor.h"
#include "THapradConfig.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "THapradUtils.h"
#include "haprad_constants.h"
#include "square_power.h"


TDelta::TDelta(const TRadCor* rc)
  : fVR(0), fInf(0), fVac(0)
{
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();
}



TDelta::~TDelta()
{
    // Do nothing
}



void TDelta::Evaluate(void)
{
    Double_t m;

    switch(fConfig->PolarizationType()) {
        case 1:
            m = kMassElectron;
            break;
        case 2:
            m = kMassMuon;
            break;
        default:
            m = kMassElectron;
    }

    Double_t S_ = fInv->X() + fInv->Q2() - fHadKin->V2();
    Double_t X_ = fInv->S() - fInv->Q2() - fHadKin->V1();

    Double_t l_m  = TMath::Log(fInv->Q2() / SQ(m));
    Double_t Li_2 = HapradUtils::fspen(1 - fHadKin->Px2() * fInv->Q2() / (S_ * X_));

    fVR  = 1.5 * l_m - 2 - 0.5 * SQ(TMath::Log(X_ / S_)) + Li_2 - SQ(kPi) / 6;
    fInf = (l_m - 1) * TMath::Log(SQ(fHadKin->Px2() - kMassC2) / S_ / X_);
    fVac = HapradUtils::vacpol(fInv->Q2());

    fVR  *= kAlpha / kPi;
    fInf *= kAlpha / kPi;
    fVac *= kAlpha / kPi;
}
