#include "TLorentzInvariants.h"
#include "HapradErrors.h"
#include "THapradConfig.h"
#include "TKinematicalVariables.h"
#include "THadronKinematics.h"
#include "THapradException.h"
#include "haprad_constants.h"
#include "square_power.h"

#include <iostream>


TLorentzInvariants::TLorentzInvariants(THapradConfig* config,
                                       TKinematicalVariables* kin)
 : fConfig(config), fKin(kin), fHadKin(0),
   fS(0), fX(0), fSx(0), fSp(0), fQ2(0)
{
    // Do nothing
}



TLorentzInvariants::~TLorentzInvariants()
{
    // Do nothing
}



void TLorentzInvariants::Evaluate(void)
{
    // Evaluate all variables.

    SetS();
    SetQ2();
    SetX();
    SetSx();
    SetSp();
    SetW2();
    SetLambdas();
}



void TLorentzInvariants::SetS()
{
    fS = 2 * kMassProton * fKin->E();
}



void TLorentzInvariants::SetX()
{
    fX  = fS * (1 - fKin->Y());
}



void TLorentzInvariants::SetSx()
{
    fSx = fS - fX;
}



void TLorentzInvariants::SetSp()
{
    fSp = fS + fX;
}



void TLorentzInvariants::SetQ2()
{
    Double_t Y;

    if (fKin->Y() >= 0) {
        Y = - fKin->Y();
        fQ2 = fS * fKin->X() * Y;
    } else {
        fQ2 = - fKin->Y();
        Y = fQ2 / (fS * fKin->X());
    }

    fKin->SetY(Y);

    Double_t y_max = 1 / (1 + SQ(kMassProton) * fKin->X() / fS);
    Double_t y_min = (kMassC2 - SQ(kMassProton)) / (fS * (1 - fKin->X()));

    if (fKin->Y() > y_max || fKin->Y() < y_min ||
                fKin->X() > 1 || fKin->X() < 0) {
        std::cout << "    y: " << fKin->Y()
                  << "    x: " << fKin->X()
                  << std::endl;
        throw TKinematicException();
    }
}



void TLorentzInvariants::SetW2()
{
    fW2 = fS - fX - fQ2 + SQ(kMassProton);
}



void TLorentzInvariants::SetLambdas()
{
    Double_t M = kMassProton;
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

    fLambdaS = SQ(fS)  - 4 * SQ(m) * SQ(M);
    fLambdaX = SQ(fX)  - 4 * SQ(m) * SQ(M);
    fLambdaQ = SQ(fSx) + 4 * SQ(M) * fQ2;
    fLambdaM = SQ(fQ2) + 4 * SQ(m) * fQ2;

    if (fLambdaS < 0) HAPRAD_WARN_MSG("Invariants", "lambda_s < 0");
    if (fLambdaX < 0) HAPRAD_WARN_MSG("Invariants", "lambda_x < 0");
    if (fLambdaQ < 0) HAPRAD_WARN_MSG("Invariants", "lambda_q < 0");
    if (fLambdaM < 0) HAPRAD_WARN_MSG("Invariants", "lambda_m < 0");

    fSqrtLs = TMath::Sqrt(TMath::Max(0.,fLambdaS));
    fSqrtLq = TMath::Sqrt(TMath::Max(0.,fLambdaQ));
    fSqrtLx = TMath::Sqrt(TMath::Max(0.,fLambdaX));
}
