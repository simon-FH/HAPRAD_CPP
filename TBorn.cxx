#include "TBorn.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "haprad_constants.h"
#include "square_power.h"

TBorn::TBorn(const TRadCor* rc)
  : fH(rc)
{
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();

    const Double_t& M   = kMassProton;
    const Double_t& m_h = kMassDetectedHadron;

    fThetaB[0] = fInv->Q2();
    fThetaB[1] = (fInv->S() * fInv->X() - SQ(M) * fInv->Q2()) / 2.;
    fThetaB[2] = (fHadKin->V1() * fHadKin->V2() - SQ(m_h) * fInv->Q2()) / 2;
    fThetaB[3] = (fHadKin->V2() * fInv->S() + fHadKin->V1() * fInv->X() -
                            fKin->Z() * fInv->Q2() * fInv->Sx()) / 2;
}



TBorn::~TBorn()
{

}



Double_t TBorn::Evaluate(void)
{
    fH.Evaluate(0.,0.,0.);

    Double_t sum = 0.;
    for (Int_t i = 0; i < 4; ++i) {
        sum = sum + fThetaB[i] * fH(i);
    }


    return 2 * sum / SQ(fInv->Q2()) ;
}
