#include "TRV2TR.h"
#include "TRadCor.h"
#include "THapradConfig.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "TSffun.h"
#include "TThetaMatrix.h"
#include "haprad_constants.h"
#include "square_power.h"
#include "TMath.h"
#include <iostream>


TRV2TR::TRV2TR(const TRadCor* rc)
{
    fRC     = rc;
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();

    M   = kMassProton;
    M2  = SQ(kMassProton);
    mh2 = SQ(kMassDetectedHadron);
    mu2 = SQ(kMassUndetectedHadron);

    tau_max = (fInv->Sx() + fInv->SqrtLq()) / (2. * SQ(M));
    tau_min = - fInv->Q2() / SQ(M) / tau_max;
}



TRV2TR::~TRV2TR()
{

}



ROOT::Math::IBaseFunctionMultiDim* TRV2TR::Clone() const
{
    return 0;
}



unsigned int TRV2TR::NDim() const
{
    return 2;
}



double TRV2TR::DoEval(const double *x) const
{
    double tau   = x[0];
    double phi_k = x[1];

    double vv = (1 - fKin->Z()) * fInv->Sx() + fKin->T() + M2 - mu2;
    double mu = (fHadKin->Eh() -
                 fHadKin->Pl() * (fInv->Sx() - tau * 2 * M2) / fInv->SqrtLq() -
                 2 * M2 * fHadKin->Pt() *
                        TMath::Sqrt((tau - tau_min) * (tau_max - tau)) *
                        TMath::Cos(fKin->PhiH() - phi_k) / fInv->SqrtLq()) / M;

    double fwiw  = 1. + tau - mu;
    double R     = vv / fwiw;
    double tldQ2 = fInv->Q2() + R * tau;
    double tldQ4 = SQ(tldQ2);
    double tldW2 = fInv->W2() - R * (1. + tau);
    double tldT  = fKin->T() - R * (tau - mu);

    TSffun tldH(fRC);
    tldH.Evaluate(tldQ2,tldW2,tldT);

    TThetaMatrix theta(fRC);
    theta.Evaluate(tau,mu,2,phi_k);

    double podinlz = 0.;
    for (int i = 0; i < 4; ++i) {
        double factor = tldH(i) / tldQ4;
        double pres = theta(i,0) / R + theta(i,1) + theta(i,2) * R;
        podinlz += pres * factor;
    }

    return podinlz / fwiw;
}
