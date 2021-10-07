#include "TSffun.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "TMath.h"
#include "TExclusiveModel.h"
#include "haprad_constants.h"
#include "square_power.h"
#include <iostream>


TSffun::TSffun(const TRadCor* rc)
{
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();
}

TSffun::~TSffun()
{
    // Do Nothing, default destructor
}

void TSffun::Evaluate(Double_t Q2, Double_t w2, Double_t t)
{
    using namespace TMath;

    Double_t st, sl, stt, slt, sltp, sfm10;
    Double_t M_p = kMassProton;
    Double_t m = kMassDetectedHadron;
    Double_t M_u = kMassUndetectedHadron;

    Double_t Sx = w2 + Q2 - SQ(M_p);
    Double_t SqrtW2 = Sqrt(w2);
    Double_t LambdaQ = SQ(Sx) + 4 * SQ(M_p) * Q2;
    Double_t Sx_t = Sx + t + SQ(M_p) - SQ(M_u);
    Double_t tq = t + Q2 - SQ(m);

    Double_t sffun_cmp = SQ(w2 - SQ(M_u) - SQ(m)) - 4 * SQ(M_u) * SQ(m);

    if (sffun_cmp < 0)
        std::cout << "sffun: SqrtLw=NaN " << sffun_cmp << std::endl;
    

    Double_t SqrtLw = Sqrt(Max(0., sffun_cmp));

    Double_t ssffun_cmp2 = Q2 * SQ(Sx_t) - Sx_t * Sx * tq - SQ(M_p) * SQ(tq) - SQ(m) * LambdaQ;

    if (ssffun_cmp2 < 0)
        std::cout << "ssffun: qll=NaN " << ssffun_cmp2 << std::endl;
    

    Double_t SqrtLl = Sqrt(Max(0., ssffun_cmp2));
    Double_t cspion = (2 * tq * w2 + (Sx - 2 * Q2) * (w2 + SQ(m) - SQ(M_u))) / SqrtLw / Sqrt(LambdaQ);

//  Exclusive peak model (cross sections sigma_L,T,LT... from MAID2003)

    ExclusiveModel(Q2, SqrtW2, cspion, st, sl, stt, slt, sltp);

// Structure functions

    if (LambdaQ > 0 && SqrtLl > 0 && SqrtLw > 0) {
        Double_t sfm20;
        Double_t sfm2tl;
        Double_t sfm4tl;
        Double_t sfm5tl;
        Double_t sfm4tt;
        Double_t sfm3tt;
        Double_t sfm2tt;
        Double_t Coetr;
        sfm10 = st - stt;
        sfm20 = 4. * (st + sl) * Q2 / LambdaQ;
        sfm2tl = 2. * slt * Sqrt(Q2) * (-Sx * tq + 2. * Q2 * Sx_t) / (LambdaQ * SqrtLl);
        sfm4tl = -slt * Sqrt(Q2) / SqrtLl;
        sfm4tt = -2. * stt * (-Sx * tq + 2. * Q2 * Sx_t) / SQ(SqrtLl);
        sfm3tt = 2. * stt * LambdaQ / SQ(SqrtLl);
        sfm2tt = 2. * stt * (SQ(-Sx * tq + 2. * Q2 * Sx_t) - 2. * Q2 * SQ(SqrtLl)) / (LambdaQ * SQ(SqrtLl));
        sfm5tl = -sltp * Sqrt(Q2) / SqrtLl;

        Coetr = 16. * kPi * (w2 - SQ(M_p)) * w2 / (kAlpha * SqrtLw) / kBarn * 1000.;
        fArray[0] = Coetr * sfm10;
        fArray[1] = Coetr * (sfm20 + sfm2tl + sfm2tt);
        fArray[2] = Coetr * sfm3tt;
        fArray[3] = Coetr * (sfm4tl + sfm4tt);

    } else {
        fArray[0] = 0;
        fArray[1] = 0;
        fArray[2] = 0;
        fArray[3] = 0;
    }
}
