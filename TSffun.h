#ifndef TSFFUN_H
#define TSFFUN_H

#include "TROOT.h"

class TRadCor;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;

class TSffun {
public:
    TSffun(const TRadCor *rc);
    ~TSffun();

    void        Evaluate(Double_t Q2, Double_t w2, Double_t t);

    Double_t&   operator() (Int_t i) { return fArray[i]; };
    Double_t    operator() (Int_t i) const { return fArray[i]; };

private:
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;
    const THadronKinematics*        fHadKin;

    Double_t    fArray[4];
};


#endif
