#ifndef TSTRUCTFUNCTIONARRAY_H
#define TSTRUCTFUNCTIONARRAY_H

#include "TROOT.h"


class TRadCor;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;


class TStructFunctionArray {
public:
    TStructFunctionArray(const TRadCor* rc);
    ~TStructFunctionArray();

    void        Evaluate(Double_t tau, Double_t mu, Double_t R);

    Double_t&   operator() (Int_t i) { return fArray[i]; };
    Double_t    operator() (Int_t i) const { return fArray[i]; };

private:
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;
    const THadronKinematics*        fHadKin;

    Double_t    fArray[4];
};

#endif
