#ifndef THAPRADUTILS_H
#define THAPRADUTILS_H

#include "TROOT.h"

namespace HapradUtils {

    Double_t vacpol(const Double_t Q2);
    Double_t fspen(const Double_t x);
    Double_t fspens(const Double_t x);
    Double_t dfint(Int_t narg, double *arg, Int_t *nent, Double_t *ent, Double_t *table);
}

#endif
