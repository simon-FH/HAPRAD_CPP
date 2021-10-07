#ifndef TSEMIINCLUSIVEMODEL_H
#define TSEMIINCLUSIVEMODEL_H

#include "TROOT.h"

namespace HapradUtils {

    void SemiInclusiveModel(Double_t q2, Double_t X,
                            Double_t Y, Double_t Z,
                            Double_t pt2, Double_t mx2,
                            Double_t pl, Double_t& A,
                            Double_t& Ac, Double_t& Acc);

}

#endif
