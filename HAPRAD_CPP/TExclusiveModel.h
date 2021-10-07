#ifndef TEXCLUSIVEMODEL_H
#define TEXCLUSIVEMODEL_H

const Int_t nq = 18;
const Int_t nw = 47;
const Int_t nt = 61;

void ExclusiveModel(Double_t q2m, Double_t wm, Double_t csthcm, Double_t &st,
                    Double_t &sl, Double_t &stt, Double_t &stl, Double_t &stlp);

#endif
