#include "THapradUtils.h"
#include "Partons.h"

#include "haprad_constants.h"
#include "square_power.h"


namespace HapradUtils {


Double_t vacpol(const Double_t Q2)
{
    Double_t leptonMass[3] = { 0.26110E-6,
                               0.111637E-1,
                               3.18301 };

    Double_t suml = 0;
    for (Int_t i = 0; i < 3; ++i) {
        Double_t a2    = 2. * leptonMass[i];
        Double_t sqlmi = TMath::Sqrt(SQ(Q2) + 2. * a2 * Q2);
        Double_t allmi = TMath::Log((sqlmi + Q2) / (sqlmi - Q2)) / sqlmi;

        suml = suml + 2. * (Q2 + a2) * allmi / 3. - 10. / 9. +
                    4. * a2 * (1. - a2 * allmi) / 3. / Q2;
    }

    Double_t a, b, c;

    if (Q2 < 1) {
        a = -1.345E-9;
        b = -2.302E-3;
        c =  4.091;
    } else if (Q2 < 64) {
        a = -1.512E-3;
        b = -2.822E-3;
        c =  1.218;
    } else {
        a = -1.1344E-3;
        b = -3.0680E-3;
        c =  9.9992E-1;
    }

    Double_t sumh;
    sumh = - (a + b * TMath::Log(1. + c * Q2)) * 2. * kPi / kAlpha;

    return suml + sumh;
}



Double_t fspen(const Double_t x)
{
    const Double_t f1 = 1.644934;
    Double_t result;

    if (x < -1) {
        Double_t logprod;
        logprod = TMath::Log(1 - x) * TMath::Log(SQ(x) / (1 - x));
        result = - 0.5 * logprod - f1 + fspens(1 / (1 - x));
    } else if (x < 0) {
        Double_t log2;
        log2 = TMath::Log(1 - x) * TMath::Log(1 - x);
        result = - 0.5 * log2 - fspens(x / (x - 1));
    } else if (x <= 0.5) {
        result = fspens(x);
    } else if (x <= 1) {
        Double_t logprod;
        logprod = TMath::Log(x) * TMath::Log(1 - x + 1E-10);
        result = f1 - logprod - fspens(1 - x);
    } else if (x <= 2) {
        Double_t logprod;
        logprod = TMath::Log(x) * TMath::Log((x - 1) * (x - 1) / x);
        result = f1 - 0.5  * logprod + fspens(1 - 1 / x);
    } else {
        Double_t log2;
        log2 = TMath::Log(x) * TMath::Log(x);
        result = 2 * f1 - 0.5 * log2 - fspens(1 / x);
    }

    return result;
}



Double_t fspens(const Double_t x)
{
    Double_t f   = 0;
    Double_t a   = 1;
    Double_t an  = 0;
    Double_t tch = 1E-16;
    Double_t b;

    do {
        an = an + 1;
        a  = a * x;
        b  = a / SQ(an);
        f  = f + b;
    } while (b - tch > 0);

    return f;
}



Double_t dfint(Int_t narg, double *arg, Int_t *nent, Double_t *ent, Double_t *table)
{
    Double_t kd = 1;
    Double_t m  = 1;
    Int_t    ja = 1;
    Int_t    jb;
    Int_t    k;

    Int_t    ncomb[narg];
    Double_t d[narg];
    Int_t    ient[narg];

    Double_t result;

    for (int i = 1; i <= narg; ++i) {
        ncomb[i-1] = 1;
        jb = ja - 1 +  nent[i-1];

        Int_t j = ja;
        while (j <= jb) {
            if (arg[i-1] <= ent[j-1]) goto exit_inner_loop;
            j++;
        }
        j = jb;
exit_inner_loop:
        if (j == ja) ++j;
        d[i-1] = (ent[j-1] - arg[i-1]) / (ent[j-1] - ent[j-2]);
        ient[i-1] = j - ja;
        kd = kd + ient[i-1] * m;
        m = m * nent[i-1];
        ja = jb + 1;
    }

    result = 0;

    Double_t fac;
    Int_t    iadr;
    Int_t    ifadr;
    Int_t    il;

do_again:
    fac   = 1;
    iadr  = kd;
    ifadr = 1;

    for (Int_t i = 1; i <= narg; ++i) {
        if (ncomb[i-1] == 0) {
            fac = fac * d[i-1];
            iadr  = iadr - ifadr;
        } else {
            fac = fac * (1 -  d[i-1]);
        }
        ifadr = ifadr * nent[i-1];
    }

    result = result + fac * table[iadr-1];

    il = narg;
    while (ncomb[il-1] == 0) {
        --il;
        if (il == 0) return result;
    }

    ncomb[il-1] = 0;
    if (il == narg) goto do_again;
    k = il + 1;
    while (k <= narg) {
        ncomb[k-1] = 1;
        k++;
    }
    goto do_again;
}

} // end namespace HapradUtils
