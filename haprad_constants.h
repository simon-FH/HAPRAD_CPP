#ifndef HAPRAD_CONSTANTS
#define HAPRAD_CONSTANTS

#include "TMath.h"

const Double_t kBarn = 0.389379E6;//HH: (h_bar*c)^2=0.389mb
const Double_t kAlpha = 0.729735E-2;//HH: alfa=1/137
const Double_t kPi = TMath::Pi();

const Double_t kMassProton = 0.938272;
const Double_t kMassNeutron = 0.93956536;
const Double_t kMassElectron = 0.000511;
const Double_t kMassMuon = 0.10565;
const Double_t kMassPion = 0.1395675;

const Double_t kMassDetectedHadron   = 0.1395675;
const Double_t kMassUndetectedHadron = 0.938272;//HH: 0.938272->0.93956536

const Double_t kMassC2 = TMath::Power(kMassProton + kMassDetectedHadron,2);//HH: kMassProton->kMassNeutron

const Double_t kRadianDeg = 57.2957795131;
const Double_t kEpsMachine = 1.1E-15;

#endif
