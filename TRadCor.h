#ifndef TRADCOR_H
#define TRADCOR_H

#include "TROOT.h"

class THapradConfig;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;
class TSemiInclusiveModel;

class TRadCor {
 public:
  TRadCor();
  ~TRadCor();

  // UNITS: E, Q2, p_t, maxMx2 in GeV powers; phi in DEGREES (it is divided by
  // kRadianDeg internally). EG2 ntuples store PhiPQ in degrees, RG-E stores
  // phi_PQ in RADIANS -- passing the latter straight through shrinks every
  // angle by 180/pi with no warning.
  void CalculateRCFactor(Double_t E, Double_t x, Double_t Q2, Double_t z, Double_t p_t, Double_t phi, Double_t maxMx2, Double_t targProp);
  Double_t GetRCFactor(Double_t E, Double_t x, Double_t Q2, Double_t z, Double_t p_t, Double_t phi, Double_t maxMx2, Double_t targProp);

  // The Born cross section alone: same kinematics and normalisation as
  // CalculateRCFactor(), but no delta factors and no tail integrals, so it is
  // cheap. Returns 0 (and sets the kinematic-error state) if the point is
  // rejected. Same units and conventions as CalculateRCFactor; phi in degrees.
  Double_t CalculateBorn(Double_t E, Double_t x, Double_t Q2, Double_t z, Double_t p_t, Double_t phi);

  Double_t GetFactor1(void);
  Double_t GetFactor2(void);
  Double_t GetFactor3(void);

  Double_t GetSigBorn() { return sigma_born; }
  Double_t GetSigObs() { return sig_obs; }
  Double_t GetTail(Int_t i) { return tai[i]; }

  // Semi-inclusive structure-function table (see TSemiInclusiveModel.h for the
  // format). Without one, A = Ac = Acc = 0, so sigma_Born = 0 and every RC
  // factor is returned as 0. Returns false if the file cannot be used.
  Bool_t LoadSemiInclusiveTable(const char* path);
  const TSemiInclusiveModel* GetSemiInclusiveModel(void) const { return fSIModel; };

  void RegisteredLepton(Int_t type = 1);
  void IntegratePhiRad(Int_t type = 0);
  void IntegratePhiHad(Int_t type = 0);
  void SetPolarization(Int_t type = 0);

  const THapradConfig* GetConfig(void) const { return fConfig; };
  const TKinematicalVariables* GetKinematicalVariables(void) const { return fKin; };
  const TLorentzInvariants* GetLorentzInvariants(void) const { return fInv; };
  const THadronKinematics* GetHadronKinematics(void) const { return fHadKin; };

 private:
  void Initialization(void);
  void SPhiH(void);
  Double_t RadiativeTail(void);
  Double_t ExclusiveRadiativeTail(void);

  THapradConfig* fConfig;
  TKinematicalVariables* fKin;
  TLorentzInvariants* fInv;
  THadronKinematics* fHadKin;
  TSemiInclusiveModel* fSIModel;

  // Target Z and A proportion NAZ = Z/A
  Double_t NAZ;

  // Results
  Double_t sigma_born;  // sigma_0
  Double_t sig_obs;     // sigma_{obs}
  Double_t tai[2];

  // Integration
  Double_t N;  // Normalization factor

  // Kinematical error or not
  Bool_t fKinError;

  // File with parameters found or not
  Bool_t fParametersError;
};

#endif
