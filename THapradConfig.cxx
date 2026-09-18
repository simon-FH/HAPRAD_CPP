#include "THapradConfig.h"

// fLepton was never initialised. Nothing read it -- every site that needed the
// lepton mass switched on PolarizationType() instead (see below) -- so this was
// latent rather than active. 1 = electron, matching `ilep = 1` in fhaprad.f.
THapradConfig::THapradConfig()
    : fPolType(0), fLepton(1), fPhiHad(0), fPhiRad(0), fEpsPhiR(0.1), fEpsTau(0.001), fEpsRR(0.01) {
  // Do nothing
}
