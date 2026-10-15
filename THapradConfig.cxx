#include "THapradConfig.h"

// Defaults: unpolarised target; electron (1, as `ilep = 1` in fhaprad.f); no
// integration over phi_had; phi_rad approximated rather than integrated.
THapradConfig::THapradConfig()
    : fPolType(0), fLepton(1), fPhiHad(0), fPhiRad(0), fEpsPhiR(0.1), fEpsTau(0.001), fEpsRR(0.01) {
  // Do nothing
}
