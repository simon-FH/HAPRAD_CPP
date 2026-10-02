// rc_point -- evaluate ONE kinematic point with the C++ HAPRAD and dump both
// the Lorentz/hadron kinematics (Tier 1) and the radiative-correction results
// (Tier 3) in a machine-readable form.
//
// Every reportable quantity is printed as
//
//     @<name> <value>
//
// so the surrounding library chatter (which goes to stdout/stderr uninvited)
// can be filtered out by the comparison script without ambiguity.
//
// The kinematics are obtained by driving TLorentzInvariants / THadronKinematics
// directly -- the same objects, constructed the same way, that TRadCor builds
// internally. They cannot be read back from TRadCor after the fact because
// CalculateRCFactor() deletes them before returning.

#include "TRadCor.h"
#include "THapradConfig.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "THapradException.h"
#include "TRV2TR.h"
#include "TSemiInclusiveModel.h"
#include "haprad_constants.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

static void emit(const char* key, double value) { printf("@%-10s % .14e\n", key, value); }

int main(int argc, char** argv) {
  if (argc < 7) {
    fprintf(stderr,
            "usage: rc_point E x Q2 z pt phi_deg [Mx2_threshold]\n"
            "  Mx2_threshold defaults to 0, which disables TRadCor's missing-mass\n"
            "  gate so that the point is always evaluated. This matches the\n"
            "  Fortran reference driver, which calls ihaprad() directly. Pass\n"
            "  (m_n+m_pi)^2 = 1.16168 to reproduce what GetRC uses in production.\n");
    return 2;
  }

  const double E = atof(argv[1]);
  const double x = atof(argv[2]);
  const double Q2 = atof(argv[3]);
  const double z = atof(argv[4]);
  const double pt = atof(argv[5]);
  const double phi = atof(argv[6]);
  const double m2th = (argc > 7) ? atof(argv[7]) : 0.0;

  // ---------------------------------------------------------------- Tier 1
  int kin_ok = 1;
  {
    THapradConfig cfg;
    TKinematicalVariables kin(x, -Q2, z, pt, phi / kRadianDeg, E);
    TLorentzInvariants inv(&cfg, &kin);
    THadronKinematics had(&cfg, &kin, &inv);

    try {
      inv.Evaluate();
      had.Evaluate();
    } catch (TKinematicException&) {
      kin_ok = 0;
    }

    printf("### BEGIN KIN\n");
    emit("kin_ok", kin_ok);
    if (kin_ok) {
      // Lorentz invariants
      emit("S", inv.S());
      emit("X", inv.X());
      emit("Sx", inv.Sx());
      emit("Sp", inv.Sp());
      emit("Q2", inv.Q2());
      emit("W2", inv.W2());
      emit("lambda_s", inv.LambdaS());
      emit("lambda_x", inv.LambdaX());
      emit("lambda_q", inv.LambdaQ());
      emit("lambda_m", inv.LambdaM());
      emit("sqrt_lq", inv.SqrtLq());
      // y is overwritten by SetQ2(); report what the code actually used
      emit("y", kin.Y());
      emit("t", kin.T());
      // Hadron kinematics
      emit("nu", had.Nu());
      emit("Eh", had.Eh());
      emit("ph", had.Ph());
      emit("pl", had.Pl());
      emit("pt", had.Pt());
      emit("sqnuq", had.SqNuQ());
      emit("px2", had.Px2());
      emit("V1", had.V1());
      emit("V2", had.V2());
      // tau limits, as TQQTPhi/TRV2TR compute them
      const double M2 = kMassProton * kMassProton;
      const double tau_max = (inv.Sx() + inv.SqrtLq()) / (2. * M2);
      emit("tau_max", tau_max);
      emit("tau_min", -inv.Q2() / M2 / tau_max);
    }
    printf("### END KIN\n");
  }

  // ---------------------------------------------------------------- Tier 3
  printf("### BEGIN RES\n");
  {
    TRadCor rc;
    // Structure-function table: HAPRAD_SI_TABLE=<path>. Without one the
    // amplitudes are zero, and so are sigma_Born and every RC factor.
    if (const char* table = getenv("HAPRAD_SI_TABLE")) {
      if (!rc.LoadSemiInclusiveTable(table)) return 3;
      // HAPRAD_SI_NEAREST=1: nearest-cell lookup instead of interpolation, for
      // comparison only.
      if (getenv("HAPRAD_SI_NEAREST"))
        const_cast<TSemiInclusiveModel*>(rc.GetSemiInclusiveModel())->SetInterpolation(false);
    }
    rc.CalculateRCFactor(E, x, Q2, z, pt, phi, m2th, 0.5);

    // Tier 2 diagnostic: dump the exclusive-tail integrand rv2tr(tau, phi_k)
    // itself, so it can be compared against the Fortran pointwise. This
    // separates a mis-ported integrand from a quadrature problem. It relies on
    // TRadCor keeping fKin/fInv/fHadKin alive after the call.
    if (getenv("HAPRAD_DUMP_EXC") && kin_ok) {
      const double M2 = kMassProton * kMassProton;
      const TLorentzInvariants* iv = rc.GetLorentzInvariants();
      const double tau_max = (iv->Sx() + iv->SqrtLq()) / (2. * M2);
      const double tau_min = -iv->Q2() / M2 / tau_max;

      TRV2TR f(&rc);
      printf("### BEGIN EXC\n");
      const int NT = 9, NP = 5;
      for (int i = 1; i <= NT; ++i) {
        for (int j = 0; j < NP; ++j) {
          double arg[2];
          arg[0] = tau_min + (tau_max - tau_min) * i / double(NT + 1);
          arg[1] = 2. * kPi * j / double(NP);
          char key[32];
          snprintf(key, sizeof(key), "t%dp%d", i, j);
          printf("@%-10s % .14e   tau=% .10e phi=% .10e\n", key, f(arg), arg[0], arg[1]);
        }
      }
      printf("### END EXC\n");
    }

    emit("kin_ok", kin_ok);
    emit("sib", rc.GetSigBorn());
    emit("sig_obs", rc.GetSigObs());
    emit("tai_in", rc.GetTail(0));
    emit("tai_ex", rc.GetTail(1));
    emit("f1", rc.GetFactor1());
    emit("f2", rc.GetFactor2());
    emit("f3", rc.GetFactor3());
    // Table diagnostics for this call: how many structure-function lookups the
    // integrals made, and how many fell outside the table or in an empty cell.
    if (const TSemiInclusiveModel* m = rc.GetSemiInclusiveModel()) {
      emit("sf_lookups", m->NLookups());
      emit("sf_oor", m->NOutOfRange());
      emit("sf_empty", m->NEmptyCell());
      emit("sf_unphys", m->NUnphysical());
    }
  }
  printf("### END RES\n");

  return 0;
}
