C     haprad2_point -- evaluate ONE kinematic point with the ORIGINAL HAPRAD 2.0
C     Fortran and dump the same quantities, in the same "@name value" format,
C     as the C++ driver rc_point.
C
C     This replaces rcdat.f as the main program. It calls ihaprad() directly
C     rather than going through fhaprad(), so that (a) no 1e-3 rescaling is
C     applied -- the C++ GetSigBorn() is unscaled -- and (b) the missing-mass
C     gate in fhaprad does not fire, matching rc_point's default Mx2 threshold
C     of 0. Kinematic rejection is still done by ihaprad itself, identically
C     to the C++ TKinematicException path.

      PROGRAM haprad2_point
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'phi.inc'

      DOUBLE PRECISION Ebeam, xb, q2in, zin, ptin, phiin
      DOUBLE PRECISION sib, sig, delinf, delta, tai(3)
      DOUBLE PRECISION epsphir, epstau, epsrr, raddeg
      DOUBLE PRECISION tau_mx, tau_mn
      INTEGER ilepl, iphi_radl, nargs, kin_ok
      CHARACTER*64 arg

      DATA raddeg/57.2957795131d0/

      nargs = IARGC()
      IF (nargs .LT. 6) THEN
         WRITE (0,*) 'usage: haprad2_point E x Q2 z pt phi_deg'
         CALL EXIT(2)
      ENDIF
      CALL GETARG(1, arg)
      READ (arg,*) Ebeam
      CALL GETARG(2, arg)
      READ (arg,*) xb
      CALL GETARG(3, arg)
      READ (arg,*) q2in
      CALL GETARG(4, arg)
      READ (arg,*) zin
      CALL GETARG(5, arg)
      READ (arg,*) ptin
      CALL GETARG(6, arg)
      READ (arg,*) phiin

C     Masses, matching haprad_constants.h on the C++ side:
C       kMassDetectedHadron   = 0.1395675
C       kMassUndetectedHadron = 0.938272  (= amp)
      amhh = 0.1395675d0
      amhu = amp

C     Same integration settings fhaprad.f uses, and the same ones
C     THapradConfig defaults to (EpsTau 1e-3, EpsRR 1e-2); note that the C++
C     THapradConfig default for EpsPhiR is 0.1, but it is unused when
C     iphi_rad = 0.
      ilepl     = 1
      iphi_radl = 0
      epsphir   = 0.01d0
      epstau    = 0.001d0
      epsrr     = 0.01d0

C     Sentinel: ihaprad returns early, leaving sib untouched, when it rejects
C     the point. There is no status flag to inspect.
      sib = -987654321d0
      sig = 0d0
      tai(1) = 0d0
      tai(2) = 0d0
      tai(3) = 0d0

      CALL ihaprad(Ebeam, ilepl, iphi_radl, epsphir, epstau, epsrr,
     &             xb, -q2in, zin, ptin, phiin/raddeg,
     &             sib, sig, delinf, delta, tai)

      kin_ok = 1
      IF (sib .EQ. -987654321d0) THEN
         kin_ok = 0
         sib = 0d0
      ENDIF

      WRITE (*,'(a)') '### BEGIN KIN'
      CALL emitd('kin_ok', DBLE(kin_ok))
      IF (kin_ok .EQ. 1) THEN
C        /sxy/ : note Fortran 'x' is the invariant X = 2*k2*p, and 'y' is Q^2.
         CALL emitd('S', s)
         CALL emitd('X', x)
         CALL emitd('Sx', sx)
         CALL emitd('Sp', sxp)
         CALL emitd('Q2', y)
         CALL emitd('W2', w2)
         CALL emitd('lambda_s', als)
         CALL emitd('lambda_x', alx)
         CALL emitd('lambda_q', aly)
         CALL emitd('lambda_m', alm)
         CALL emitd('sqrt_lq', sqly)
         CALL emitd('y', ys)
         CALL emitd('t', tdif)
C        /phi/
         CALL emitd('nu', anu)
         CALL emitd('Eh', ehad)
         CALL emitd('ph', pph)
         CALL emitd('pl', plh)
         CALL emitd('pt', pth)
         CALL emitd('sqnuq', SQRT(anu**2 + y))
         CALL emitd('px2', p22)
         CALL emitd('V1', vv10)
         CALL emitd('V2', vv20)
         tau_mx = (sx + sqly) / (2d0 * amp2)
         tau_mn = -y / amp2 / tau_mx
         CALL emitd('tau_max', tau_mx)
         CALL emitd('tau_min', tau_mn)
      ENDIF
      WRITE (*,'(a)') '### END KIN'

      WRITE (*,'(a)') '### BEGIN RES'
      CALL emitd('kin_ok', DBLE(kin_ok))
      CALL emitd('sib', sib)
C     ihaprad returns sig already including both tails; the C++ GetSigObs()
C     excludes them, so subtract to get the comparable quantity.
      CALL emitd('sig_obs', sig - tai(1) - tai(2))
      CALL emitd('tai_in', tai(1))
      CALL emitd('tai_ex', tai(2))
      IF (sib .NE. 0d0) THEN
         CALL emitd('f1', (sig - tai(2)) / sib)
         CALL emitd('f2', sig / sib)
         CALL emitd('f3', (sig - tai(2) + 0.5d0*tai(2)) / sib)
      ELSE
         CALL emitd('f1', 0d0)
         CALL emitd('f2', 0d0)
         CALL emitd('f3', 0d0)
      ENDIF
      WRITE (*,'(a)') '### END RES'

      END

      SUBROUTINE emitd(name, value)
      IMPLICIT NONE
      CHARACTER*(*) name
      DOUBLE PRECISION value
      CHARACTER*10 pad
      pad = name
      WRITE (*,'(a1,a10,1x,e22.14)') '@', pad, value
      RETURN
      END
