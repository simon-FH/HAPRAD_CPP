C     born_harmonics -- the Born cross section of HAPRAD 2.0 at many points, and
C     its three azimuthal harmonics.
C
C     For PLAN.md Phase 2.3: building a structure-function table from HAPRAD
C     2.0's own PDF x FF model, so that the C++ and the FORTRAN can be compared
C     with the SAME physics underneath.
C
C     BORNSET is the setup part of ihaprad() and BORNPHI the vv10/vv20 part of
C     sphih(), copied with the print statements removed; everything after them
C     (deltas, qqt -- the tail integrals) is left out, which is what makes this
C     fast. The copy is checked against the full calculation in haprad2_point.
C
C     The Born cross section is exactly B0 + Bc cos(phi) + Bcc cos(2 phi): the
C     structure functions do not depend on phi and the kinematic coefficients
C     are at most quadratic in cos(phi). Twelve equally spaced angles therefore
C     give the three harmonics exactly.
C
C     stdin : one point per line,   j0 j1 j2 j3  E  x  Q2  z  p_t
C     stdout: one result per line,  j0 j1 j2 j3  ok  B0  Bc  Bcc
C     (j0..j3 are cell indices, passed through untouched; i1/i2 are taken
C     by haprad_consts.inc.)
C             ok = 0 when ihaprad's kinematic checks would reject the point.

      PROGRAM born_harmonics
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INTEGER j0, j1, j2, j3, ok, k, nphi
      PARAMETER (nphi = 12)
      DOUBLE PRECISION eb, xx, qq, zz, pp, ph, sib, b0, bc, bcc, twopi

C     haprad_consts.inc truncates pi to 3.1415926; use the exact value for the
C     sampling angles so the projection is exact.
      twopi = 8.d0 * ATAN(1.d0)

C     Same masses as haprad2_point and the C++ haprad_constants.h.
      amhh = 0.1395675d0
      amhu = amp

 10   READ (*, *, END = 99) j0, j1, j2, j3, eb, xx, qq, zz, pp
      CALL bornset(eb, xx, qq, zz, pp, ok)
      b0 = 0.d0
      bc = 0.d0
      bcc = 0.d0
      IF (ok .EQ. 1) THEN
         DO k = 0, nphi - 1
            ph = twopi * DBLE(k) / DBLE(nphi)
            CALL bornphi(ph, sib)
            b0 = b0 + sib
            bc = bc + sib * COS(ph)
            bcc = bcc + sib * COS(2.d0 * ph)
         ENDDO
         b0 = b0 / DBLE(nphi)
         bc = 2.d0 * bc / DBLE(nphi)
         bcc = 2.d0 * bcc / DBLE(nphi)
      ENDIF
      WRITE (*, '(4(i9,1x),i2,3(1x,e24.16))') j0, j1, j2, j3, ok,
     &      b0, bc, bcc
      GOTO 10
 99   CONTINUE
      END


C     Setup part of ihaprad() for a p_t input, prints removed. ok = 1 if the
C     point passes the same kinematic checks.
      SUBROUTINE bornset(bmom, xmas, q2m, zmas, ptmas, ok)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'tail.inc'
      INCLUDE 'phi.inc'
      INCLUDE 'epsmarch.inc'
      DOUBLE PRECISION bmom, xmas, q2m, zmas, ptmas
      DOUBLE PRECISION snuc, yma, ymi, sqnuq, p22max, tdmax
      INTEGER ok

      ok = 0
      IF (ptmas .LT. 0.d0) RETURN

      ipol = 0
      iphi_had = 0
      iphi_rad = 0
      ilep = 1
      isf1 = 1
      isf2 = isf20
      isf3 = 1
      un = 1.d0
      pl = 0.d0
      pn = 0.d0
      qn = 0.d0

      snuc = 2.d0 * amp * bmom
      xs = xmas
      zdif = zmas
      tdif = ptmas
      y = q2m
      ys = y / (snuc * xs)
      yma = 1.d0 / (1.d0 + amp ** 2 * xs / snuc)
      amc2 = (amp + amhh) ** 2
      ymi = (amc2 - amp ** 2) / (snuc * (1.d0 - xs))
      IF (ys .GT. yma .OR. ys .LT. ymi .OR.
     &    xs .GT. 1.d0 .OR. xs .LT. 0.d0) RETURN

      CALL conkin(snuc)
      ehad = anu * zdif
      sqnuq = SQRT(anu ** 2 + y)
      IF (ehad .LT. amhh) RETURN
      pph = SQRT(ehad ** 2 - amhh ** 2)

      pth = tdif
      IF (pph .LT. pth) RETURN
      plh = SQRT(pph ** 2 - pth ** 2)
      IF (pph .GT. pth) THEN
         an = an * sqly / 2.d0 / amp / plh
      ELSE
         an = 0.d0
      ENDIF
      tdif = amhh ** 2 - y + 2.d0 * (sqnuq * plh - anu * ehad)

      p22max = w2 - amhh ** 2
      p22 = amp2 + sx * (1.d0 - zdif) + tdif
      IF (p22 .LT. amc2 .OR. p22 .GT. p22max) RETURN

      tdmin = amhh ** 2 - y + 2.d0 * ( sqnuq * pph - anu * ehad)
      tdmax = amhh ** 2 - y + 2.d0 * (-sqnuq * pph - anu * ehad)
      IF ((tdif - tdmin) .GT. epsmarch .OR. tdif .LT. tdmax) RETURN

      ok = 1
      RETURN
      END


C     The vv10/vv20 part of sphih() at angle phih, then bornin().
      SUBROUTINE bornphi(phih, sib)
      IMPLICIT NONE
      INCLUDE 'haprad_consts.inc'
      INCLUDE 'sxy.inc'
      INCLUDE 'phi.inc'
      DOUBLE PRECISION phih, sib, costs, sints, costx, sintx, d

      phidif = phih
      d = s * x * y - amp2 * y ** 2 - aml2 * aly
      costs = (s * (s - x) + 2.d0 * amp2 * y) / sqls / sqly
      costx = (x * (s - x) - 2.d0 * amp2 * y) / sqlx / sqly
      IF (d .GT. 0.d0) THEN
         sints = 2.d0 * amp * SQRT(d) / sqls / sqly
         sintx = 2.d0 * amp * SQRT(d) / sqlx / sqly
      ELSE
         sints = 0.d0
         sintx = 0.d0
      ENDIF
      vv10 = (s * ehad - sqls * (costs * plh + sints * pth *
     &        COS(phidif))) / amp
      vv20 = (x * ehad - sqlx * (costx * plh + sintx * pth *
     &        COS(phidif))) / amp
      CALL bornin(sib)
      RETURN
      END
