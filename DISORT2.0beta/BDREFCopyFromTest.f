c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: BDREF.f,v 2.1 2000/03/27 21:40:51 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      REAL FUNCTION  BDREF( WVNMLO, WVNMHI, MU, MUP, DPHI )

c      Supplies surface bi-directional reflectivity.
c
c      This is only a "stub" version. The user must replace this
c      by his/her own BDREF function.
c
c
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c              angles) are positive.
c
c  INPUT:
c
c    WVNMLO : Lower wavenumber (inv cm) of spectral interval
c
c    WVNMHI : Upper wavenumber (inv cm) of spectral interval
c
c    MU     : Cosine of angle of reflection (positive)
c
c    MUP    : Cosine of angle of incidence (positive)
c
c    DPHI   : Difference of azimuth angles of incidence and reflection
c                (radians)
c
c
c   Called by- DREF, SURFAC

c +-------------------------------------------------------------------+
c
c     .. Scalar Arguments ..

      REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
c     ..
c     .. Local Scalars ..

      INTEGER   IREF
      REAL      B0, B, CTHETA, GAMMA, H0, H, HH, P, THETA, W
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC COS, SQRT
c     ..

      IREF = 1

      IF ( IREF.EQ.1 ) THEN

c                              ** Hapke's BRDF model (times Pi/Mu0)
c                              ** (Hapke, B., Theory of reflectance
c                              ** and emittance spectroscopy, Cambridge
c                              ** University Press, 1993, Eq. 8.89 on
c                              ** page 233. Parameters are from
c                              ** Fig. 8.15 on page 231, expect for w.)

         CTHETA = MU * MUP + (1.-MU**2)**.5 * (1.-MUP**2)**.5
     &            * COS( DPHI )
         THETA = ACOS( CTHETA )

         P    = 1. + 0.5 * CTHETA

         HH   = 0.06
         B0   = 1.0
         B    = B0 * HH / ( HH + TAN( THETA/2.) )

         W = 0.6
         GAMMA = SQRT( 1. - W )
         H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
         H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )

         BDREF = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

      ELSE

         CALL ERRMSG( 'BDREF--Need to supply surface BDRF model',
     &                 .TRUE.)

      END IF

      RETURN
      END
