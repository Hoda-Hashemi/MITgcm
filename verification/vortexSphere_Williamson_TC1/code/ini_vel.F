#include "MYPACKAGE_OPTIONS.h"

CBOP
C     !ROUTINE: TC1_PRESCRIBE_VELOCITY
C     !INTERFACE:
      SUBROUTINE TC1_PRESCRIBE_VELOCITY( bi, bj, myThid )

C     !DESCRIPTION:
C     Prescribe the Williamson TC1 solid-body velocity on the local tile.
C     The flow is written from a streamfunction so it stays divergence-free
C     on the staggered grid used by MITgcm.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "MYPACKAGE.h"

C     !INPUT PARAMETERS:
      INTEGER bi, bj, myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i, j, k
      _RL alpha
      _RL u0
      _RL dLonDeg, dLatDeg, dLonRad, dLatRad
      _RL lonDeg, latDeg, lonW, lonE, latS, latN, dx, dy
      _RL psiN, psiS, psiW, psiE
      _RL tiny

      u0 = 2.D0 * PI * rSphere / (12.D0 * 86400.D0)
      alpha = myPa_param1
      dLonDeg = 360.D0 / DBLE(Nx)
      dLatDeg = 180.D0 / DBLE(Ny)
      dLonRad = dLonDeg * deg2rad
      dLatRad = dLatDeg * deg2rad
      dy = rSphere * dLatRad
      tiny = 1.D-14

      DO k = 1, Nr
       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
          lonDeg = xG(i,j,bi,bj)
          latDeg = yC(i,j,bi,bj)
          latS = latDeg - 0.5D0 * dLatDeg
          latN = latDeg + 0.5D0 * dLatDeg
          psiN = u0 * rSphere * (
     &      COS(alpha) * SIN(latN*deg2rad)
     &    - SIN(alpha) * COS(latN*deg2rad) * COS(lonDeg*deg2rad) )
          psiS = u0 * rSphere * (
     &      COS(alpha) * SIN(latS*deg2rad)
     &    - SIN(alpha) * COS(latS*deg2rad) * COS(lonDeg*deg2rad) )
          uVel(i,j,k,bi,bj) =
     &      ( psiN - psiS ) / dy
        ENDDO
       ENDDO

       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
          lonW = xC(i,j,bi,bj) - 0.5D0 * dLonDeg
          lonE = xC(i,j,bi,bj) + 0.5D0 * dLonDeg
          IF ( lonW .LT. 0.D0 ) lonW = lonW + 360.D0
          IF ( lonE .GE. 360.D0 ) lonE = lonE - 360.D0
          latDeg = yG(i,j,bi,bj)
          dx = rSphere * COS(latDeg*deg2rad) * dLonRad
          IF ( ABS(dx) .GT. tiny ) THEN
            psiW = u0 * rSphere * (
     &        COS(alpha) * SIN(latDeg*deg2rad)
     &      - SIN(alpha) * COS(latDeg*deg2rad) * COS(lonW*deg2rad) )
            psiE = u0 * rSphere * (
     &        COS(alpha) * SIN(latDeg*deg2rad)
     &      - SIN(alpha) * COS(latDeg*deg2rad) * COS(lonE*deg2rad) )
            vVel(i,j,k,bi,bj) =
     &        ( psiW - psiE ) / dx
          ELSE
            vVel(i,j,k,bi,bj) = 0.D0
          ENDIF
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END
