#include "MYPACKAGE_OPTIONS.h"

CBOP
C     !ROUTINE: TC1_PRESCRIBE_VELOCITY
C     !INTERFACE:
      SUBROUTINE TC1_PRESCRIBE_VELOCITY( bi, bj, myThid )

C     !DESCRIPTION:
C     Prescribe the Williamson TC1 solid-body velocity on the local tile.

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
      _RL lonRad, latRad

      u0 = 2.D0 * PI * rSphere / (12.D0 * 86400.D0)
      alpha = myPa_param1

      DO k = 1, Nr
       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
          latRad = yC(i,j,bi,bj) * deg2rad
          lonRad = xG(i,j,bi,bj) * deg2rad
          uVel(i,j,k,bi,bj) = u0 * (
     &      COS(latRad) * COS(alpha)
     &    + SIN(latRad) * COS(lonRad) * SIN(alpha) )
        ENDDO
       ENDDO

       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
          latRad = yG(i,j,bi,bj) * deg2rad
          lonRad = xC(i,j,bi,bj) * deg2rad
          vVel(i,j,k,bi,bj) = -u0 * SIN(lonRad) * SIN(alpha)
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END
