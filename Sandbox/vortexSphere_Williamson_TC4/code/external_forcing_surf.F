#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: EXTERNAL_FORCING_SURF
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_SURF(
     I             iMin, iMax, jMin, jMax,
     I             myTime, myIter, myThid )
C     !DESCRIPTION:
C     TC4 analytic translating-low forcing for one-layer z-coordinate runs.
CEOP

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "FFIELDS.h"
#include "GRID.h"
#include "SURFACE.h"

      INTEGER iMin, iMax
      INTEGER jMin, jMax
      _RL myTime
      INTEGER myIter
      INTEGER myThid

      INTEGER bi, bj
      INTEGER i, j
      _RL forcU, forcV, forcH
      _RL unusedU, unusedV
      _RL thick

      thick = drF(1)

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          surfaceForcingT(i,j,bi,bj) = 0. _d 0
          surfaceForcingS(i,j,bi,bj) = 0. _d 0
          Qnet(i,j,bi,bj) = 0.
          Qsw(i,j,bi,bj) = 0.
          saltFlux(i,j,bi,bj) = 0.
         ENDDO
        ENDDO

C       Continuity source: MITgcm uses dEta/dt = -div - EmPmR/rhoConst.
C       Williamson TC4 needs d h/dt + div(h u) = F_h, so EmPmR = -rho F_h.
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          CALL TC4_FORCING_POINT(
     I         xC(i,j,bi,bj), yC(i,j,bi,bj), myTime,
     O         unusedU, unusedV, forcH )
          EmPmR(i,j,bi,bj) =
     &      -rhoConst*forcH*maskInC(i,j,bi,bj)
         ENDDO
        ENDDO

        IF ( staggerTimeStep ) THEN
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           PmEpR(i,j,bi,bj) = -EmPmR(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDIF

C       Momentum sources are applied as surfaceForcing/drF by APPLY_FORCING.
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          CALL TC4_FORCING_POINT(
     I         xG(i,j,bi,bj), yC(i,j,bi,bj), myTime,
     O         forcU, unusedV, forcH )
          surfaceForcingU(i,j,bi,bj) =
     &      forcU*thick*maskInW(i,j,bi,bj)

          CALL TC4_FORCING_POINT(
     I         xC(i,j,bi,bj), yG(i,j,bi,bj), myTime,
     O         unusedU, forcV, forcH )
          surfaceForcingV(i,j,bi,bj) =
     &      forcV*thick*maskInS(i,j,bi,bj)
         ENDDO
        ENDDO

       ENDDO
      ENDDO

      RETURN
      END
