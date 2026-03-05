#include "CPP_OPTIONS.h"

CBOP
C--   !ROUTINE: HASHEMI_CG2D_STARTUP_SOLVE
      SUBROUTINE HASHEMI_CG2D_STARTUP_SOLVE( myThid )
      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "CG2D.h"

      INTEGER myThid
      INTEGER numIters, nIterMin
      _RL firstResidual, minResidualSq, lastResidual
      _RL cg2d_b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL cg2d_x(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      ! _BEGIN_MASTER(myThid)
      !   WRITE(standardMessageUnit,'(A,1PE15.6)') 'omega = ', omega
      !   WRITE(standardMessageUnit,'(A,1PE15.6)') 'rSphere = ', rSphere
      !   WRITE(standardMessageUnit,'(A,1PE15.6)') 'gravity = ', gravity
      ! _END_MASTER(myThid)

      numIters = cg2dMaxIters
      IF ( numIters .LE. 0 ) numIters = 1000
      nIterMin = -1

C--   SOLVERS AND PRECONDITIONERS (FOR OUR EQUAITON)

      CALL HASHEMI_FILL_CG2D_SYSTEM( cg2d_b, cg2d_x, myThid )

C--   Select one preconditioner builder:
C--   Old MITgcm-style preconditioner (kept for reference):

C--   True CG-ADI path: Krylov iteration uses ADI line solves as
C--   preconditioning operator M^{-1}(r), not only static pC/pW/pS weights.
!       CALL HASHEMI_SOLVE_CG_ADI(
!      U     cg2d_b, cg2d_x,
!      O     firstResidual, minResidualSq, lastResidual,
!      U     numIters, nIterMin,
!      I     myThid )

C--   Baseline MITgcm CG2D path (kept for reference):
      CALL HASHEMI_BUILD_CG2D_PRECOND( myThid )
      CALL CG2D( cg2d_b, cg2d_x, firstResidual, minResidualSq,
     & lastResidual, numIters, nIterMin, myThid )



      ! _BEGIN_MASTER(myThid)
      !  WRITE(standardMessageUnit,'(A,1PE15.6)') 'Initial Res:', firstResidual
      !  WRITE(standardMessageUnit,'(A,1PE15.6)') 'Final Res:  ', lastResidual
      !  WRITE(standardMessageUnit,'(A,I5)')      'Iters:      ', numIters
      ! _END_MASTER(myThid)

C--   Fill MITgcm state velocities from solved psi (diagnostic mapping).
      CALL HASHEMI_FILL_UV_FROM_PSI( cg2d_x, myThid )
      CALL HASHEMI_WRITE_UV_CSV( myThid )

      CALL HASHEMI_WRITE_PSI_CSV( cg2d_x, myThid )

      RETURN
      END

C--   !ROUTINE: HASHEMI_FILL_CG2D_SYSTEM
      SUBROUTINE HASHEMI_FILL_CG2D_SYSTEM( cg2d_b, cg2d_x, myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "CG2D.h"

      INTEGER myThid
      _RL cg2d_b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL cg2d_x(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      INTEGER bi, bj, iLat, iLon, jG
      _RL dLat, dLat2, lat, sinLat, cosLat, cosLat2
      _RL latMHalf, latPHalf, cosLatMHalf, cosLatPHalf
      _RL R2, coriolisTerm, coriolisTerm2
      _RL dC, dW, dE, dS, dN
      _RL coefWest, coefEast, coefSouth, coefNorth
      LOGICAL isOceanC, isOceanW, isOceanE, isOceanS, isOceanN

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO iLat=1-OLy,sNy+OLy
         DO iLon=1-OLx,sNx+OLx
            !!initialize the coefficents
          aW2d(iLon,iLat,bi,bj) = zeroRL
          aS2d(iLon,iLat,bi,bj) = zeroRL
          aC2d(iLon,iLat,bi,bj) = zeroRL
          cg2d_b(iLon,iLat,bi,bj) = zeroRL
          cg2d_x(iLon,iLat,bi,bj) = zeroRL
         ENDDO
        ENDDO

        DO iLat=1,sNy
         jG = myYGlobalLo + (bj-1)*sNy + iLat - 1
         DO iLon=1,sNx
          lat = yC(iLon,iLat,bi,bj) * deg2rad
          sinLat = SIN( lat )
          cosLat = COS( lat )
          cosLat2 = cosLat*cosLat

      !     dLat = delY(jG) * deg2rad
          dLat = dyF(iLon,iLat,bi,bj) / rSphere
          dLat2 = dLat*dLat
          latMHalf = lat - halfRL*dLat
          latPHalf = lat + halfRL*dLat

C--         Pole handling:
C--         At the first/last global latitude row, set cos(phi_{j-1/2}) or
C--         cos(phi_{j+1/2}) to zero so meridional face fluxes vanish. This is a
C--         valid FV pole decoupling for lat-lon grids (zero normal flux at poles).
            !!poles decoupling
          IF ( jG .EQ. 1 ) THEN
             cosLatMHalf = zeroRL
          ELSE
             cosLatMHalf = COS(latMHalf)
          ENDIF
          IF ( jG .EQ. Ny ) THEN
             cosLatPHalf = zeroRL
          ELSE
             cosLatPHalf = COS(latPHalf)
          ENDIF
          !!coriolis
          coriolisTerm = twoRL * omega * sinLat
          coriolisTerm2 = coriolisTerm*coriolisTerm
          !!ALTERED assigning depths
      !     dC = MAX( zeroRL, -R_low(iLon  ,iLat  ,bi,bj) )
      !     dW = MAX( zeroRL, -R_low(iLon-1,iLat  ,bi,bj) )
      !     dE = MAX( zeroRL, -R_low(iLon+1,iLat  ,bi,bj) )
      !     dS = MAX( zeroRL, -R_low(iLon  ,iLat-1,bi,bj) )
      !     dN = MAX( zeroRL, -R_low(iLon  ,iLat+1,bi,bj) )

          dC = oneRL
          dW = oneRL
          dE = oneRL
          dS = oneRL
          dN = oneRL

C--       Land/ocean handling:
C--       Ocean points is from bathymetry depth d=MAX(0,-R_low). Dry points are
C--       forced to an identity row (aC2d=1, RHS=0), so psi stays fixed at 0.
          isOceanC = dC .GT. zeroRL
          isOceanW = dW .GT. zeroRL
          isOceanE = dE .GT. zeroRL
          isOceanS = dS .GT. zeroRL
          isOceanN = dN .GT. zeroRL

C--       initializing the coefficients
          coefWest = zeroRL
          coefEast = zeroRL
          coefSouth = zeroRL
          coefNorth = zeroRL

          IF ( isOceanC .AND. isOceanW ) coefWest = oneRL/cosLat 
          IF ( isOceanC .AND. isOceanE ) coefEast = oneRL/cosLat 
          IF ( isOceanC .AND. isOceanS ) coefSouth = cosLatMHalf
          IF ( isOceanC .AND. isOceanN ) coefNorth = cosLatPHalf

          IF ( isOceanC ) THEN
           aW2d(iLon  ,iLat  ,bi,bj) = coefWest
           aW2d(iLon+1,iLat  ,bi,bj) = coefEast
           aS2d(iLon  ,iLat  ,bi,bj) = coefSouth
           aS2d(iLon  ,iLat+1,bi,bj) = coefNorth

           R2 = rSphere*rSphere
C--        Main diagonal is filled in a second pass from final face-coeffs
           aC2d(iLon,iLat,bi,bj) = zeroRL

           cg2d_b(iLon,iLat,bi,bj) = - coriolisTerm
     &      * (R2*dLat2*cosLat)
          ELSE
           aC2d(iLon,iLat,bi,bj) = oneRL 
           cg2d_b(iLon,iLat,bi,bj) = zeroRL
          ENDIF

         ENDDO
        ENDDO

C--     Build main diagonal from final off-diagonal values.
        R2 = rSphere*rSphere
        DO iLat=1,sNy
         jG = myYGlobalLo + (bj-1)*sNy + iLat - 1
         DO iLon=1,sNx
            !!ALTERED!!!!
      !     dC = MAX( zeroRL, -R_low(iLon  ,iLat  ,bi,bj) )

          dC = oneRL
          isOceanC = dC .GT. zeroRL
          IF ( isOceanC ) THEN
           lat = yC(iLon,iLat,bi,bj) * deg2rad
           sinLat = SIN( lat )
           cosLat = COS( lat )
           
           cosLat2 = cosLat*cosLat

           dLat = dyF(iLon,iLat,bi,bj) / rSphere
           dLat2 = dLat*dLat

           coriolisTerm = twoRL * omega * sinLat
           coriolisTerm2 = coriolisTerm*coriolisTerm

           aC2d(iLon,iLat,bi,bj) = -(
     &       aW2d(iLon  ,iLat  ,bi,bj) + aW2d(iLon+1,iLat  ,bi,bj)
     &      +aS2d(iLon  ,iLat  ,bi,bj) + aS2d(iLon  ,iLat+1,bi,bj) )
     &      - (coriolisTerm2/(gravity*dC)) * (R2*dLat2*cosLat)

          ELSE
           aC2d(iLon,iLat,bi,bj) = oneRL
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL EXCH_UV_XY_RS( aW2d, aS2d, .FALSE., myThid )
      CALL EXCH_XY_RS( aC2d, myThid )
      RETURN
      END

C--     !!ROUTINE: HASHEMI_BUILD_CG2D_PRECOND
C--    This routine builds a custom preconditioner (pC, pW, pS) specifically 
C--    for the Hashemi psi-equation. 
C--
C--    In standard MITgcm, "ini_cg2d/update_cg2d" builds a preconditioner 
C--    based on the pressure matrix. Because we have replaced the matrix 
C--    with custom coefficients (aW2d, aS2d, aC2d), we must rebuild the 
C--    preconditioner to match. This ensures the CG2D solver "sees" 
C--    our custom operator, preventing divergence or slow convergence.

      SUBROUTINE HASHEMI_BUILD_CG2D_PRECOND( myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "CG2D.h"
      INTEGER myThid, bi, bj, i, j
      _RS aC, aCw, aCs
      _RL tol

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          pC(i,j,bi,bj) = zeroRL
          pW(i,j,bi,bj) = zeroRL
          pS(i,j,bi,bj) = zeroRL
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      _BEGIN_MASTER(myThid)
      cg2dNorm = oneRL
      cg2dNormaliseRHS = .FALSE.
      tol = cg2dTargetResidual
      IF ( tol .LE. zeroRL ) tol = 1.D-13
      cg2dTolerance_sq = tol*tol
      _END_MASTER(myThid)

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          aC  = aC2d( i, j, bi,bj)
          aCs = aC2d( i,j-1,bi,bj)
          aCw = aC2d(i-1,j, bi,bj)
          IF ( aC .EQ. zeroRL ) THEN
           pC(i,j,bi,bj) = oneRL
          ELSE
           pC(i,j,bi,bj) = oneRL / aC
          ENDIF
          IF ( aC + aCw .EQ. zeroRL ) THEN
           pW(i,j,bi,bj) = zeroRL
          ELSE
           pW(i,j,bi,bj) = -aW2d(i,j,bi,bj) / 
     &          ((cg2dpcOffDFac*(aCw+aC))**2)
          ENDIF
          IF ( aC + aCs .EQ. zeroRL ) THEN
           pS(i,j,bi,bj) = zeroRL
          ELSE
           pS(i,j,bi,bj) = -aS2d(i,j,bi,bj) / 
     &          ((cg2dpcOffDFac*(aCs+aC))**2)
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL EXCH_XY_RS( pC, myThid )
      CALL EXCH_UV_XY_RS( pW, pS, .FALSE., myThid )
      RETURN
      END


C--   !ROCKET: CALCULATING THE VELOCITIES:                          
C--   Map solved streamfunction psi (cg2d_x) to MITgcm uVel/vVel fields.
C--
C--   C-grid note:
C--   - uVel is zonal velocity on W-faces, approximated from -d(psi)/dy.
C--   - vVel is meridional velocity on S-faces, approximated from +d(psi)/dx.
C--   - We use local metric distances (dyC, dxC) and apply face masks
C--     (_maskW, _maskS) so land faces remain zero.
      SUBROUTINE HASHEMI_FILL_UV_FROM_PSI( cg2d_x, myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

      INTEGER myThid
      _RL cg2d_x(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER bi, bj, i, j, k, jG
      _RL dPsiDy, dPsiDx, denY, denX
      LOGICAL printedDenY, wetS0, wetSm, wetSp

C--   Fill overlap of psi first (needed for i-1/i+1 and j-1/j+1 stencils).  Because calculating a derivative at point requires looking at neighbors (like or ), the code first ensures that data from neighboring "tiles" (CPU cores) is shared so the edges of each grid patch are up to date.
      CALL EXCH_XY_RL( cg2d_x, myThid )
      printedDenY = .FALSE.

C--   Reset velocity arrays in this thread tile range.
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           uVel(i,j,k,bi,bj) = zeroRL
           vVel(i,j,k,bi,bj) = zeroRL
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C--   Compute k=1 velocities from psi.
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         jG = myYGlobalLo + (bj-1)*sNy + j - 1
         DO i=1,sNx
C--       Meridional derivative for u = -d(psi)/dy (metrics are included):
C--        one-sided at global poles, centered otherwise.
          IF ( jG .EQ. 1 ) THEN !!south pole - Forward Difference
           denY = dyC(i,j+1,bi,bj)
           IF ( denY .GT. zeroRL ) THEN
            dPsiDy = ( cg2d_x(i,j+1,bi,bj) - cg2d_x(i,j,bi,bj) )/denY
           ELSE
            dPsiDy = zeroRL
           ENDIF
          ELSEIF ( jG .EQ. Ny ) THEN  !!north pole - Backward Difference
           denY = dyC(i,j,bi,bj)
           IF ( denY .GT. zeroRL ) THEN
            dPsiDy = ( cg2d_x(i,j,bi,bj) - cg2d_x(i,j-1,bi,bj) )/denY
           ELSE
            dPsiDy = zeroRL
           ENDIF
          ELSE
           denY = dyC(i,j,bi,bj) + dyC(i,j+1,bi,bj)

           IF ( denY .GT. zeroRL ) THEN !!interior points - Central Difference
C--         BUGTRACKING: Print one representative denY sample (metric scaling check).
            IF ( .NOT. printedDenY ) THEN
             WRITE(standardMessageUnit,
     & '(A,2I6,3(1X,1PE15.7))')
     & 'BUGTRACKING denY(i,j),dyC(j),dyC(j+1):',
     & i,j, denY, dyC(i,j,bi,bj), dyC(i,j+1,bi,bj)
             printedDenY = .TRUE.
            ENDIF
C--         BUGTRACKING: Wet-aware centered stencil test for coastal-noise diagnosis.
            IF ( _maskW(i,j+1,1,bi,bj) .EQ. oneRL .AND.
     &           _maskW(i,j-1,1,bi,bj) .EQ. oneRL ) THEN
             dPsiDy = ( cg2d_x(i,j+1,bi,bj) - cg2d_x(i,j-1,bi,bj) )
     &                /denY
            ELSE
             dPsiDy = zeroRL
            ENDIF
           ELSE
            dPsiDy = zeroRL
           ENDIF
          ENDIF

C--        Zonal derivative for v = +d(psi)/dx (periodic handled via overlap).
C--        Wet-aware stencil: centered if both x-neighbor S-faces are wet,
C--        otherwise one-sided; if no wet support, set to zero.
          wetS0 = _maskS(i  ,j,1,bi,bj) .EQ. oneRL
          wetSm = _maskS(i-1,j,1,bi,bj) .EQ. oneRL
          wetSp = _maskS(i+1,j,1,bi,bj) .EQ. oneRL

          IF ( .NOT. wetS0 ) THEN
           dPsiDx = zeroRL
          ELSEIF ( wetSm .AND. wetSp ) THEN
           denX = dxC(i,j,bi,bj) + dxC(i+1,j,bi,bj)
           IF ( denX .GT. zeroRL ) THEN
            dPsiDx = ( cg2d_x(i+1,j,bi,bj) - cg2d_x(i-1,j,bi,bj) )/denX
           ELSE
            dPsiDx = zeroRL
           ENDIF
          ELSEIF ( wetSp ) THEN
           denX = dxC(i+1,j,bi,bj)
           IF ( denX .GT. zeroRL ) THEN
            dPsiDx = ( cg2d_x(i+1,j,bi,bj) - cg2d_x(i,j,bi,bj) )/denX
           ELSE
            dPsiDx = zeroRL
           ENDIF
          ELSEIF ( wetSm ) THEN
           denX = dxC(i,j,bi,bj)
           IF ( denX .GT. zeroRL ) THEN
            dPsiDx = ( cg2d_x(i,j,bi,bj) - cg2d_x(i-1,j,bi,bj) )/denX
           ELSE
            dPsiDx = zeroRL
           ENDIF
          ELSE
           dPsiDx = zeroRL
          ENDIF

          uVel(i,j,1,bi,bj) = -dPsiDy * _maskW(i,j,1,bi,bj)
          vVel(i,j,1,bi,bj) =  dPsiDx * _maskS(i,j,1,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C--   Exchange halo for staggered velocity fields.
      CALL EXCH_UV_3D_RL( uVel, vVel, .TRUE., Nr, myThid )
      RETURN
      END


C--   !ROUTINE: HASHEMI_SOLVE_CG_ADI
C--   Conjugate-gradient solver with ADI line-solve preconditioning.
C--   This is the "true CG-ADI" path: at each Krylov iteration, apply
C--   M^{-1}(r) via alternating implicit row/column tridiagonal solves.
      SUBROUTINE HASHEMI_SOLVE_CG_ADI(
     U     cg2d_b, cg2d_x,
     O     firstResidual, minResidualSq, lastResidual,
     U     numIters, nIterMin,
     I     myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "CG2D.h"
      _RL cg2d_b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL cg2d_x(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL firstResidual, minResidualSq, lastResidual
      INTEGER numIters, nIterMin, myThid

      INTEGER bi, bj, i, j, it, actualIts
      _RL r(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL z(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL q(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL resSq, rhsSq, rz, rzNew, pA
      _RL alphaCG, betaCG, tol, stopSq, eps
      LOGICAL printResidual
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      eps = 1.D-30
      tol = cg2dTargetResidual
      IF ( tol .LE. zeroRL ) tol = 1.D-13

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          r(i,j,bi,bj) = zeroRL
          z(i,j,bi,bj) = zeroRL
          p(i,j,bi,bj) = zeroRL
          q(i,j,bi,bj) = zeroRL
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL EXCH_XY_RL( cg2d_x, myThid )
      CALL HASHEMI_APPLY_A(
     I      cg2d_x,
     O      q,
     I      myThid )

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          r(i,j,bi,bj) = cg2d_b(i,j,bi,bj) - q(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL HASHEMI_DOT_RL(
     I      r, r,
     O      resSq,
     I      myThid )
      CALL HASHEMI_DOT_RL(
     I      cg2d_b, cg2d_b,
     O      rhsSq,
     I      myThid )

      firstResidual = SQRT(MAX(resSq,zeroRL))
      lastResidual = firstResidual
      minResidualSq = resSq
      IF ( nIterMin .GE. 0 ) nIterMin = 0

C--   Relative stopping criterion with safe absolute floor.
      stopSq = MAX( tol*tol*MAX(rhsSq,oneRL), tol*tol )

      printResidual = .FALSE.
      IF ( debugLevel .GE. debLevZero ) THEN
       _BEGIN_MASTER(myThid)
        printResidual = printResidualFreq.GE.1
        WRITE(msgBuf,'(A,1PE21.14,A,1PE21.14,A,1PE21.14)')
     &   ' hashemi_cgadi: init res=', firstResidual,
     &   ' ; rhs=', SQRT(MAX(rhsSq,zeroRL)),
     &   ' ; tol=', SQRT(stopSq)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
       _END_MASTER(myThid)
      ENDIF

      IF ( resSq .LE. stopSq ) THEN
       numIters = 0
       RETURN
      ENDIF

      CALL HASHEMI_APPLY_ADI_PRECOND(
     I      r,
     O      z,
     I      myThid )
      CALL HASHEMI_DOT_RL(
     I      r, z,
     O      rz,
     I      myThid )
      IF ( ABS(rz) .LE. eps ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1,sNy
          DO i=1,sNx
           z(i,j,bi,bj) = r(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       rz = resSq
      ENDIF

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          p(i,j,bi,bj) = z(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      actualIts = 0

      DO it=1,numIters
       CALL EXCH_XY_RL( p, myThid )
       CALL HASHEMI_APPLY_A(
     I        p,
     O        q,
     I        myThid )
       CALL HASHEMI_DOT_RL(
     I        p, q,
     O        pA,
     I        myThid )
       IF ( ABS(pA) .LE. eps ) GOTO 11

       alphaCG = rz/pA
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1,sNy
          DO i=1,sNx
           cg2d_x(i,j,bi,bj) = cg2d_x(i,j,bi,bj) + alphaCG*p(i,j,bi,bj)
           r(i,j,bi,bj) = r(i,j,bi,bj) - alphaCG*q(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO

       CALL HASHEMI_DOT_RL(
     I        r, r,
     O        resSq,
     I        myThid )
       lastResidual = SQRT(MAX(resSq,zeroRL))
       actualIts = it

       IF ( nIterMin .GE. 0 ) THEN
        IF ( resSq .LT. minResidualSq ) THEN
         minResidualSq = resSq
         nIterMin = it
        ENDIF
       ENDIF

       IF ( printResidual ) THEN
        IF ( MOD(it-1,printResidualFreq).EQ.0 ) THEN
         WRITE(msgBuf,'(A,I6,A,1PE21.14)')
     &    ' hashemi_cgadi: iter=', it, ' ; resid.= ', lastResidual
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ENDIF
       ENDIF

       IF ( resSq .LE. stopSq ) GOTO 11

       CALL HASHEMI_APPLY_ADI_PRECOND(
     I        r,
     O        z,
     I        myThid )
       CALL HASHEMI_DOT_RL(
     I        r, z,
     O        rzNew,
     I        myThid )
       IF ( ABS(rz) .LE. eps ) GOTO 11
       betaCG = rzNew/rz
       rz = rzNew

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1,sNy
          DO i=1,sNx
           p(i,j,bi,bj) = z(i,j,bi,bj) + betaCG*p(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

 11   CONTINUE
      numIters = actualIts
      IF ( nIterMin .LT. 0 ) minResidualSq = resSq
      RETURN
      END

C--   !ROUTINE: HASHEMI_APPLY_A
C--   Apply matrix A to a field vector: y = A*x.
      SUBROUTINE HASHEMI_APPLY_A(
     I     xFld,
     O     yFld,
     I     myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "CG2D.h"
      _RL xFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid
      INTEGER bi, bj, i, j

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          yFld(i,j,bi,bj) = zeroRL
         ENDDO
        ENDDO
        DO j=1,sNy
         DO i=1,sNx
          yFld(i,j,bi,bj) =
     &      aW2d(i  ,j  ,bi,bj)*xFld(i-1,j  ,bi,bj)
     &     +aW2d(i+1,j  ,bi,bj)*xFld(i+1,j  ,bi,bj)
     &     +aS2d(i  ,j  ,bi,bj)*xFld(i  ,j-1,bi,bj)
     &     +aS2d(i  ,j+1,bi,bj)*xFld(i  ,j+1,bi,bj)
     &     +aC2d(i  ,j  ,bi,bj)*xFld(i  ,j  ,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END

C--   !ROUTINE: HASHEMI_DOT_RL
C--   Global dot product over interior cell range.
      SUBROUTINE HASHEMI_DOT_RL(
     I     xFld, yFld,
     O     dotVal,
     I     myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
      _RL xFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dotVal
      INTEGER myThid
      INTEGER bi, bj, i, j
      _RL dotTile(nSx,nSy)

      DO bj=1,nSy
       DO bi=1,nSx
        dotTile(bi,bj) = zeroRL
       ENDDO
      ENDDO

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          dotTile(bi,bj) = dotTile(bi,bj)
     &                   + xFld(i,j,bi,bj)*yFld(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL GLOBAL_SUM_TILE_RL( dotTile, dotVal, myThid )
      RETURN
      END

C--   !ROUTINE: HASHEMI_APPLY_ADI_PRECOND
C--   Apply one ADI preconditioning operation: z = M^{-1} r.
C--   Sweep-1: implicit along longitude rows.
C--   Sweep-2: implicit along latitude columns.
      SUBROUTINE HASHEMI_APPLY_ADI_PRECOND(
     I     rFld,
     O     zFld,
     I     myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "CG2D.h"
      _RL rFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid

      INTEGER bi, bj, i, j
      _RL zTmp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL rowL(1:sNx), rowD(1:sNx), rowU(1:sNx)
      _RL rowR(1:sNx), rowX(1:sNx)
      _RL colL(1:sNy), colD(1:sNy), colU(1:sNy)
      _RL colR(1:sNy), colX(1:sNy)
      _RL eps, den

      eps = 1.D-30
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          zTmp(i,j,bi,bj) = zeroRL
          zFld(i,j,bi,bj) = zeroRL
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C--   Row sweep
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          rowL(i) = aW2d(i  ,j,bi,bj)
          rowD(i) = aC2d(i  ,j,bi,bj)
          rowU(i) = aW2d(i+1,j,bi,bj)
          den = rowD(i)
          IF ( ABS(den) .LT. eps ) THEN
           IF ( den .GE. 0.0d0 ) THEN
            rowD(i) = eps
           ELSE
            rowD(i) = -eps
           ENDIF
          ENDIF
          rowR(i) = rFld(i,j,bi,bj)
         ENDDO
         rowR(1) = rowR(1) - rowL(1)*zTmp(0,j,bi,bj)
         rowR(sNx) = rowR(sNx) - rowU(sNx)*zTmp(sNx+1,j,bi,bj)
         rowL(1) = zeroRL
         rowU(sNx) = zeroRL
         CALL HASHEMI_TRIDIAG_THOMAS_X(
     I       rowL, rowD, rowU, rowR,
     O       rowX )
         DO i=1,sNx
          zTmp(i,j,bi,bj) = rowX(i)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL EXCH_XY_RL( zTmp, myThid )

C--   Column sweep
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO i=1,sNx
         DO j=1,sNy
          colL(j) = aS2d(i,j  ,bi,bj)
          colD(j) = aC2d(i,j  ,bi,bj)
          colU(j) = aS2d(i,j+1,bi,bj)
          den = colD(j)
          IF ( ABS(den) .LT. eps ) THEN
           IF ( den .GE. 0.0d0 ) THEN
            colD(j) = eps
           ELSE
            colD(j) = -eps
           ENDIF
          ENDIF
          colR(j) = rFld(i,j,bi,bj)
     &            - aW2d(i  ,j,bi,bj)*zTmp(i-1,j,bi,bj)
     &            - aW2d(i+1,j,bi,bj)*zTmp(i+1,j,bi,bj)
         ENDDO
         colR(1) = colR(1) - colL(1)*zTmp(i,0,bi,bj)
         colR(sNy) = colR(sNy) - colU(sNy)*zTmp(i,sNy+1,bi,bj)
         colL(1) = zeroRL
         colU(sNy) = zeroRL
         CALL HASHEMI_TRIDIAG_THOMAS_Y(
     I       colL, colD, colU, colR,
     O       colX )
         DO j=1,sNy
          zFld(i,j,bi,bj) = colX(j)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL EXCH_XY_RL( zFld, myThid )
      RETURN
      END

C--   !ROUTINE: HASHEMI_SOLVE_ADI_LINE
C--   Alternating-direction line-relaxation solver:
C--   1) implicit in longitude (row tridiagonal solves),
C--   2) implicit in latitude  (column tridiagonal solves).
C--   This follows the same algorithmic idea as Python ADI line sweeps,
C--   and is provided as a solver path for direct comparison.
      SUBROUTINE HASHEMI_SOLVE_ADI_LINE(
     U     cg2d_b, cg2d_x,
     O     firstResidual, lastResidual,
     U     numIters,
     I     myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "CG2D.h"
      _RL cg2d_b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL cg2d_x(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL firstResidual, lastResidual
      INTEGER numIters
      INTEGER myThid

      INTEGER bi, bj, i, j, it, actualIts
      _RL rhsInf, resInf, relRes, tol
      _RL rowL(1:sNx), rowD(1:sNx), rowU(1:sNx)
      _RL rowR(1:sNx), rowX(1:sNx)
      _RL colL(1:sNy), colD(1:sNy), colU(1:sNy)
      _RL colR(1:sNy), colX(1:sNy)
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      LOGICAL printResidual

      tol = cg2dTargetResidual
      IF ( tol .LE. zeroRL ) tol = 1.D-13

      CALL EXCH_XY_RL( cg2d_x, myThid )
      CALL HASHEMI_COMPUTE_RESIDUAL_INF(
     I      cg2d_b, cg2d_x,
     O      resInf, rhsInf,
     I      myThid )
      firstResidual = resInf
      lastResidual = resInf
      actualIts = 0

      printResidual = .FALSE.
      IF ( debugLevel .GE. debLevZero ) THEN
       _BEGIN_MASTER(myThid)
        printResidual = printResidualFreq.GE.1
        WRITE(msgBuf,'(A,1PE21.14,A,1PE21.14)')
     &   ' hashemi_adi: init residual=', resInf, ' ; rhsInf=', rhsInf
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
       _END_MASTER(myThid)
      ENDIF

      relRes = oneRL
      IF ( rhsInf .GT. zeroRL ) relRes = resInf/rhsInf
      IF ( relRes .LE. tol ) THEN
       numIters = 0
       RETURN
      ENDIF

      DO it=1,numIters
C--   Row sweep (longitude implicit, latitude explicit)
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1,sNy
          DO i=1,sNx
           rowL(i) = aW2d(i  ,j,bi,bj)
           rowD(i) = aC2d(i  ,j,bi,bj)
           rowU(i) = aW2d(i+1,j,bi,bj)
           rowR(i) = cg2d_b(i,j,bi,bj)
     &             - aS2d(i,j  ,bi,bj)*cg2d_x(i,j-1,bi,bj)
     &             - aS2d(i,j+1,bi,bj)*cg2d_x(i,j+1,bi,bj)
          ENDDO
          rowR(1) = rowR(1) - rowL(1)*cg2d_x(0     ,j,bi,bj)
          rowR(sNx) = rowR(sNx) - rowU(sNx)*cg2d_x(sNx+1,j,bi,bj)
          rowL(1) = zeroRL
          rowU(sNx) = zeroRL
          CALL HASHEMI_TRIDIAG_THOMAS_X(
     I         rowL, rowD, rowU, rowR,
     O         rowX )
          DO i=1,sNx
           cg2d_x(i,j,bi,bj) = rowX(i)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       CALL EXCH_XY_RL( cg2d_x, myThid )

C--   Column sweep (latitude implicit, longitude explicit)
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO i=1,sNx
          DO j=1,sNy
           colL(j) = aS2d(i,j  ,bi,bj)
           colD(j) = aC2d(i,j  ,bi,bj)
           colU(j) = aS2d(i,j+1,bi,bj)
           colR(j) = cg2d_b(i,j,bi,bj)
     &             - aW2d(i  ,j,bi,bj)*cg2d_x(i-1,j,bi,bj)
     &             - aW2d(i+1,j,bi,bj)*cg2d_x(i+1,j,bi,bj)
          ENDDO
          colR(1) = colR(1) - colL(1)*cg2d_x(i,0     ,bi,bj)
          colR(sNy) = colR(sNy) - colU(sNy)*cg2d_x(i,sNy+1,bi,bj)
          colL(1) = zeroRL
          colU(sNy) = zeroRL
          CALL HASHEMI_TRIDIAG_THOMAS_Y(
     I         colL, colD, colU, colR,
     O         colX )
          DO j=1,sNy
           cg2d_x(i,j,bi,bj) = colX(j)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       CALL EXCH_XY_RL( cg2d_x, myThid )

       CALL HASHEMI_COMPUTE_RESIDUAL_INF(
     I      cg2d_b, cg2d_x,
     O      resInf, rhsInf,
     I      myThid )
       lastResidual = resInf
       actualIts = it
       relRes = oneRL
       IF ( rhsInf .GT. zeroRL ) relRes = resInf/rhsInf

       IF ( printResidual ) THEN
        IF ( MOD(it-1,printResidualFreq).EQ.0 ) THEN
         WRITE(msgBuf,'(A,I6,A,1PE21.14,A,1PE21.14)')
     &    ' hashemi_adi: iter=', it, ' ; abs=', resInf,
     &    ' ; rel=', relRes
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        ENDIF
       ENDIF

       IF ( relRes .LE. tol ) EXIT
      ENDDO

      numIters = actualIts
      RETURN
      END

C--   !ROUTINE: HASHEMI_COMPUTE_RESIDUAL_INF
C--   Compute infinity-norm residual: || b - A*x ||_inf and ||b||_inf.
      SUBROUTINE HASHEMI_COMPUTE_RESIDUAL_INF(
     I     cg2d_b, cg2d_x,
     O     resInf, rhsInf,
     I     myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "CG2D.h"
      _RL cg2d_b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL cg2d_x(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL resInf, rhsInf
      INTEGER myThid
      INTEGER bi, bj, i, j
      _RL ax, rr, babs

      resInf = zeroRL
      rhsInf = zeroRL
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          ax = aW2d(i  ,j  ,bi,bj)*cg2d_x(i-1,j  ,bi,bj)
     &       +aW2d(i+1,j  ,bi,bj)*cg2d_x(i+1,j  ,bi,bj)
     &       +aS2d(i  ,j  ,bi,bj)*cg2d_x(i  ,j-1,bi,bj)
     &       +aS2d(i  ,j+1,bi,bj)*cg2d_x(i  ,j+1,bi,bj)
     &       +aC2d(i  ,j  ,bi,bj)*cg2d_x(i  ,j  ,bi,bj)
          rr = ABS( cg2d_b(i,j,bi,bj) - ax )
          babs = ABS( cg2d_b(i,j,bi,bj) )
          IF ( rr .GT. resInf ) resInf = rr
          IF ( babs .GT. rhsInf ) rhsInf = babs
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      _GLOBAL_MAX_RL( resInf, myThid )
      _GLOBAL_MAX_RL( rhsInf, myThid )
      RETURN
      END

C--   !ROUTINE: HASHEMI_TRIDIAG_THOMAS_X
C--   Thomas algorithm for X-line solves (length sNx).
      SUBROUTINE HASHEMI_TRIDIAG_THOMAS_X(
     I     aL, aD, aU, rhs,
     O     sol )
      IMPLICIT NONE
#include "SIZE.h"
      _RL aL(1:sNx), aD(1:sNx), aU(1:sNx)
      _RL rhs(1:sNx), sol(1:sNx)
      INTEGER i
      _RL cP(1:sNx), dP(1:sNx), den, eps

      eps = 1.D-30
      den = aD(1)
      IF ( ABS(den) .LT. eps ) THEN
       IF ( den .GE. 0.0d0 ) THEN
        den = eps
       ELSE
        den = -eps
       ENDIF
      ENDIF
      cP(1) = aU(1)/den
      dP(1) = rhs(1)/den

      DO i=2,sNx
       den = aD(i) - aL(i)*cP(i-1)
       IF ( ABS(den) .LT. eps ) THEN
        IF ( den .GE. 0.0d0 ) THEN
         den = eps
        ELSE
         den = -eps
        ENDIF
       ENDIF
       IF ( i .LT. sNx ) cP(i) = aU(i)/den
       dP(i) = ( rhs(i) - aL(i)*dP(i-1) )/den
      ENDDO

      sol(sNx) = dP(sNx)
      DO i=sNx-1,1,-1
       sol(i) = dP(i) - cP(i)*sol(i+1)
      ENDDO
      RETURN
      END

C--   !ROUTINE: HASHEMI_TRIDIAG_THOMAS_Y
C--   Thomas algorithm for Y-line solves (length sNy).
      SUBROUTINE HASHEMI_TRIDIAG_THOMAS_Y(
     I     aL, aD, aU, rhs,
     O     sol )
      IMPLICIT NONE
#include "SIZE.h"
      _RL aL(1:sNy), aD(1:sNy), aU(1:sNy)
      _RL rhs(1:sNy), sol(1:sNy)
      INTEGER j
      _RL cP(1:sNy), dP(1:sNy), den, eps

      eps = 1.D-30
      den = aD(1)
      IF ( ABS(den) .LT. eps ) THEN
       IF ( den .GE. 0.0d0 ) THEN
        den = eps
       ELSE
        den = -eps
       ENDIF
      ENDIF
      cP(1) = aU(1)/den
      dP(1) = rhs(1)/den

      DO j=2,sNy
       den = aD(j) - aL(j)*cP(j-1)
       IF ( ABS(den) .LT. eps ) THEN
        IF ( den .GE. 0.0d0 ) THEN
         den = eps
        ELSE
         den = -eps
        ENDIF
       ENDIF
       IF ( j .LT. sNy ) cP(j) = aU(j)/den
       dP(j) = ( rhs(j) - aL(j)*dP(j-1) )/den
      ENDDO

      sol(sNy) = dP(sNy)
      DO j=sNy-1,1,-1
       sol(j) = dP(j) - cP(j)*sol(j+1)
      ENDDO
      RETURN
      END

C--   !ROUTINE: HASHEMI_BUILD_CG2D_PRECOND_ADI
C--   ADI-inspired factorized sparse preconditioner compatible with CG2D.
C--   (CG2D applies preconditioning as a 5-point stencil via pC/pW/pS.)
      SUBROUTINE HASHEMI_BUILD_CG2D_PRECOND_ADI( myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "CG2D.h"
      INTEGER myThid, bi, bj, i, j
      _RL tol
      _RL aC, aCw, aCs
      _RL aWloc, aEloc, aSloc, aNloc
      _RL diagAbs, dirX, dirY
      _RL denW, denS, tauX, tauY
      _RL adiTheta, adiShift, epsPre

C--   Reset preconditioner arrays
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          pC(i,j,bi,bj) = zeroRL
          pW(i,j,bi,bj) = zeroRL
          pS(i,j,bi,bj) = zeroRL
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      _BEGIN_MASTER(myThid)
      cg2dNorm = oneRL
      cg2dNormaliseRHS = .FALSE.
      tol = cg2dTargetResidual
      IF ( tol .LE. zeroRL ) tol = 1.D-13
      cg2dTolerance_sq = tol*tol
      _END_MASTER(myThid)

C--   Blend weight in [0.25,0.85] for robust directional damping.
      adiTheta = cg2dpcOffDFac
      IF ( adiTheta .LT. 0.25d0 ) adiTheta = 0.25d0
      IF ( adiTheta .GT. 0.85d0 ) adiTheta = 0.85d0
C--   Small shift stabilises nearly singular wet cells.
      adiShift = 0.10d0
      epsPre = 1.D-30

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          aC = aC2d(i,j,bi,bj)
          IF ( ABS(aC) .LE. epsPre ) THEN
           pC(i,j,bi,bj) = oneRL
           pW(i,j,bi,bj) = zeroRL
           pS(i,j,bi,bj) = zeroRL
          ELSE
           aCw = aC2d(i-1,j,bi,bj)
           aCs = aC2d(i,j-1,bi,bj)
           aWloc = aW2d(i,j,bi,bj)
           aEloc = aW2d(i+1,j,bi,bj)
           aSloc = aS2d(i,j,bi,bj)
           aNloc = aS2d(i,j+1,bi,bj)

           pC(i,j,bi,bj) = oneRL / aC

           diagAbs = ABS(aC)
           dirX = ABS(aWloc) + ABS(aEloc)
           dirY = ABS(aSloc) + ABS(aNloc)
           tauX = adiShift*MAX(diagAbs,dirX)
           tauY = adiShift*MAX(diagAbs,dirY)

           denW = adiTheta*diagAbs
     &          + (oneRL-adiTheta)*ABS(aCw) + tauX
           denS = adiTheta*diagAbs
     &          + (oneRL-adiTheta)*ABS(aCs) + tauY

           IF ( denW .GT. epsPre ) THEN
            pW(i,j,bi,bj) = -aWloc / ( denW*diagAbs )
           ELSE
            pW(i,j,bi,bj) = zeroRL
           ENDIF

           IF ( denS .GT. epsPre ) THEN
            pS(i,j,bi,bj) = -aSloc / ( denS*diagAbs )
           ELSE
            pS(i,j,bi,bj) = zeroRL
           ENDIF
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL EXCH_XY_RS( pC, myThid )
      CALL EXCH_UV_XY_RS( pW, pS, .FALSE., myThid )
      RETURN
      END


C-- ROCKET: sacing subroutines                                           

C--   !ROUTINE: HASHEMI_WRITE_UV_CSV
C--   Write diagnosed k=1 uVel/vVel fields to per-rank CSV for plotting.
      SUBROUTINE HASHEMI_WRITE_UV_CSV( myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"

      INTEGER myThid
      INTEGER bi, bj, i, j, iLat, iLonG, ioUnit
      _RL uLoc, vLoc, speed
      CHARACTER*(MAX_LEN_FNAM) outFile

      WRITE(outFile,'(A,I8.8,A)')
     & '../results/solutionVel/VelOcean.', myProcId, '.csv'

      ioUnit = 170 + myProcId
      OPEN(UNIT=ioUnit,FILE=outFile,STATUS='REPLACE')
      WRITE(ioUnit,'(A)') 'iLat,iLon,uVel,vVel,speed'

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         iLat = myYGlobalLo + (bj-1)*sNy + j - 1
         DO i=1,sNx
          iLonG = myXGlobalLo + (bi-1)*sNx + i - 1
          uLoc = uVel(i,j,1,bi,bj)
          vLoc = vVel(i,j,1,bi,bj)
          speed = SQRT( uLoc*uLoc + vLoc*vLoc )
          WRITE(ioUnit,
     & '(I0,",",I0,",",ES24.16,",",ES24.16,",",ES24.16)')
     &     iLat, iLonG, uLoc, vLoc, speed
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CLOSE(ioUnit)
      RETURN
      END

C--   !ROUTINE: HASHEMI_WRITE_PSI_CSV
      SUBROUTINE HASHEMI_WRITE_PSI_CSV( cg2d_x, myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"

      INTEGER myThid
      _RL cg2d_x(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER bi, bj, i, j, iLat, iLonG, ioUnit
      CHARACTER*(MAX_LEN_FNAM) outFile

C--   Every MPI process creates its own filename based on myProcId
      WRITE(outFile,'(A,I8.8,A)')
     & '../results/solutionPsi/PsiOcean.', myProcId, '.csv'
      
      ioUnit = 70 + myProcId
      OPEN(UNIT=ioUnit,FILE=outFile,STATUS='REPLACE')

C--   Write Header to match your old format
      WRITE(ioUnit,'(A)') 'iLat,iLon,psi'

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         ! Calculate Global Latitude (1 to 720)
         iLat = myYGlobalLo + (bj-1)*sNy + j - 1
         DO i=1,sNx
          ! Calculate Global Longitude (1 to 1440)
          iLonG = myXGlobalLo + (bi-1)*sNx + i - 1
          
          ! Write row exactly like your old code
          WRITE(ioUnit,'(I0,",",I0,",",ES24.16)')
     &      iLat, iLonG, cg2d_x(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CLOSE(ioUnit)
      RETURN
      END
