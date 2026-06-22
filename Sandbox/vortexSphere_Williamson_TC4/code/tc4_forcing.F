#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: TC4_FORCING_POINT
C     !INTERFACE:
      SUBROUTINE TC4_FORCING_POINT(
     I     lonDeg, latDeg, myTime,
     O     forcU, forcV, forcH )
C     !DESCRIPTION:
C     Williamson TC4 analytic source terms at a single spherical point.
C     lonDeg and latDeg are MITgcm grid coordinates in degrees.
CEOP
      IMPLICIT NONE

      _RL lonDeg, latDeg, myTime
      _RL forcU, forcV, forcH

      _RL pi, deg2radLoc
      _RL radius, omega, grav
      _RL gh0, lat0, lon0, sigma, f0, psi0
      PARAMETER ( pi = 3.14159265358979323844D0 )
      PARAMETER ( deg2radLoc = pi/180. _d 0 )
      PARAMETER ( radius = 6371000. _d 0 )
      PARAMETER ( omega = 7.292 _d -5 )
      PARAMETER ( grav = 9.81 _d 0 )
      PARAMETER ( gh0 = 1. _d 5 )
      PARAMETER ( lat0 = pi/4. _d 0 )
      PARAMETER ( lon0 = -7. _d 0*pi/12. _d 0 )
      PARAMETER ( sigma = 12.74244 _d 0*12.74244 _d 0 )
      PARAMETER ( f0 = 2. _d 0*omega*SIN(pi/4. _d 0) )
      PARAMETER ( psi0 = -0.03 _d 0*gh0/f0 )

      _RL lon, lat, sinLat, cosLat, tanLat, coriolis
      _RL sin2, cos2, uBar, duBar, ghBarLat
      _RL psi, psiLon, psiLat, psiLL, psiLT, psiTT
      _RL u0, phaseSpeed
      _RL uVel, vVel, uLon, uLat, uTime
      _RL vLon, vLat, vTime, ghLon, ghLat, ghTime
      _RL totalDU, totalDV

      CALL TC4_GET_U0( u0 )
      phaseSpeed = u0/radius

      lon = lonDeg*deg2radLoc
      lat = latDeg*deg2radLoc
      sinLat = SIN(lat)
      cosLat = COS(lat)

      IF ( ABS(cosLat) .LT. 1. _d -12 ) THEN
       forcU = 0. _d 0
       forcV = 0. _d 0
       forcH = 0. _d 0
       RETURN
      ENDIF

      tanLat = sinLat/cosLat
      coriolis = 2. _d 0*omega*sinLat

      sin2 = SIN(2. _d 0*lat)
      cos2 = COS(2. _d 0*lat)
      uBar = u0*sin2**14
      duBar = 28. _d 0*u0*sin2**13*cos2
      ghBarLat = -(radius*coriolis + uBar*tanLat)*uBar

      CALL TC4_PSI_DERIVS(
     I     lon, lat, myTime, phaseSpeed, lon0, lat0,
     I     psi0, sigma,
     O     psi, psiLon, psiLat, psiLL, psiLT, psiTT )

      uVel = uBar - psiLat/radius
      vVel = psiLon/(radius*cosLat)

      uLon = -psiLT/radius
      uLat = duBar - psiTT/radius
      uTime = phaseSpeed*psiLT/radius

      vLon = psiLL/(radius*cosLat)
      vLat = psiLT/(radius*cosLat)
     &     + psiLon*sinLat/(radius*cosLat*cosLat)
      vTime = -phaseSpeed*psiLL/(radius*cosLat)

      ghLon = f0*psiLon
      ghLat = ghBarLat + f0*psiLat
      ghTime = -f0*phaseSpeed*psiLon

      totalDU = uTime
     &        + uVel*uLon/(radius*cosLat)
     &        + vVel*uLat/radius
      totalDV = vTime
     &        + uVel*vLon/(radius*cosLat)
     &        + vVel*vLat/radius

      forcU = totalDU - uVel*vVel*tanLat/radius
     &      + ghLon/(radius*cosLat) - coriolis*vVel
      forcV = totalDV + uVel*uVel*tanLat/radius
     &      + ghLat/radius + coriolis*uVel
      forcH = ( ghTime + uVel*ghLon/(radius*cosLat)
     &        + vVel*ghLat/radius )/grav

      RETURN
      END

CBOP
C     !ROUTINE: TC4_GET_U0
C     !INTERFACE:
      SUBROUTINE TC4_GET_U0( u0 )
C     !DESCRIPTION:
C     Runtime TC4 background jet speed. Defaults to the Williamson c1 case.
CEOP
      IMPLICIT NONE

      _RL u0
      _RL savedU0
      CHARACTER*32 env
      LOGICAL initialized
      SAVE savedU0, initialized
      DATA initialized /.FALSE./

      IF ( .NOT. initialized ) THEN
       savedU0 = 20. _d 0
       env = ' '
       CALL GETENV( 'TC4_U0_VALUE', env )
       IF ( env .NE. ' ' ) THEN
        READ(env,*,ERR=10) savedU0
       ENDIF
 10    CONTINUE
       initialized = .TRUE.
      ENDIF

      u0 = savedU0

      RETURN
      END

CBOP
C     !ROUTINE: TC4_PSI_DERIVS
C     !INTERFACE:
      SUBROUTINE TC4_PSI_DERIVS(
     I     lon, lat, myTime, phaseSpeed, lon0, lat0, psi0, sigma,
     O     psi, psiLon, psiLat, psiLL, psiLT, psiTT )
C     !DESCRIPTION:
C     First and second derivatives of the TC4 translating-low streamfunction.
CEOP
      IMPLICIT NONE

      _RL lon, lat, myTime, phaseSpeed, lon0, lat0, psi0, sigma
      _RL psi, psiLon, psiLat, psiLL, psiLT, psiTT

      _RL phase, sinLat, cosLat, sin0, cos0
      _RL cosPhase, sinPhase, cVal, denom, scale, secondScale
      _RL cLon, cLat, cLL, cLT, cTT

      phase = lon - phaseSpeed*myTime - lon0
      sinLat = SIN(lat)
      cosLat = COS(lat)
      sin0 = SIN(lat0)
      cos0 = COS(lat0)
      cosPhase = COS(phase)
      sinPhase = SIN(phase)

      cVal = sin0*sinLat + cos0*cosLat*cosPhase
      denom = 1. _d 0 + cVal
      IF ( denom .LT. 1. _d -12 ) THEN
       psi = 0. _d 0
       psiLon = 0. _d 0
       psiLat = 0. _d 0
       psiLL = 0. _d 0
       psiLT = 0. _d 0
       psiTT = 0. _d 0
       RETURN
      ENDIF

      psi = psi0*EXP(-sigma*(1. _d 0-cVal)/denom)
      scale = 2. _d 0*sigma/(denom*denom)
      secondScale = scale*scale - 2. _d 0*scale/denom

      cLon = -cos0*cosLat*sinPhase
      cLat = sin0*cosLat - cos0*sinLat*cosPhase
      cLL = -cos0*cosLat*cosPhase
      cLT = cos0*sinLat*sinPhase
      cTT = -sin0*sinLat - cos0*cosLat*cosPhase

      psiLon = psi*scale*cLon
      psiLat = psi*scale*cLat
      psiLL = psi*(scale*cLL + cLon*cLon*secondScale)
      psiLT = psi*(scale*cLT + cLon*cLat*secondScale)
      psiTT = psi*(scale*cTT + cLat*cLat*secondScale)

      RETURN
      END
