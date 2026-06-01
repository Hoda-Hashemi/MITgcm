#include "MYPACKAGE_OPTIONS.h"

CBOP
C !ROUTINE: MYPACKAGE_CALC_RHS

C !INTERFACE: ==========================================================
      SUBROUTINE MYPACKAGE_CALC_RHS(
     I                     myTime, myIter, myThid )

C !DESCRIPTION:
C Prescribe the TC1 velocity field every time step. The flow is frozen
C in time here; only the free surface is advanced by continuity.

C !USES: ===============================================================
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "MYPACKAGE.h"

C !INPUT PARAMETERS: ===================================================
C myTime    :: Current time in simulation
C myIter    :: Current time-step number in simulation
C myThid    :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

#ifdef ALLOW_MYPACKAGE
C !LOCAL VARIABLES: ====================================================
      INTEGER bi, bj
CEOP

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
         CALL TC1_PRESCRIBE_VELOCITY( bi, bj, myThid )
       ENDDO
      ENDDO

#endif /* ALLOW_MYPACKAGE */

      RETURN
      END
