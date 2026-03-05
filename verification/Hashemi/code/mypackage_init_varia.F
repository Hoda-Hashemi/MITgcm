#include "MYPACKAGE_OPTIONS.h"

      SUBROUTINE MYPACKAGE_INIT_VARIA( myThid )
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     === Routine arguments ===
      INTEGER myThid

C     === Local variables ===
C     Use the official constants from EEPARAMS.h instead of 
C     re-defining them, to avoid the syntax errors you see above.
      
      LOGICAL startupSolved
      SAVE startupSolved
      DATA startupSolved /.FALSE./

      IF ( startupSolved ) RETURN

      CALL HASHEMI_CG2D_STARTUP_SOLVE( myThid )
      startupSolved = .TRUE.

      RETURN
      END

