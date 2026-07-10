C    !===  GRID SPACING & DOMAIN SIZE =====================================
      INTEGER sNx
      INTEGER sNy
      INTEGER OLx
      INTEGER OLy
      INTEGER nSx
      INTEGER nSy
      INTEGER nPx
      INTEGER nPy
      INTEGER Nx
      INTEGER Ny
      INTEGER Nr

      PARAMETER ( sNx =  240 )        ! 1440 longitude cells / 6 ranks
      PARAMETER ( sNy =  144 )        ! 720 latitude cells / 5 ranks
      PARAMETER ( OLx =    2 )        ! overlap in x-dir
      PARAMETER ( OLy =    2 )        ! overlap in y-dir
      PARAMETER ( nSx =    1 )
      PARAMETER ( nSy =    1 )
      PARAMETER ( nPx =    6 )
      PARAMETER ( nPy =    5 )
      PARAMETER ( Nx  = sNx*nSx*nPx )
      PARAMETER ( Ny  = sNy*nSy*nPy )
      PARAMETER ( Nr  =    1 )        ! vertical levels

      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,
     &            MAX_OLY = OLy )
