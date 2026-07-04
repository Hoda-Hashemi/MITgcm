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

      PARAMETER ( sNx =   16 )        ! CS32 tile width
      PARAMETER ( sNy =    8 )        ! CS32 tile height
      PARAMETER ( OLx =    2 )        ! overlap in x-dir
      PARAMETER ( OLy =    2 )        ! overlap in y-dir
      PARAMETER ( nSx =    1 )
      PARAMETER ( nSy =    1 )
      PARAMETER ( nPx =   48 )
      PARAMETER ( nPy =    1 )
      PARAMETER ( Nx  = sNx*nSx*nPx )
      PARAMETER ( Ny  = sNy*nSy*nPy )
      PARAMETER ( Nr  =    1 )        ! vertical levels

      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,
     &            MAX_OLY = OLy )
