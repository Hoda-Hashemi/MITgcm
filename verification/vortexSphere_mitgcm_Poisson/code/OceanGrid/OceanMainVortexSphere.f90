program HodamainVortexSphere
	use simulationParameters
	use OceanGridDataModule
	use OceanAdaptiveGridDataModule
	use coriolisTermDataModule
	use outputModule
	use initializeModule
	use gridOperationsModule
	use poissonSolverModule
	use basicDataStructures
	use AnalysisModule
	! use petscModule
    !$ use omp_lib
	implicit none

	integer :: i, kkk
	integer :: width
	!COMMENT: fancy colors in the terminal
	character(len=*), parameter :: reset = char(27)//'[0m'
	character(len=*), parameter :: green  = char(27)//"[38;2;18;237;58m"  ! #12ED3A
	character(len=*), parameter :: orange  = char(27)//"[38;2;255;202;0m" ! #FFCA00
	character(len=*), parameter :: cyan  = char(27)//"[38;2;6;249;242m"   ! #06F9F2 
	character(len=*), parameter :: red  = char(27)//"[38;2;249;6;13m"     ! #F9060D
	character(len=*), parameter :: neon5  = char(27)//"[38;2;243;12;152m" ! #F30C98
	character(len=*), parameter :: pink  = char(27)//"[38;2;223;9;246m"   ! #DF09F6 
	character(len=*), parameter :: lime  = char(27)//"[38;2;175;236;19m"  ! #AFEC13  
	character(len=*), parameter :: sky = char(27)//"[38;2;0;142;255m"     ! #008EFF  
    
	call system('uname > ../results/os.txt')
	open(UNIT=2,FILE='../results/os.txt')
	read(2,*)osName
	close(2)
	width = 142 
	
    !COMMENT: Constructing the Grid Over the Earth
	print '(A)', sky // repeat('─', width) // reset
	print '(A)', sky // ' Ocean: Constructing a grid (nθ x nƛ)' //reset
	print '(A)', sky // repeat('─', width) // reset
	call constructOceanGrid()
	! call constructGrid()
	write(*,*)
	! stop

    !COMMENT: Reading Bathymetry and Interpolating to GRID
	if(include_bathymetry)then
		print '(A)', orange // repeat('─', width) // reset
		print '(A)', orange // ' Reading Bathymetry Data etopo5 from NOAA ' //reset
		print '(A)', orange // repeat('─', width) // reset
		call readBathymetry()

		if(save_bathymetry) then
			print '(A)', lime // repeat('─', width) // reset
			print '(A)', lime // ' Ocean: Saving Bathymetry into file '  //reset
			print '(A)', lime // repeat('─', width) // reset
			call saveBathymetry()
		end if

		print '(A)', lime // repeat('─', width) // reset
		print '(A)', lime // ' Ocean: Interpolate Bathymetry To Ocean Grid '  //reset
		print '(A)', lime // repeat('─', width) // reset
		call interpolateBathymetryToOceanGrid()
		! call interpolateBathymetryToGrid()
	end if
	write(*,*) 
	! stop

    !COMMENT: Building Psi One Coefficients Matrix ψ₁
	print '(A)', lime // repeat('─', width) // reset
	print '(A)', lime // ' Ocean: Builing Coefficient Matrix for ψ₁:   ∇ψ₁ - (1/Ld²) ψ₁ = - f(θ)  ' //reset
	print '(A)', lime // repeat('─', width) // reset
	! call OceanbuildCoefficientsMatrixPsiOne()
	! stop
	! call buildCoefficientsMatrix()
	write(*,*)

	!COMMENT: Solving for ψ₁
	print '(A)', lime // repeat('─', width) // reset
	print '(A)', lime // ' Ocean: Solving for ψ₁ ' //reset
	print '(A)', lime // repeat('─', width) // reset
	! call solvePoissonForPsiOne()
	! call solvePoissonForPsi()
	write(*,*)

	call ocean_solve_poisson_cg2d(1)
	
	!COMMENT: Computing Velocity field for uψ₁ = uψ₁(θ,ƛ)
	! print '(A)', lime // repeat('─', width) // reset
	! print '(A)', lime // ' Ocean: Calculating u₁(θ,ƛ) ' //reset
	! print '(A)', lime // repeat('─', width) // reset	
	! call VelocityFieldsOneOnGrid()

	!COMMENT: Saving the MESH
	print '(A)', orange // repeat('─', width) // reset
	print '(A)', orange // ' Ocean: Saving the Mesh '//reset
	print '(A)', orange // repeat('─', width) // reset
	call saveOceanMesh()
	! kkk = 1
	! call saveMesh(kkk)
    write(*,*)

	!COMMENT: Done statement
	print '(A)', green // repeat('─', width) // reset
	print '(A)', green // ' Ocean: Task Completed Successfully ' // reset
	print '(A)', green // repeat('─', width) // reset

	! deallocate(psi1onGrid)
	! deallocate(velocityPsiOneOnGrid)

end program HodamainVortexSphere

