module basicDataStructures
	implicit none
	save
	integer, parameter :: QR_K = selected_real_kind (32)
	integer, parameter :: QI_K = selected_int_kind (16)
	real (kind=QR_K), parameter :: pi=real(4,QR_K)*atan(real(1.0,QR_K))
	real (kind=QR_K), parameter :: half=real(0.5,QR_K)
	real (kind=QR_K), parameter :: quarter=real(0.25,QR_K)
	real (kind=QR_K), parameter :: zero=real(0,QR_K)
	real (kind=QR_K), parameter :: one=real(1,QR_K)
	real (kind=QR_K), parameter :: two=real(2,QR_K)
	real (kind=QR_K), parameter :: three=real(3,QR_K)
	real (kind=QR_K), parameter :: four=real(4,QR_K)
	real (kind=QR_K), parameter :: five=real(5,QR_K)
	real (kind=QR_K), parameter :: eight=real(8,QR_K)
	real (kind=QR_K), parameter :: ten=real(10,QR_K)
	real (kind=QR_K), parameter :: twelve=real(12,QR_K)
	real (kind=QR_K), parameter :: twenty=real(20,QR_K)
	real (kind=QR_K), parameter :: bandie=real(24,QR_K)

	real (kind=QR_K), parameter :: hundred=real(100,QR_K)

	real (kind=QR_K), parameter :: thousand=real(1000,QR_K)
	real (kind=QR_K), parameter :: twothousand=real(2000,QR_K)
	real (kind=QR_K), parameter :: fivethousand=real(5000,QR_K)

	real (kind=QR_K), parameter :: tenthousand=real(10000,QR_K)
	real (kind=QR_K), parameter :: hundredthousand=real(100000,QR_K)
	real (kind=QR_K), parameter :: million=real(1000000,QR_K)
end module basicDataStructures

module simulationParameters
	use basicDataStructures
	implicit none
	save
	!ALTERED: free surface is set to TRUE now. !.TRUE. !!.FALSE.                                     
	logical, parameter :: include_hs = .TRUE. 
	real (kind=QR_K), parameter :: R = real(6371000.0,QR_K) !meters
!	real (kind=QR_K) :: Ld = real(1000000.0,QR_K) ! meters ;!Ld = real(1000000.0,QR_K) ! meters
	real (kind=QR_K), parameter :: OMEGA=two*pi/real(24*3600,QR_K)

	real (kind=QR_K), parameter :: GRAVITY=real(9.81,QR_K)

	!ADDED                                                          
	real (kind=QR_K), parameter :: f2OMEGA= two*OMEGA
	real (kind=QR_K), parameter :: d0Constant = thousand
	real(kind=QR_K), parameter :: sinTheta02 = sin(pi*quarter)* sin(pi*quarter)

	CHARACTER(len=256) :: osName
!	real (kind=QR_K), parameter :: Ld = real(1000000.0,QR_K) ! meters ;!Ld = real(1000000.0,QR_K) ! meters
!	real (kind=QR_K), parameter :: a_sim = real(4.1D7,QR_K) ! a is in m^2/s

	!ALTERED:  ! huge(0.0) !R !R/ten  !R/hundred !R/thousand   !thousand!hundred !real(1000000.0,QR_K) ! meters ;!Ld = real(1000000.0,QR_K) ! meters
	! real (kind=QR_K) :: Ld0 = thousand
	 
	!ADDED: I want to evaluate it and not take it as constant
	real (kind=QR_K) :: Ld0 = hundred
	! sqrt(GRAVITY*D)/f2OMEGA
	
	real (kind=QR_K), parameter :: a_sim = real(1.0D6,QR_K) ! a is in m^2/s
	integer, parameter :: m_sim=4
	real (kind=QR_K), parameter ::  gammaCutoffFactor=real(10,QR_K)
	real (kind=QR_K), parameter ::  deltat = 1.0D0*3600.0D0;

	integer, parameter :: regridFrequency = 5
	integer, parameter :: saveMeshFrequency = 1000000
	integer, parameter :: saveParticlesFrequency = 1

	!ALTERED: include_bathymetry is set to FALSE now. !.TRUE. !.FALSE.                                  
	logical, parameter :: include_bathymetry =.TRUE.
	logical, parameter :: save_bathymetry =.TRUE. 
	!ADDED: SCALING BOTH SIDES OR NOT
	logical, parameter :: ScaleEquation = .TRUE. 
	! logical, parameter :: UseNewSolver = .FALSE.

end module simulationParameters

module greenFunctionDataModule
	use basicDataStructures
	implicit none
	save
	type blockDataStructure
		integer (kind=QR_K) :: iStartZone, iEndZone, nBlockDivisions
		integer*1, allocatable :: iZone(:)
	end type blockDataStructure
real (kind=QR_K) :: cosGammaMAX
real (kind=QR_K) :: cosGammaMIN
	type greenFunctionDataStructure
!		real (kind=QR_K) :: dGamma
		real (kind=QR_K) :: dCosGamma
		integer :: nZones
		integer (kind=QR_K), allocatable :: nCellsZone(:), startingIndex(:)
		real (kind=QR_K), allocatable :: greenTable(:),greenTableA(:),dCosGammaA(:),cosGamma1(:),cosGamma2(:)
		real (kind=QR_K), allocatable :: greenVTable(:) !, cosGammaArrFull(:)
		integer :: greenTableNEntries
		integer :: nBlocks
		type(blockDataStructure), allocatable :: blocks(:)
	end type greenFunctionDataStructure
	type(greenFunctionDataStructure) :: greenFunctionData
	integer,parameter :: ndGammaPerGammaStar = 100000 !used in constructing the Green function table
end module greenFunctionDataModule

module coriolisTermDataModule
	use basicDataStructures
	implicit none
	save
	type coriolisTermDataStructure
		real (kind=QR_K) :: dtheta
		real (kind=QR_K), allocatable :: coriolisTermTable(:,:)
		integer :: ntheta = 360*4
	end type coriolisTermDataStructure
	integer :: nPhiQuad = 8
	type(coriolisTermDataStructure) :: coriolisTermData
end module coriolisTermDataModule

module sphericalGridDataModule
	use basicDataStructures
	implicit none
	save
	type sphericalGridDataStructure
		integer :: nLat =  720 !2160
		integer :: nPhi
		real (kind=QR_K) :: dLat, dphi, d0
		real (kind=QR_K), allocatable, dimension(:) :: sinLat,cosLat,sinPhi,cosPhi, dA
		real (kind=QR_K), allocatable, dimension(:) :: sinLatC,cosLatC,sinPhiC,cosPhiC
		real (kind=QR_K), allocatable, dimension(:,:) :: d
	end type sphericalGridDataStructure
	type(sphericalGridDataStructure) :: sphericalGrid
	real (kind=QR_K), allocatable :: QOnGrid(:,:)
	real (kind=QR_K), allocatable :: hsPrevOnGrid(:,:)
	real (kind=QR_K), allocatable :: divOnGrid(:,:)
	real (kind=QR_K), allocatable :: psiOnGrid(:,:)
	real (kind=QR_K), allocatable :: phiOnGrid(:,:)
	real (kind=QR_K), allocatable :: velocityOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityVorOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityDivOnGrid(:,:,:)
	integer :: nActiveCells
	LOGICAL, allocatable :: isActiveCell(:,:)
	integer (KIND=2), allocatable :: activeCells(:,:)
end module sphericalGridDataModule

module bathymetryDataModule
	implicit none
	save
	type bathymetryDataStructure
		integer :: nLong=4320
		integer :: nLat=2161
		INTEGER (KIND=2), allocatable :: z(:,:)
	end type bathymetryDataStructure
	type(bathymetryDataStructure) :: bathymetry
end module bathymetryDataModule

module particlesDataModule
	use basicDataStructures
	implicit none
	save
	integer, parameter :: nParticlesMAX=10000000
	integer :: nParticles
	type particle
		real (kind=QR_K) :: lat, phi, Q, strengthVor
		real (kind=QR_K) :: sinLat,cosLat,sinPhi,cosPhi
		real (kind=QR_K) :: uLat,uLon
		real (kind=QR_K) :: h,d,hPrev,hs,hsPrev,strengthDiv
		integer :: index, indexZone
	end type particle

	type elementShell
		type(particle), allocatable :: el(:)
	end type elementShell
	
	type(elementShell), allocatable :: particles(:)
end module particlesDataModule

module sparseMatricesDataModule
implicit none
save
!double precision :: nmaxOvernElements = 0
type sprsmatrix
	integer*4, allocatable :: i(:)
	integer*4, allocatable :: j(:)
	integer*8, allocatable :: nentries(:)
	double precision, allocatable :: value(:)
	integer :: m,n,rowIndexed
	integer*8 :: length
end type sprsmatrix

type sprsmatrixA
	integer*4, allocatable :: j(:)
	integer*8, allocatable :: nentries(:)
	REAL*4, allocatable :: value(:)
	! double precision, allocatable :: value(:)
	integer :: m,n,rowIndexed
	integer*8 :: length
end type sprsmatrixA

type sprsrow
	integer :: nentries
	integer, allocatable :: j(:)
	double precision, allocatable :: value(:)
end type sprsrow

type sprsmatrixnew
	type(sprsrow), allocatable :: sprsrows(:)
	integer :: length
end type sprsmatrixnew

type(sprsmatrixA) :: AC
type(sprsmatrixA) :: BC

end module sparseMatricesDataModule

module GMRESSettings
implicit none
save
integer,parameter :: GMRES_ninner = 20
integer,parameter :: GMRES_nouter = 1000
double precision :: GMRES_tol=1.0D-8
double precision,parameter :: GMRES_mu=0.0D0
integer,parameter :: GMRES_constrained=1
integer,parameter :: GMRES_save = 0
double precision,parameter :: GMRES_threshold=5.0E-008
end module GMRESSettings

module QuickSearchVariables
	use basicDataStructures
	implicit none
	save
	real(kind=QR_K) :: phiMinQS,phiMaxQS ! The X boundaries for the box that contains all the vorticies
	real(kind=QR_K) :: latMinQS,latMaxQS ! The Y boundaries for the box that contains all the vorticies

	real(kind=QR_K) :: zoneSize  ! The width size of the unit zone

	integer :: nPhiZones            ! The number of unit zones along X direction
	integer :: nLatZones            ! The number of unit zones along Y direction
	!integer :: nBlck                ! The number of unit zones in the (YZ) Plane
	integer :: nZones              ! The total number of unit zones

	type zone
		integer :: nEls             ! The number of elements contained in each unit zone
		integer :: startLoc       ! The starting location in Renumbering the elements contained in each unit zone
	end type zone

	integer,dimension (9) :: neighZoneRelInd        ! Used in mapping from a unit zone to its 27 neighbours
	type(zone), allocatable, dimension(:) :: zones     ! Unit Zones

end module QuickSearchVariables

module quadIntegrationDataModule
	use basicDataStructures
implicit none
save
real(kind=QR_K) :: weight9(9), dx9(9), dy9(9)
real(kind=QR_K) :: weight4(4), dx4(4), dy4(4)
real(kind=QR_K) :: weight1(1) !, dx1, dy1

end module quadIntegrationDataModule

module smoothedGreenDataModule
	use basicDataStructures
	implicit none
	save
	real (kind=QR_K) :: sigmaFraction=one/ten
	type smoothedGreenDataStructure
		real (kind=QR_K) :: dtheta
		real (kind=QR_K), allocatable :: smoothedGreenTermTable(:,:)
		integer :: ntheta = 360*4
	end type smoothedGreenDataStructure
	integer :: nPhiQuad = 8
	type(smoothedGreenDataStructure) :: smoothedGreenData
end module smoothedGreenDataModule

!ADDED MODULE: OceanGrid Module stored interpolated ocean particles data.                      
module OceanGridDataModule
	use basicDataStructures
	implicit none
	save
	type OceanGridDataStructure
	!ALTERED value of nLat: !3!720 !!720,2161 !!1762 !!2161                                    
		integer :: nLat = 720 !200 !720
		integer :: nPhi 
		real (kind=QR_K) :: dLat, dphi, d0
		real (kind=QR_K), allocatable, dimension(:) :: sinLat,cosLat,sinPhi,cosPhi
		real (kind=QR_K), allocatable, dimension(:) :: sinLatC,cosLatC,sinPhiC,cosPhiC
		real (kind=QR_K), allocatable, dimension(:,:) :: d
		real (kind=QR_K), allocatable, dimension(:,:) :: f0
		real (kind=QR_K), allocatable, dimension(:,:) :: hs1
		real (kind=QR_K), allocatable, dimension(:,:) :: h1
	end type OceanGridDataStructure
	type(OceanGridDataStructure) :: OceanGrid

	! real(kind=QR_K), allocatable :: gradH1Lambda(:,:), gradH1Theta(:,:)
	! real(kind=QR_K), allocatable :: dotVH1(:,:), QonGrid(:,:)

	real (kind=QR_K), allocatable :: psi1onGrid(:,:)
	real (kind=QR_K), allocatable :: psi1onBand(:,:)
	real (kind=QR_K), allocatable :: velocitypsi1Band(:,:)
	real (kind=QR_K), allocatable :: psi1onNhemi(:,:)
	real (kind=QR_K), allocatable :: psi1onShemi(:,:)
	real (kind=QR_K), allocatable :: phi1onGrid(:,:)

	!TEMP:saving terms
	real (kind=QR_K), allocatable :: Term1(:,:),Term2(:,:),Term3(:,:),Term4(:,:), Term5(:,:)

	real (kind=QR_K), allocatable :: GradH1(:,:,:)
	real (kind=QR_K), allocatable :: uDotGradH1(:,:)

	real (kind=QR_K), allocatable :: VelocityMagnitude_PhiOne(:,:),VelocityMagnitude_PsiOne(:,:),VelocityMagnitude_TotOne(:,:)

	real (kind=QR_K), allocatable :: velocityOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityPsiOneOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityPhiOneOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityVorOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityDivOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: zetaLaplacianPsi(:,:)
	real (kind=QR_K), allocatable :: zetaCurlVelocity(:,:)
	real (kind=QR_K), allocatable :: zetaCurlVelocityPsi(:,:)
	
	real (kind=QR_K), allocatable :: zetaMinusF(:,:)

end module OceanGridDataModule

!ADDED Module: OceanAdaptiveGrid. ROLE: used to store values for when Qd is the source term          
module OceanAdaptiveGridDataModule
	use basicDataStructures
	implicit none
	save
	type OceanAdaptiveGridDataStructure
		!ALTERED value of nLat: 3 , 720                                               
		integer :: nLat =720 
		integer :: nPhi 
		! real(kind= QR_K) :
		real (kind=QR_K) :: dLat, dphi, d0
		real (kind=QR_K), allocatable, dimension(:) :: sinLat,cosLat,sinPhi,cosPhi, dA
		real (kind=QR_K), allocatable, dimension(:) :: sinLatC,cosLatC,sinPhiC,cosPhiC
		real (kind=QR_K), allocatable, dimension(:,:) :: d
		real (kind=QR_K), allocatable, dimension(:,:) :: hs2
		real (kind=QR_K), allocatable, dimension(:,:) :: h2
		real (kind=QR_K), allocatable :: Q(:,:)
		real (kind=QR_K), allocatable :: strengthQ(:,:)
	end type OceanAdaptiveGridDataStructure
	type(OceanAdaptiveGridDataStructure) :: OceanAdaptiveGrid

	real (kind=QR_K), allocatable :: psi2onGrid(:,:)
	real (kind=QR_K), allocatable :: psi2onBand(:,:)
	real (kind=QR_K), allocatable :: psi2onNhemi(:,:)
	real (kind=QR_K), allocatable :: psi2onShemi(:,:)
	real (kind=QR_K), allocatable :: phi2onGrid(:,:)

	real (kind=QR_K), allocatable :: VelocityMagnitude_PhiTwo(:,:),VelocityMagnitude_PsiTwo(:,:),VelocityMagnitude_TotTwo(:,:)

	real (kind=QR_K), allocatable :: velocityTwoOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityPsiTwoOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityPhiTwoOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityVorTwoOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: velocityDivTwoOnGrid(:,:,:)
	real (kind=QR_K), allocatable :: zetaLaplacianPsiTwo(:,:)
	real (kind=QR_K), allocatable :: zetaCurlVelocityTwo(:,:)
	real (kind=QR_K), allocatable :: zetaCurlVelocityPsiTwo(:,:)

end module OceanAdaptiveGridDataModule

!ADDED Module BUT not used: windDataModule. ROLE: introduce wind data from Copernicus                                                        
module windDataModule
	implicit none
	save
	type windDataStructure
		integer :: nLong=1770
		integer :: nLat=720
		INTEGER (KIND=2), allocatable :: velocity(:,:)
	end type windDataStructure
	type(windDataStructure) :: wind
end module windDataModule

!ADDED Module not used: OceanParticlesData. ROLE: storing ocean particles for when using free-grid approach                                           
! module OceanparticlesDataModule
! 	use basicDataStructures
! 	implicit none
! 	save
! 	integer, parameter :: nParticlesMAX=6300000 !!6.3 million particles max
! 	integer :: nParticles
! 	type Oceanparticle
! 		real (kind=QR_K) :: lat, phi
! 		real (kind=QR_K) :: sinLat,cosLat,sinPhi,cosPhi
! 		real (kind=QR_K) :: uLat,uLon
! 		real (kind=QR_K) :: d
! 	end type Oceanparticle
! 	type elementShell
! 		type(Oceanparticle), allocatable :: el(:)
! 	end type elementShell
! 	type(elementShell), allocatable :: Oceanparticles(:)
! end module OceanparticlesDataModule

