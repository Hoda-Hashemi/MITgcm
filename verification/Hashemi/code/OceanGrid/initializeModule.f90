module initializeModule

	implicit none
	!COMMENT: Original Legacy Code — Some obsolete write statements or commented-out lines may remain.
	public :: initialize
	public :: initializeMarshall
	public :: readBathymetry
	public :: initializeQuadIntegrationParameters

	!TOFIX: Added needs fixing
	! public :: initialize_Q_AdaptiveGrid

contains
	subroutine initialize()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		!		use greenFunctionModule
		use particlesDataModule
		use outputModule
		implicit none
		real (kind=QR_K) :: a, Lat, phi, sinLat,cosLat,cosLatPowerm
		integer :: iLat, iphi, m, i
		real (kind=QR_K) :: aOverR2,aOverLd2,factor,coriolisTerm,R2,Q
		real (kind=QR_K) :: dA, cosd,totalStrength,sinPhi0,cosPhi0,sinPhi,cosPhi
		real (kind=QR_K) :: d2, lat0, phi0, sigma0,cosLat0,sinLat0,cosDPhi,d
		real (kind=QR_K) :: weight

		lat0=pi/real(6.0,QR_K); sinLat0=sin(Lat0); cosLat0=cos(Lat0)
		phi0=pi;; sinPhi0=sin(Phi0); cosPhi0=cos(Phi0)
		sigma0=real(0.05,QR_K);


		if(allocated(particles))deallocate(particles)
		allocate(particles(nParticlesMAX))
		do i=1,nParticlesMAX
			if(allocated(particles(i)%el))deallocate(particles(i)%el)
		end do

		nParticles=0
		!Lat is latitude
		a = a_sim ! a is in m^2/s
		m = m_sim
		R2=R*R

		aOverR2=a/(R2)
		aOverLd2=a/(Ld0*Ld0)
		factor = -(real((1+m)*(2+m),QR_K)*aOverR2+aOverLd2) ! factor has the unit /sec

		!!put particles cells centers
		do iLat=1,sphericalGrid%nLat
			Lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K) !Lat here is latitude
			sinLat=sphericalGrid%sinLatC(iLat) !sin(Lat)
			cosLat=sphericalGrid%cosLatC(iLat) !cos(Lat);
			coriolisTerm = two*OMEGA*sinLat !sin(Lat) !coriolisTerm has the units of /sec
			cosLatPowerm=cosLat**m
			dA = sphericalGrid%dA(ilat)

			! write(*,*) 'dA= ', dA 
			! write(log_unit,*) dA

			do iphi=1,sphericalGrid%nPhi
				phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi-1,QR_K)
				sinPhi=sphericalGrid%sinPhiC(iPhi) !sin(Lat)
				cosPhi=sphericalGrid%cosPhiC(iPhi) !cos(Lat);

				!! Great-circle distance: which is the angle d between two points on the sphere
				cosDPhi = cosPhi*cosPhi0+sinPhi*sinPhi0 !cos(phi-phi0)
				cosDPhi = one

				cosd=sinLat*sinLat0+cosLat*cosLat0*cosDPhi

				if(cosd.gt.one)cosd=one
				if(cosd.lt.-one)cosd=-one

				d = acos(cosd)

				if(d.lt.zero) d=d+two*pi
				d=d/sigma0


				if(d.le.real(2.25,QR_K).and.sphericalGrid%d(iLat,iPhi).gt.zero)then
					d2=d*d
					weight = sphericalGrid%d(iLat,iPhi)/sphericalGrid%d0
					Q = 0 !(one*two*OMEGA*weight)*exp(-d2) !/(pi*sigma02)
					nParticles=nParticles+1
					allocate(particles(nParticles)%el(1))
					particles(nParticles)%el(1)%Lat=Lat+zero*two*sphericalGrid%dLat*sin(eight*phi)
					particles(nParticles)%el(1)%Phi=phi
					particles(nParticles)%el(1)%Q = Q !Q has the units of /sec
					particles(nParticles)%el(1)%strengthVor=Q*dA !strength has the units of m^2/sec
					particles(nParticles)%el(1)%strengthDiv=zero
					totalStrength = totalStrength + particles(nParticles)%el(1)%strengthVor
					particles(nParticles)%el(1)%sinLat=sin(Lat)
					particles(nParticles)%el(1)%cosLat=cos(Lat)
					particles(nParticles)%el(1)%sinPhi=sin(phi)
					particles(nParticles)%el(1)%cosPhi=cos(phi)
					particles(nParticles)%el(1)%uLat=zero
					particles(nParticles)%el(1)%uLon=zero
					particles(nParticles)%el(1)%hPrev = 1.0D0 !THIS IS A NON_ZERO DUMMY VALUE TO AVOID TAKING LOG OF ZERO
				end if
			end do
		end do
		return
	end subroutine initialize
	
	!COMMENT: MARSHALL'S INITIALIZATION
	subroutine initializeMarshall()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		!		use greenFunctionModule
		use particlesDataModule
		use outputModule
		implicit none
		real (kind=QR_K) :: a, Lat, phi, sinLat,cosLat,cosLatPowerm
		integer :: iLat, iphi, m, i
		real (kind=QR_K) :: aOverR2,aOverLd2,factor,coriolisTerm,psiOvera,Qminusf,R2,Q
		real (kind=QR_K) :: dA, cosd,totalStrength,sinPhi0,cosPhi0,sinPhi,cosPhi
		real (kind=QR_K) :: d0,d,d1,d2, lat0, phi0, sigma0,S2,cosLat0,sinLat0,cosDPhi,sinDPhi,denom,numer
		real (kind=QR_K) :: weight,Qmax,Qmin,colatitude,randFactor

		lat0=pi/real(6.0,QR_K); sinLat0=sin(Lat0); cosLat0=cos(Lat0)
		phi0=pi;; sinPhi0=sin(Phi0); cosPhi0=cos(Phi0)
		sigma0=real(0.05,QR_K);


		if(allocated(particles))deallocate(particles)
		allocate(particles(nParticlesMAX))
		do i=1,nParticlesMAX
			if(allocated(particles(i)%el))deallocate(particles(i)%el)
		end do

		call srand(0)

		nParticles=0
		!Lat is latitude
		a = a_sim ! a is in m^2/s
		m = m_sim
		R2=R*R
		aOverR2=a/(R2)
		aOverLd2=a/(Ld0*Ld0)
		factor = -(real((1+m)*(2+m),QR_K)*aOverR2+aOverLd2) ! factor has the unit /sec
		Qmax=-1000000000.0
		Qmin=1000000000.0

		!put particles cells centers
		do iLat=1,sphericalGrid%nLat
			Lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K) !Lat here is latitude
			sinLat=sphericalGrid%sinLatC(iLat) !sin(Lat)
			cosLat=sphericalGrid%cosLatC(iLat) !cos(Lat);
			coriolisTerm = two*OMEGA*sinLat !sin(Lat) !coriolisTerm has the units of /sec
			cosLatPowerm=cosLat**m
			dA = sphericalGrid%dA(ilat)
			colatitude=half*pi-Lat;
			do iphi=1,sphericalGrid%nPhi
				phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi-1,QR_K)
				psiOvera=-sinLat*cosLatPowerm*cos(real(m,QR_K)*phi) !psi has units of m^2/s, psiOvera is unitless

				randFactor = one+(one-two*rand())/ten
				Qminusf=randFactor*factor*psiOvera !Qminusf has the units of /sec
				Q=Qminusf+coriolisTerm  !Q has the units of /sec

			!if (iLat.eq.sphericalGrid%nLat/2)write(*,*)iPhi,Q

				nParticles=nParticles+1
				allocate(particles(nParticles)%el(1))
				particles(nParticles)%el(1)%Lat=Lat
				particles(nParticles)%el(1)%Phi=phi
				particles(nParticles)%el(1)%Q=Q !Q has the units of /sec

				!particles(nParticles)%el(1)%Q=(factor*R*R/real(24,QR_K))*((sin(-colatitude&
				!&+half*sphericalGrid%dLat)**6)*(sin(-four*phi+two*sphericalGrid%dphi)&
				!&+sin(four*phi+two*sphericalGrid%dphi))&
				!&-(sin(colatitude+half*sphericalGrid%dLat)**6)*(sin(-four*phi+two*sphericalGrid%dphi)&
				!&+sin(four*phi+two*sphericalGrid%dphi)))

				if(abs(Q).gt.Qmax)Qmax=abs(Q)
				if(abs(Q).lt.Qmin)Qmin=abs(Q)

				particles(nParticles)%el(1)%strengthVor=Q*dA !strength has the units of m^2/sec


				particles(nParticles)%el(1)%strengthDiv=zero
				totalStrength = totalStrength + particles(nParticles)%el(1)%strengthVor
				particles(nParticles)%el(1)%sinLat=sin(Lat)
				particles(nParticles)%el(1)%cosLat=cos(Lat)
				particles(nParticles)%el(1)%sinPhi=sin(phi)
				particles(nParticles)%el(1)%cosPhi=cos(phi)
				particles(nParticles)%el(1)%uLat=zero
				particles(nParticles)%el(1)%uLon=zero
				particles(nParticles)%el(1)%hPrev = 1.0D0 !THIS IS A NON_ZERO DUMMY VALUE TO AVOID TAKING LOG OF ZERO
			end do
		end do

		write(*,*)'number of particles = ',nParticles
		write(*,*)'totalStrength=',totalStrength
		write(*,*)'Qmax,Qmin=',Qmax,Qmin,two*OMEGA

		return
	end subroutine initializeMarshall

	!COMMENT: Reading Bathymetry 
	subroutine readBathymetry()
		use bathymetryDataModule
		implicit none
		integer :: iLon,iLat
		logical :: file_exists

		if(.NOT.allocated(bathymetry%z))allocate(bathymetry%z(bathymetry%nLong,bathymetry%nLat))
		write(*,*)'Reading Bathymetry from ../resources/topography/etopo5.csv'

		INQUIRE(FILE='../resources/topography/etopo5.csv', EXIST=file_exists)
		if(file_exists)then
			open(UNIT=2,FILE='../resources/topography/etopo5.csv')
			read(2,*)
				do iLat=1,bathymetry%nLat
					do iLon=1,bathymetry%nLong
					read(2,*)bathymetry%z(iLon,iLat)
					!				write(*,*)topo
					!			write(*,*)greenFunctionData%greenTableNEntries
				end do
			end do
			close(2)
		else
			write(*,*)'bathymetry file does not exist'
			stop
		end if
		write(*,*)'*** Done reading Bathymetry ***'
		return
	end subroutine readBathymetry

	subroutine initializeQuadIntegrationParameters()
		use quadIntegrationDataModule
		implicit none

    	!Niq=9
		weight9(1)=real(16,QR_K)/real(81,QR_K); dx9(1)=zero; dy9(1)=zero;
		weight9(2)=real(25,QR_K)/real(324,QR_K); dx9(2)=-sqrt(three/five); dy9(2)=-sqrt(three/five);
		weight9(3)=real(25,QR_K)/real(324,QR_K); dx9(3)=sqrt(three/five); dy9(3)=-sqrt(three/five);
		weight9(4)=real(25,QR_K)/real(324,QR_K); dx9(4)=sqrt(three/five); dy9(4)=sqrt(three/five);
		weight9(5)=real(25,QR_K)/real(324,QR_K); dx9(5)=-sqrt(three/five); dy9(5)=sqrt(three/five);
		weight9(6)=real(10,QR_K)/real(81,QR_K); dx9(6)=zero; dy9(6)=sqrt(three/five);
		weight9(7)=real(10,QR_K)/real(81,QR_K); dx9(7)=zero; dy9(7)=-sqrt(three/five);
		weight9(8)=real(10,QR_K)/real(81,QR_K); dx9(8)=-sqrt(three/five); dy9(8)=zero
		weight9(9)=real(10,QR_K)/real(81,QR_K); dx9(9)=sqrt(three/five); dy9(9)=zero
		dx9=half*dx9; dy9=half*dy9;

		!Niq=4
		weight4(1)=real(1,QR_K)/real(4,QR_K); dx4(1)=-sqrt(one/three); dy4(1)=-sqrt(one/three);
		weight4(2)=real(1,QR_K)/real(4,QR_K); dx4(2)=-sqrt(one/three); dy4(2)=sqrt(one/three);
		weight4(3)=real(1,QR_K)/real(4,QR_K); dx4(3)=sqrt(one/three); dy4(3)=-sqrt(one/three);
		weight4(4)=real(1,QR_K)/real(4,QR_K); dx4(4)=sqrt(one/three); dy4(4)=sqrt(one/three);
		dx4=half*dx4; dy4=half*dy4;

		!Niq=1;
		weight1(1)=one
		return
	end subroutine initializeQuadIntegrationParameters
	
	!TOFIX: (I do not advise using this subroutine yet still in progress and i was experimenting so...)
	! subroutine initialize_Q_AdaptiveGrid()
	! 	use OceanAdaptiveGridDataModule
	! 	use simulationParameters
	! 	implicit none

	! 	integer :: iLat, iPhi, jLat, jPhi, jPhiWrap
	! 	integer :: nBlock, iLatB, iLatT, iPhiL, iPhiR,nBuffersOcean
	! 	logical :: bufferIsAllOcean, bufferFilled
	! 	real (kind=QR_K) :: Queow

	! 	if(allocated(OceanAdaptiveGrid%Q))deallocate(OceanAdaptiveGrid%Q)
	! 	if(.NOT.allocated(OceanAdaptiveGrid%Q))allocate(OceanAdaptiveGrid%Q(OceanAdaptiveGrid%nLat,OceanAdaptiveGrid%nPhi))
	! 	OceanAdaptiveGrid%Q=zero

	! 	if(allocated(OceanAdaptiveGrid%strengthQ))deallocate(OceanAdaptiveGrid%strengthQ)
	! 	if(.NOT.allocated(OceanAdaptiveGrid%strengthQ))allocate(OceanAdaptiveGrid%strengthQ(OceanAdaptiveGrid%nLat,OceanAdaptiveGrid%nPhi))
	! 	OceanAdaptiveGrid%strengthQ=zero

	! 	Queow = twothousand

	! 	nBlock = 19 !!ODD
	! 	write(*,*) 'Checking buffer of size =', nBlock, 'x ', nBlock
		
	! 	bufferFilled = .false. 
		
	! 	nBuffersOcean=0
		
	! 	!* Loop over grid
	! 	do iLat = 1, OceanAdaptiveGrid%nLat
	! 		if (bufferFilled) exit
	! 		do iPhi = 1, OceanAdaptiveGrid%nPhi
	! 			if (bufferFilled) exit
	! 			!* enter the loop if the core is ocean
	! 			if (OceanAdaptiveGrid%d(iLat,iPhi).gt.zero)then
	! 				!*Define buffer limits
	! 				iLatB = max(iLat - nBlock/2, 1)
	! 				iLatT = min(iLat + nBlock/2, OceanAdaptiveGrid%nLat)
	! 				iPhiL = iPhi - nBlock/2
	! 				iPhiR = iPhi + nBlock/2
		
	! 				bufferIsAllOcean = .true.
		
	! 				!* Loop over the buffer
	! 				do jLat = iLatB, iLatT
	! 					do jPhi = iPhiL, iPhiR
	! 						jPhiWrap = jPhi
	! 						if (jPhiWrap < 1) jPhiWrap = jPhiWrap + OceanAdaptiveGrid%nPhi
	! 						if (jPhiWrap > OceanAdaptiveGrid%nPhi) jPhiWrap = jPhiWrap - OceanAdaptiveGrid%nPhi

	! 						if (OceanAdaptiveGrid%d(jLat, jPhiWrap)<= 0) then
	! 							bufferIsAllOcean = .false.
	! 							exit
	! 						end if
	! 					end do
	! 					if (.not. bufferIsAllOcean) exit
	! 				end do

	! 				if (bufferIsAllOcean) then
	! 					nBuffersOcean=nBuffersOcean+1
	! 					write(*,*) 'Found a buffer at ilat = ', iLat, 'and iLon = ', iPhi ,'and boundaries of (iLatB,iLatT,iPhiL,iPhiR) = ', iLatB,iLatT,iPhiL,iPhiR

	! 					! Fill Q in the buffer
	! 					do jLat = iLatB, iLatT
	! 						do jPhi = iPhiL, iPhiR
	! 							jPhiWrap = jPhi
	! 							if (jPhiWrap < 1) jPhiWrap = jPhiWrap + OceanAdaptiveGrid%nPhi
	! 							if (jPhiWrap > OceanAdaptiveGrid%nPhi) jPhiWrap = jPhiWrap - OceanAdaptiveGrid%nPhi

	! 							OceanAdaptiveGrid%Q(jLat, jPhiWrap) = Queow
	! 							OceanAdaptiveGrid%strengthQ(jLat, jPhiWrap) = OceanAdaptiveGrid%Q(jLat, jPhiWrap)*OceanAdaptiveGrid%dA(jLat)
	! 							write(*,*) 'Q ===== ', OceanAdaptiveGrid%Q(jLat, jPhiWrap)
	! 						end do
	! 					end do
	! 					bufferFilled = .true.!  ! mark buffer as filled
	! 					exit 

	! 				end if

	! 			end if

	! 		end do
	! 	end do

	! 	write(*,*) 'Buffer check complete with nBuffersOcean = ', nBuffersOcean
	! end subroutine initialize_Q_AdaptiveGrid

end module initializeModule





