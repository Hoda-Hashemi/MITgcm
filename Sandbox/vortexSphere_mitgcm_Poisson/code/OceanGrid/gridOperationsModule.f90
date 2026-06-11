module gridOperationsModule
	!COMMENT: the subroutines are oreded the same order here.
	implicit none
	!COMMENT: Original Legacy Code — Some obsolete write statements or commented-out lines may remain.
    public :: constructGrid
	public :: ReorderNew
	public :: reorderOnGrid
	public :: interpolateFromParticlesToGridSerial
	public :: interpolate_Q_D_FromParticlesToGrid
	public :: interpolate_Q_hs_FromParticlesToGrid
	public :: interpolate_hs_FromParticlesToGrid
	public :: interpolate_Q_FromParticlesToGrid
	public :: interpolate_D_FromParticlesToGrid
	public :: compute_d_h_h_s_div_AtParticlesLocations
	public :: compute_d_AtParticlesLocations
	!	public :: computeh0hAtParticlesLocations
	public :: initialize_h_d_hPrev_AtParticlesLocations
	public :: interpolate_V_FromGridToParticles
	public :: regrid
	public :: markActiveGridCells
	public :: interpolateBathymetryToGrid
	public :: setQonGridToZero

	!COMMENT: The following subroutines are added and not all of them are finalized. 
	public :: interpolateBathymetryToOceanGrid !!FINALIZED
	public :: constructOceanGrid !!FINALIZED
	
	public :: constructAdaptiveGrid
	public :: interpolateBathymetryToAdaptiveGrid
	! public :: AdvectionQAdaptiveGrid

contains

	subroutine constructGrid()
		use sphericalGridDataModule 
		use simulationParameters
		implicit none
		integer :: iphi,ilat
		real (kind=QR_K) :: lat, phi, colatitude

		!		sphericalGrid%dLat=min(half*Ld0/R,pi/real(sphericalGrid%nLat,QR_K))
		!		sphericalGrid%nLat = int(pi/sphericalGrid%dLat+real(1.0D-6,QR_K))

		!		sphericalGrid%nLat = 360
		!		sphericalGrid%nPhi=2*sphericalGrid%nLat
		sphericalGrid%nPhi=2*sphericalGrid%nLat
		sphericalGrid%dLat=pi/real(sphericalGrid%nLat,QR_K)
		sphericalGrid%dphi=two*pi/real(sphericalGrid%nPhi,QR_K)

		write(*,*)'ntheta , nPhi=',sphericalGrid%nLat,sphericalGrid%nPhi

		if(.NOT.allocated(sphericalGrid%sinLat))allocate(sphericalGrid%sinLat(sphericalGrid%nLat+1))
		if(.NOT.allocated(sphericalGrid%cosLat))allocate(sphericalGrid%cosLat(sphericalGrid%nLat+1))
		if(.NOT.allocated(sphericalGrid%sinLatC))allocate(sphericalGrid%sinLatC(sphericalGrid%nLat))
		if(.NOT.allocated(sphericalGrid%cosLatC))allocate(sphericalGrid%cosLatC(sphericalGrid%nLat))
		if(.NOT.allocated(sphericalGrid%dA))allocate(sphericalGrid%dA(sphericalGrid%nLat))
		if(.NOT.allocated(sphericalGrid%sinPhi))allocate(sphericalGrid%sinPhi(sphericalGrid%nPhi))
		if(.NOT.allocated(sphericalGrid%cosPhi))allocate(sphericalGrid%cosPhi(sphericalGrid%nPhi))
		if(.NOT.allocated(sphericalGrid%sinPhiC))allocate(sphericalGrid%sinPhiC(sphericalGrid%nPhi))
		if(.NOT.allocated(sphericalGrid%cosPhiC))allocate(sphericalGrid%cosPhiC(sphericalGrid%nPhi))

		!!cell vertices
		do ilat=1,sphericalGrid%nLat+1
			lat=-half*pi+sphericalGrid%dLat*real(ilat-1,QR_K)
			sphericalGrid%sinLat(ilat)=sin(lat)
			sphericalGrid%cosLat(ilat)=cos(lat)
		end do

		do iphi=1,sphericalGrid%nPhi
			phi=sphericalGrid%dphi*real(iphi-1,QR_K)
			sphericalGrid%sinPhi(iphi)=sin(phi)
			sphericalGrid%cosPhi(iphi)=cos(phi)
		end do

		!!cell centers
		do ilat=1,sphericalGrid%nLat
			lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(ilat-1,QR_K)
			sphericalGrid%sinLatC(ilat)=sin(lat)
			sphericalGrid%cosLatC(ilat)=cos(lat)
		end do

		do iphi=1,sphericalGrid%nPhi
			phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi-1,QR_K)
			sphericalGrid%sinPhiC(iphi)=sin(phi)
			sphericalGrid%cosPhiC(iphi)=cos(phi)
		end do

		do ilat=1,sphericalGrid%nLat
			lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(ilat-1,QR_K)
			colatitude=half*pi-Lat;
			sphericalGrid%dA(ilat)=(cos(colatitude-half*sphericalGrid%dLat)-cos(colatitude+half*sphericalGrid%dLat))*sphericalGrid%dphi !*R*R
		end do

		if(.not.include_bathymetry)then
				if(allocated(sphericalGrid%d))deallocate(sphericalGrid%d)
				if(.NOT.allocated(sphericalGrid%d))allocate(sphericalGrid%d(sphericalGrid%nLat,sphericalGrid%nPhi))
				!ALTERED: filling the matrix with constant value
				sphericalGrid%d0 = ((two*OMEGA*Ld0)**2)/gravity
				do ilat=1,sphericalGrid%nLat
					do iphi=1,sphericalGrid%nPhi
						sphericalGrid%d(ilat,iPhi) = sphericalGrid%d0
					end do
				end do

				write(*,*)'sphericalGrid%d0=',sphericalGrid%d0
			end if
		!stop

		return
	end subroutine constructGrid

	subroutine ReorderNew()
		use QuickSearchVariables
		use particlesDataModule
		use sphericalGridDataModule
		use simulationParameters
		implicit none

		integer :: i, j
		integer :: m, l, k
		integer :: iPhi, iLat, izone
		real (kind=QR_K) :: phiC, latC, totalStrength
		integer, allocatable, dimension(:) :: counter

		zoneSize = gammaCutoffFactor*Ld0/R

		nLatZones = INT(pi/zoneSize)
		write(*,*)pi/DBLE(nLatZones),zoneSize
		zoneSize=pi/real(nLatZones,QR_K)

		nPhiZones = 2*nLatZones
		write(*,*)two*pi/real(nPhiZones,QR_K),zoneSize

		phiMinQS = 0
		phiMaxQS = two*pi
		latMinQS = -half*pi
		latMaxQS = half*pi

		nZones=nPhiZones*nLatZones

		if(allocated(zones))deallocate(zones)
		allocate(zones(nZones))
		allocate(counter(nZones))

		! initialize number of elements of zone i to zero
		zones(1:nZones)%nEls=0
		counter(1:nZones)=0

		! The zones are indexed from 1 to nZones marching first in the lat direction then long
		! march over the elements to determine:
		! (a) number of elements in each zone
		! (b) the zone to which element i belongs to
		do i=1,nParticles
			iPhi=int((particles(i)%el(1)%Phi-phiMinQS)/zoneSize)+1
			if(iPhi.lt.1.or.iPhi.gt.nPhiZones)then
				write(*,*)'stopping in reorder'
				write(*,*)'phi ',iPhi,nPhiZones,particles(i)%el(1)%Phi
				stop
			end if
			iLat=int((particles(i)%el(1)%lat-latMinQS)/zoneSize)+1
			if(iLat.lt.1.or.iLat.gt.nLatZones)then
				write(*,*)'stopping in reorder'
				write(*,*)'lat ',iLat,nLatZones,particles(i)%el(1)%lat
				stop
			end if
			!izone=(ix-1)*nBlck + (iz-1)*nyZones + iy
			izone=(iPhi-1)*nLatZones + iLat
			zones(izone)%nEls=zones(izone)%nEls+1
			particles(i)%el(1)%indexZone=izone
		end do

		! store the index of the fist element in each zone in startLoc
		zones(1)%startLoc=1
		do i=2,nZones
			zones(i)%startLoc=zones(i-1)%startLoc+zones(i-1)%nEls
		end do

		! reord all elements in zones
		do i=1,nParticles
			izone=particles(i)%el(1)%indexZone
			counter(izone)=counter(izone)+1
			j=zones(izone)%startLoc+counter(izone)-1
			particles(j)%el(1)%index = i
		end do

		deallocate(counter)

		k=0
		do m=-1,1
			do l=-1,1
				k=k+1
				neighZoneRelInd(k) = l + m*nLatZones
			end do
		end do

		return
	end subroutine ReorderNew

	subroutine reorderOnGrid()
		use QuickSearchVariables
		use particlesDataModule
		use sphericalGridDataModule
		use simulationParameters
		implicit none

		integer :: i, j
		integer :: m, l, k
		integer :: iPhi, iLat, izone
		real (kind=QR_K) :: phiC, latC, totalStrength
		integer, allocatable, dimension(:) :: counter

		zoneSize = sphericalGrid%dLat
		nLatZones = sphericalGrid%nLat
		nPhiZones = sphericalGrid%nPhi
		phiMinQS = zero
		latMinQS = -half*pi
		nZones=nPhiZones*nLatZones

		if(allocated(zones))deallocate(zones)
		allocate(zones(nZones))
		if(allocated(counter))deallocate(counter)
		allocate(counter(nZones))

		! initialize number of elements of zone i to zero
		zones(1:nZones)%nEls=0
		counter(1:nZones)=0

		! The zones are indexed from 1 to nZones marching first in the lat direction then long
		! march over the elements to determine:
		! (a) number of elements in each zone
		! (b) the zone to which element i belongs to
		do i=1,nParticles
			!! Addition: Reflect latitude if it exceeds bounds
			! if (particles(i)%el(1)%lat > half*pi) then
			! 	particles(i)%el(1)%lat = pi - particles(i)%el(1)%lat
			! 	particles(i)%el(1)%Phi = pi + particles(i)%el(1)%Phi
				! particles(i)%el(1)%uLat = -particles(i)%el(1)%uLat
				! particles(i)%el(1)%sinLat = sin(particles(i)%el(1)%lat)
				! particles(i)%el(1)%cosLat = cos(particles(i)%el(1)%lat)
				
			! else if (particles(i)%el(1)%lat < -half*pi) then
			! 	particles(i)%el(1)%lat = -pi - particles(i)%el(1)%lat
			! 	particles(i)%el(1)%Phi = pi + particles(i)%el(1)%Phi
				! particles(i)%el(1)%uLat = -particles(i)%el(1)%uLat
				! particles(i)%el(1)%sinLat = sin(particles(i)%el(1)%lat)
				! particles(i)%el(1)%cosLat = cos(particles(i)%el(1)%lat)
			! end if

			iPhi=int((particles(i)%el(1)%Phi-phiMinQS)/zoneSize)+1
			iLat=int((particles(i)%el(1)%lat-latMinQS)/zoneSize)+1

			!!uncommented, debugging reasons
			if(iPhi.lt.1.or.iPhi.gt.nPhiZones .or. iLat.lt.1.or.iLat.gt.nLatZones)then
				write(*,*)'stopping in reorder'
				write(*,*) 'out of bounds: i=',i , 'iPhi=',iPhi,'iLat=',iLat , &
				'Phi=',particles(i)%el(1)%Phi,'lat=',particles(i)%el(1)%lat
				
			end if
			! if(iPhi.lt.1.or.iPhi.gt.nPhiZones)then
			! 	write(*,*)'stopping in reorder'
			! 	write(*,*)'phi ',iPhi,nPhiZones,particles(i)%el(1)%Phi
			! 	stop
			! end if
			! if(iLat.lt.1.or.iLat.gt.nLatZones)then
			! 	write(*,*)'stopping in reorder'
			! 	write(*,*)'lat ',iLat,nLatZones,particles(i)%el(1)%lat
			! 	stop
			! end if
			izone=(iPhi-1)*nLatZones + iLat
			zones(izone)%nEls=zones(izone)%nEls+1
			particles(i)%el(1)%indexZone=izone
		end do

		! store the index of the fist element in each zone in startLoc
		zones(1)%startLoc=1
		do i=2,nZones
			zones(i)%startLoc=zones(i-1)%startLoc+zones(i-1)%nEls
		end do

		! reorder all elements in zones
		do i=1,nParticles
			izone=particles(i)%el(1)%indexZone
			counter(izone)=counter(izone)+1
			j=zones(izone)%startLoc+counter(izone)-1
			particles(j)%el(1)%index = i
		end do

		deallocate(counter)

		k=0
		do m=-1,1
			do l=-1,1
				k=k+1
				neighZoneRelInd(k) = l + m*nLatZones
			end do
		end do

		return
	end subroutine reorderOnGrid

	subroutine interpolateFromParticlesToGridSerial()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		!$ use omp_lib
		implicit none
		integer :: iLat, iphi, i, cutoff, sigma,iLat0,iphi0,iPhiG
		real (kind=QR_K) :: d,d2, Lat0, phi0, sigma2, coef, phi0G, Lat, phi,sinLat,cosd
		real (kind=QR_K) :: Lat0G,cosLat,cosLat0,sinLat0,cosDPhi,sinDPhi,denom,numer,totalStrength
		real (kind=QR_K) :: cosPhi,cosPhi0,sinPhi,sinPhi0,dum,cosDPhi1
		integer :: CHUNK, i1, i2, threadID, nThreads
		integer :: kooki

		if(.NOT.allocated(QOnGrid))allocate(QOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		QOnGrid=zero

		cutoff = 4 !stop at cutoff*sigma
		sigma = 2! sigma = 2.0*spacing, spacing = dtheta=dphi
		sigma2=(2.0*sphericalGrid%dLat)**2;
		coef = 1.0/(pi*sigma2)
		totalStrength = 0.0
		kooki=0

		i1=1; i2=nParticles

		!		nThreads = OMP_get_max_threads()
		!		call OMP_SET_NUM_THREADS(nThreads)
		!		CHUNK = nParticles/nThreads
		!		write(*,*)CHUNK
		!		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!		!$OMP &coef,nParticles,cutoff,sigma,sigma2,QOnGrid)
		!		threadID = omp_get_thread_num()
		!		i1 = threadID*CHUNK+1
		!		i2=i1-1+CHUNK
		!		if(threadID+1.eq.nThreads)i2=nParticles

		do i=i1,i2 !1,nParticles
			Lat0 = particles(i)%el(1)%Lat
			cosLat0=particles(i)%el(1)%cosLat
			sinLat0=particles(i)%el(1)%sinLat
			iLat0 = 1+int((Lat0+half*pi)/sphericalGrid%dLat)
			Lat0G=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat0-1,QR_K)

			phi0 = particles(i)%el(1)%Phi
			cosPhi0=particles(i)%el(1)%cosPhi
			sinPhi0=particles(i)%el(1)%sinPhi
			iphi0 = 1+int(phi0/sphericalGrid%dphi)
			phi0G= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi0-1,QR_K)

			!			if((2.0*pi-phi0).lt.0.75*sphericalGrid%dLat.and.abs(lat0-pi/6.0).lt.5.0*sphericalGrid%dLat)&
			!			&write(*,*)i,Lat0,Lat0G,Lat0-Lat0G,phi0,phi0G,phi0-phi0G,iphi0,ilat0
			!			if(phi0.lt.0.75*sphericalGrid%dLat.and.abs(lat0-pi/6.0).lt.5.0*sphericalGrid%dLat)&
			!			&write(*,*)i,Lat0,Lat0G,Lat0-Lat0G,phi0,phi0G,phi0-phi0G,iphi0,ilat0

			do iLat=max(iLat0-cutoff*sigma,1),min(iLat0+cutoff*sigma,sphericalGrid%nLat)
				Lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K) !Lat here is latitude
				cosLat=sphericalGrid%cosLatC(iLat)
				sinLat=sphericalGrid%sinLatC(iLat)
				dum = 0.0
				!S2=(sin(half*(Lat-lat0)))**2
				do iphi = -8,8
					phi=phi0G+iphi*sphericalGrid%dphi
					!iphiG = 1+int(phi/sphericalGrid%dphi)
					iPhiG = iphi0+iphi
					if(iPhiG.lt.1)iPhiG=iphiG+sphericalGrid%nPhi
					if(iPhiG.gt.sphericalGrid%nPhi)iPhiG=iphiG-sphericalGrid%nPhi

					cosPhi=sphericalGrid%cosPhiC(iPhiG)
					sinPhi=sphericalGrid%sinPhiC(iPhiG)
					cosDPhi = cosPhi*cosPhi0+sinPhi*sinPhi0 !

					!cosDPhi1 = cos(phi-phi0)
					!if(abs(cosDPhi1-cosDPhi).gt.1.0E-10)then
					!	write(*,*)sinPhi-sin(phi),cosPhi-cos(phi),sinPhi0-sin(phi0),cosPhi0-cos(phi0)
					!	stop
					!end if

					cosd=sinLat*sinLat0+cosLat*cosLat0*cosDPhi
					!					write(*,*)iLat-iLat0,iPhi,sinLat-sinLat0,cosLat-cosLat0,cosDPhi,cosd
					!					write(*,*)
					if(cosd.gt.1.0)cosd=1.0
					if(cosd.lt.-1.0)cosd=-1.0
					d = acos(cosd)

					!if(ISNAN(d))then
					!write(*,*)'XXXXXXXX',sinLat,sinLat0,cosLat,cosLat0,cosDPhi,sinLat*sinLat0+cosLat*cosLat0*cosDPhi
					!end if

					if(d.lt.0.0)d=d+two*pi
					d2=d*d/sigma2
					QOnGrid(iLat,iPhiG) = QOnGrid(iLat,iPhiG) + coef*particles(i)%el(1)%strengthVor*exp(-d2)
					dum = dum + coef*particles(i)%el(1)%strengthVor*exp(-d2)
					!if(ISNAN(QOnGrid(iLat,iPhiG)).or.ISNAN(dum))then
					!write(*,*)'&&&&&&',coef,particles(i)%el(1)%strengthVor,d2,exp(-d2)
					!write(*,*)d,cosd,sinLat,sinLat0,cosLat,cosLat0,cosDPhi
					!stop
					!end if
				end do
				totalStrength = totalStrength + dum*sphericalGrid%dA(iLat)
			end do
		end do
		!		!$OMP END PARALLEL

		!		do iLat=1,sphericalGrid%nLat
		!			dum=0.0
		!			do iPhi=1,sphericalGrid%nPhi
		!!				write(*,*)QOnGrid(iLat,iPhi)
		!				dum = dum +  QOnGrid(iLat,iPhi)
		!			end do
		!			totalStrength = totalStrength + dum*sphericalGrid%dA(iLat)
		!		end do
		! write(*,*)'totalStrength after interpolating from particles to grid =',totalStrength

		return
	end subroutine interpolateFromParticlesToGridSerial

	subroutine markActiveGridCells()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		!$ use omp_lib
		implicit none
		integer :: iLat, iphi, i, cutoff, sigma,iLat0,iphi0
		real (kind=QR_K) :: lat0, phi0
		integer :: CHUNK, i1, i2, threadID, nThreads
		integer :: ilatMax,ilatMin,iPhiMax,iPhiMin,iPhi1

		if(.NOT.allocated(isActiveCell))allocate(isActiveCell(sphericalGrid%nLat,sphericalGrid%nPhi))
		isActiveCell = .FALSE.

		cutoff = 4 !stop at cutoff*sigma
		sigma = 2! sigma = 2.0*spacing, spacing = dtheta=dphi
		nActiveCells = 0

		i1=1; i2=nParticles

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nParticles/nThreads
		write(*,*)CHUNK
		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &nParticles,cutoff,sigma,isActiveCell)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nParticles

		do i=i1,i2 !1,nParticles
			Lat0 = particles(i)%el(1)%Lat
			iLat0 = 1+int((Lat0+half*pi)/sphericalGrid%dLat)
			phi0 = particles(i)%el(1)%Phi
			iphi0 = 1+int(phi0/sphericalGrid%dphi)

			iLatMin=max(iLat0-cutoff*sigma,1)
			iLatMax=min(iLat0+cutoff*sigma,sphericalGrid%nLat)
			iPhiMin= iphi0-8;
			iPhiMax= iphi0+8;

			!write(*,*)threadID,i,i1,i2
			do iLat=iLatMin,iLatMax
				if(iLat.lt.1.or.iLat.gt.sphericalGrid%nLat)stop "ilat"
				do iPhi=iPhiMin,iPhiMax
					iPhi1=iPhi
					if(iPhi.lt.1)iPhi1=iPhi+sphericalGrid%nPhi
					if(iPhi.gt.sphericalGrid%nPhi)iPhi1=iPhi-sphericalGrid%nPhi
					!isActiveCell(iLat,iPhi1)=.TRUE.  !removed 28FEB19
					if(sphericalGrid%d(iLat,iPhi1).gt.zero)isActiveCell(iLat,iPhi1)=.TRUE. !added 28FEB19
				end do
			end do
		end do
		!$OMP END PARALLEL

		nActiveCells=0
		do ilat=1,sphericalGrid%nLat
			do iphi=1,sphericalGrid%nPhi
				if(isActiveCell(iLat,iPhi))nActiveCells=nActiveCells+1
			end do
		end do

		if(allocated(activeCells))deallocate(activeCells)
		allocate(activeCells(nActiveCells,2))
		nActiveCells=0
		do ilat=1,sphericalGrid%nLat
			do iphi=1,sphericalGrid%nPhi
				if(isActiveCell(iLat,iPhi))then
					nActiveCells=nActiveCells+1
					activeCells(nActiveCells,1)=iLat
					activeCells(nActiveCells,2)=iPhi
				end if
			end do
		end do

		if(allocated(isActiveCell))deallocate(isActiveCell)

		write(*,*)'number of active cells =',nActiveCells
		return
	end subroutine markActiveGridCells

	subroutine interpolate_Q_D_FromParticlesToGrid()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		use QuickSearchVariables
		!$ use omp_lib
		implicit none
		integer :: iLat, iphi, i, cutoff, sigma !,iLat0,iphi0,iPhiG
		real (kind=QR_K) :: cosLat,cosLat0,sinLat,sinLat0,lat0,phi0
		real (kind=QR_K) :: cosPhi,cosPhi0,sinPhi,sinPhi0,dum,cosDPhi
		real (kind=QR_K) :: sigma2,d2,d,cosd,coef
		real (kind=QR_K) :: totalStrength, dumDiv, term, totalDivStrength
		integer :: CHUNK, i1, i2, threadID, nThreads
		integer :: kooki
		integer :: iL, iP, iP1, iS, j1, j2, jj, izone

		write(*,*)'       markActiveGridCells()';call markActiveGridCells()
		write(*,*)'       reorderOnGrid()';call reorderOnGrid()

		if(.NOT.allocated(QOnGrid))allocate(QOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		QOnGrid=zero
		if(.NOT.allocated(divOnGrid))allocate(divOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		divOnGrid=zero

		cutoff = 4 !stop at cutoff*sigma
		sigma = 2! sigma = 2.0*spacing, spacing = dtheta=dphi
		sigma2=(2.0*sphericalGrid%dLat)**2;
		coef = 1.0/(pi*sigma2)
		totalStrength = 0.0
		kooki=0

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nActiveCells/nThreads
		!		write(*,*)'CHUNK',CHUNK
		i1=1; i2=nActiveCells

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &coef,nParticles,cutoff,sigma,sigma2,QOnGrid,nActiveCells,activeCells,&
		!$OMP &zones, nLatZones,nZones,divOnGrid)

		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nActiveCells

		do i=i1,i2
			!			write(*,*)threadID,i,i1,i2

			iLat = activeCells(i,1)
			iPhi = activeCells(i,2)

			cosLat=sphericalGrid%cosLatC(iLat)
			sinLat=sphericalGrid%sinLatC(iLat)
			cosPhi=sphericalGrid%cosPhiC(iPhi)
			sinPhi=sphericalGrid%sinPhiC(iPhi)
			dum = zero
			dumDiv = zero
			!if(1.eq.2)then
			do iL=max(iLat-sigma*cutoff,1),min(iLat+sigma*cutoff,sphericalGrid%nLat)
				do iP1=iPhi-sigma*cutoff,iPhi+sigma*cutoff
					iP=iP1
					if(iP1.lt.1)iP=iP1+sphericalGrid%nPhi
					if(iP1.gt.sphericalGrid%nPhi)iP=iP1-sphericalGrid%nPhi
					izone=(iP-1)*nLatZones + iL
					!					if(izone.lt.1.or.izone.gt.nZones)stop "dfdf"
					if(zones(izone)%nEls.gt.0) then
						j1 = zones(izone)%startLoc
						j2 = j1+zones(izone)%nEls-1
						do jj = j1, j2 !sources
							iS = particles(jj)%el(1)%index
							Lat0 = particles(iS)%el(1)%Lat
							cosLat0=particles(iS)%el(1)%cosLat
							sinLat0=particles(iS)%el(1)%sinLat
							phi0 = particles(iS)%el(1)%Phi
							cosPhi0=particles(iS)%el(1)%cosPhi
							sinPhi0=particles(iS)%el(1)%sinPhi
							cosDPhi = cosPhi*cosPhi0+sinPhi*sinPhi0 !
							cosd=sinLat*sinLat0+cosLat*cosLat0*cosDPhi
							if(cosd.gt.1.0)then
								cosd=1.0
							else if(cosd.lt.-1.0) then
								cosd=-1.0
							end if
							d = acos(cosd)
							if(d.lt.0.0)d=d+two*pi
							d2=d*d/sigma2
							term = coef*exp(-d2)
							dum = dum + term*particles(iS)%el(1)%strengthVor
							dumDiv = dumDiv + term*particles(iS)%el(1)%strengthDiv
							! INTERPOLATE hs
							!							hs = particles(iS)%el(1)%h-particles(iS)%el(1)%d
						end do
					end if
				end do
			end do

			QOnGrid(iLat,iPhi) = dum
			divOnGrid(iLat,iPhi) = dumDiv
			!end if
		end do
		!$OMP END PARALLEL

		!write(*,*)'here'

		!stop

		totalStrength = zero
		totalDivStrength = zero
		do i=1,nActiveCells
			iLat = activeCells(i,1)
			iPhi = activeCells(i,2)
			totalStrength = totalStrength+QOnGrid(iLat,iPhi)*sphericalGrid%dA(iLat)
			totalDivStrength = totalDivStrength+divOnGrid(iLat,iPhi)*sphericalGrid%dA(iLat)
		end do
		write(*,*)'totalVorStrength and totalDivStrength after interpolating from particles to grid =',totalStrength, totalDivStrength

		if(allocated(activeCells))deallocate(activeCells)

		return
	end subroutine interpolate_Q_D_FromParticlesToGrid
	
	subroutine interpolate_Q_hs_FromParticlesToGrid()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		use QuickSearchVariables
		!$ use omp_lib
		implicit none
		integer :: iLat, iphi, i, cutoff, sigma !,iLat0,iphi0,iPhiG
		real (kind=QR_K) :: cosLat,cosLat0,sinLat,sinLat0,lat0,phi0
		real (kind=QR_K) :: cosPhi,cosPhi0,sinPhi,sinPhi0,dum,cosDPhi
		real (kind=QR_K) :: sigma2,d2,d,cosd,coef
		real (kind=QR_K) :: totalStrength, dumDiv, term, totalDivStrength
		integer :: CHUNK, i1, i2, threadID, nThreads
		integer :: kooki
		integer :: iL, iP, iP1, iS, j1, j2, jj, izone

		write(*,*)'       markActiveGridCells()';call markActiveGridCells()
		write(*,*)'       reorderOnGrid()';call reorderOnGrid()

		if(.NOT.allocated(QOnGrid))allocate(QOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		QOnGrid=zero
		if(.NOT.allocated(hsPrevOnGrid))allocate(hsPrevOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		hsPrevOnGrid=zero

		cutoff = 4 !stop at cutoff*sigma
		sigma = 2! sigma = 2.0*spacing, spacing = dtheta=dphi
		sigma2=(2.0*sphericalGrid%dLat)**2;
		coef = 1.0/(pi*sigma2)
		totalStrength = 0.0
		kooki=0

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nActiveCells/nThreads
		!		write(*,*)'CHUNK',CHUNK
		i1=1; i2=nActiveCells

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &coef,nParticles,cutoff,sigma,sigma2,QOnGrid,nActiveCells,activeCells,&
		!$OMP &zones, nLatZones,nZones,hsPrevOnGrid)

		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nActiveCells

		do i=i1,i2
			!			write(*,*)threadID,i,i1,i2

			iLat = activeCells(i,1)
			iPhi = activeCells(i,2)

			cosLat=sphericalGrid%cosLatC(iLat)
			sinLat=sphericalGrid%sinLatC(iLat)
			cosPhi=sphericalGrid%cosPhiC(iPhi)
			sinPhi=sphericalGrid%sinPhiC(iPhi)
			dum = zero
			dumDiv = zero
			!if(1.eq.2)then
			do iL=max(iLat-sigma*cutoff,1),min(iLat+sigma*cutoff,sphericalGrid%nLat)
				do iP1=iPhi-sigma*cutoff,iPhi+sigma*cutoff
					iP=iP1
					if(iP1.lt.1)iP=iP1+sphericalGrid%nPhi
					if(iP1.gt.sphericalGrid%nPhi)iP=iP1-sphericalGrid%nPhi
					izone=(iP-1)*nLatZones + iL
					!					if(izone.lt.1.or.izone.gt.nZones)stop "dfdf"
					if(zones(izone)%nEls.gt.0) then
						j1 = zones(izone)%startLoc
						j2 = j1+zones(izone)%nEls-1
						do jj = j1, j2 !sources
							iS = particles(jj)%el(1)%index
							Lat0 = particles(iS)%el(1)%Lat
							cosLat0=particles(iS)%el(1)%cosLat
							sinLat0=particles(iS)%el(1)%sinLat
							phi0 = particles(iS)%el(1)%Phi
							cosPhi0=particles(iS)%el(1)%cosPhi
							sinPhi0=particles(iS)%el(1)%sinPhi
							cosDPhi = cosPhi*cosPhi0+sinPhi*sinPhi0 !
							cosd=sinLat*sinLat0+cosLat*cosLat0*cosDPhi
							if(cosd.gt.1.0)then
								cosd=1.0
							else if(cosd.lt.-1.0) then
								cosd=-1.0
							end if
							d = acos(cosd)
							if(d.lt.0.0)d=d+two*pi
							d2=d*d/sigma2
							term = coef*exp(-d2)
							dum = dum + term*particles(iS)%el(1)%strengthVor
							dumDiv = dumDiv + term*particles(iS)%el(1)%hs
						end do
					end if
				end do
			end do

			QOnGrid(iLat,iPhi) = dum
			hsPrevOnGrid(iLat,iPhi) = dumDiv
			!end if
		end do
		!$OMP END PARALLEL

		totalStrength = zero
		totalDivStrength = zero
		do i=1,nActiveCells
			iLat = activeCells(i,1)
			iPhi = activeCells(i,2)
			totalStrength = totalStrength+QOnGrid(iLat,iPhi)*sphericalGrid%dA(iLat)
			totalDivStrength = totalDivStrength+hsPrevOnGrid(iLat,iPhi)*sphericalGrid%dA(iLat)
		end do
		write(*,*)'totalVorStrength and totalDivStrength after interpolating from particles to grid =',totalStrength, totalDivStrength

		if(allocated(activeCells))deallocate(activeCells)

		return
	end subroutine interpolate_Q_hs_FromParticlesToGrid

	subroutine interpolate_hs_FromParticlesToGrid()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		use QuickSearchVariables
		!$ use omp_lib
		implicit none
		integer :: iLat, iphi, i, cutoff, sigma !,iLat0,iphi0,iPhiG
		real (kind=QR_K) :: cosLat,cosLat0,sinLat,sinLat0,lat0,phi0
		real (kind=QR_K) :: cosPhi,cosPhi0,sinPhi,sinPhi0,dum,cosDPhi
		real (kind=QR_K) :: sigma2,d2,d,cosd,coef
		real (kind=QR_K) :: term
		integer :: CHUNK, i1, i2, threadID, nThreads
		integer :: iL, iP, iP1, iS, j1, j2, jj, izone

		if(.NOT.allocated(hsPrevOnGrid))allocate(hsPrevOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		hsPrevOnGrid=zero

		cutoff = 4 !stop at cutoff*sigma
		sigma = 2! sigma = 2.0*spacing, spacing = dtheta=dphi
		sigma2=(2.0*sphericalGrid%dLat)**2;
		coef = 1.0/(pi*sigma2)

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nActiveCells/nThreads
		!		write(*,*)'CHUNK',CHUNK
		i1=1; i2=nActiveCells

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &coef,nParticles,cutoff,sigma,sigma2,nActiveCells,activeCells,&
		!$OMP &zones, nLatZones,nZones,hsPrevOnGrid)

		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nActiveCells

		do i=i1,i2
			!			write(*,*)threadID,i,i1,i2

			iLat = activeCells(i,1)
			iPhi = activeCells(i,2)

			cosLat=sphericalGrid%cosLatC(iLat)
			sinLat=sphericalGrid%sinLatC(iLat)
			cosPhi=sphericalGrid%cosPhiC(iPhi)
			sinPhi=sphericalGrid%sinPhiC(iPhi)
			dum = zero
			do iL=max(iLat-sigma*cutoff,1),min(iLat+sigma*cutoff,sphericalGrid%nLat)
				do iP1=iPhi-sigma*cutoff,iPhi+sigma*cutoff
					iP=iP1
					if(iP1.lt.1)iP=iP1+sphericalGrid%nPhi
					if(iP1.gt.sphericalGrid%nPhi)iP=iP1-sphericalGrid%nPhi
					izone=(iP-1)*nLatZones + iL
					if(zones(izone)%nEls.gt.0) then
						j1 = zones(izone)%startLoc
						j2 = j1+zones(izone)%nEls-1
						do jj = j1, j2 !sources
							iS = particles(jj)%el(1)%index
							Lat0 = particles(iS)%el(1)%Lat
							cosLat0=particles(iS)%el(1)%cosLat
							sinLat0=particles(iS)%el(1)%sinLat
							phi0 = particles(iS)%el(1)%Phi
							cosPhi0=particles(iS)%el(1)%cosPhi
							sinPhi0=particles(iS)%el(1)%sinPhi
							cosDPhi = cosPhi*cosPhi0+sinPhi*sinPhi0 !
							cosd=sinLat*sinLat0+cosLat*cosLat0*cosDPhi
							if(cosd.gt.1.0)then
								cosd=1.0
							else if(cosd.lt.-1.0) then
								cosd=-1.0
							end if
							d = acos(cosd)
							if(d.lt.0.0)d=d+two*pi
							d2=d*d/sigma2
							term = coef*exp(-d2)
							dum = dum + term*particles(iS)%el(1)%hs
						end do
					end if
				end do
			end do
			hsPrevOnGrid(iLat,iPhi) = dum
		end do
		!$OMP END PARALLEL

		if(allocated(activeCells))deallocate(activeCells)

		return
	end subroutine interpolate_hs_FromParticlesToGrid

	subroutine interpolate_Q_FromParticlesToGrid()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		use QuickSearchVariables
		!$ use omp_lib
		implicit none
		integer :: iLat, iphi, i, cutoff, sigma !,iLat0,iphi0,iPhiG
		real (kind=QR_K) :: cosLat,cosLat0,sinLat,sinLat0,lat0,phi0
		real (kind=QR_K) :: cosPhi,cosPhi0,sinPhi,sinPhi0,dum,cosDPhi
		real (kind=QR_K) :: sigma2,d2,d,cosd,coef
		real (kind=QR_K) :: totalStrength, term
		integer :: CHUNK, i1, i2, threadID, nThreads
		integer :: iL, iP, iP1, iS, j1, j2, jj, izone

		write(*,*)'       markActiveGridCells()';call markActiveGridCells()
		write(*,*)'       reorderOnGrid()';call reorderOnGrid()

		if(.NOT.allocated(QOnGrid))allocate(QOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))

		QOnGrid=zero

		cutoff = 4 !stop at cutoff*sigma
		sigma  = 2! sigma = 2.0*spacing, spacing = dtheta=dphi
		sigma2 = (2.0*sphericalGrid%dLat)**2;
		coef   = 1.0/(pi*sigma2)
		totalStrength = 0.0

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nActiveCells/nThreads
		i1=1; i2=nActiveCells

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &coef,nParticles,cutoff,sigma,sigma2,QOnGrid,nActiveCells,activeCells,&
		!$OMP &zones, nLatZones,nZones)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nActiveCells
		do i=i1,i2
			iLat = activeCells(i,1)
			iPhi = activeCells(i,2)
			cosLat=sphericalGrid%cosLatC(iLat)
			sinLat=sphericalGrid%sinLatC(iLat)
			cosPhi=sphericalGrid%cosPhiC(iPhi)
			sinPhi=sphericalGrid%sinPhiC(iPhi)
			dum = zero

			do iL=max(iLat-sigma*cutoff,1),min(iLat+sigma*cutoff,sphericalGrid%nLat)
				do iP1=iPhi-sigma*cutoff,iPhi+sigma*cutoff
					iP=iP1
					if(iP1.lt.1)iP=iP1+sphericalGrid%nPhi
					if(iP1.gt.sphericalGrid%nPhi)iP=iP1-sphericalGrid%nPhi
					izone=(iP-1)*nLatZones + iL
					if(zones(izone)%nEls.gt.0) then
						j1 = zones(izone)%startLoc
						j2 = j1+zones(izone)%nEls-1
						do jj = j1, j2 !sources
							iS = particles(jj)%el(1)%index
							Lat0 = particles(iS)%el(1)%Lat
							cosLat0=particles(iS)%el(1)%cosLat
							sinLat0=particles(iS)%el(1)%sinLat
							phi0 = particles(iS)%el(1)%Phi
							cosPhi0=particles(iS)%el(1)%cosPhi
							sinPhi0=particles(iS)%el(1)%sinPhi
							cosDPhi = cosPhi*cosPhi0+sinPhi*sinPhi0 !
							cosd=sinLat*sinLat0+cosLat*cosLat0*cosDPhi
							if(cosd.gt.1.0)then
								cosd=1.0
							else if(cosd.lt.-1.0) then
								cosd=-1.0
							end if
							d = acos(cosd)
							if(d.lt.0.0)d=d+two*pi
							d2=d*d/sigma2
							term = coef*exp(-d2)
							dum = dum + term*particles(iS)%el(1)%strengthVor
						end do
					end if
				end do
			end do
			QOnGrid(iLat,iPhi) = dum
			
			! write(*,*) QOnGrid(iLat,iPhi)

		end do
		!$OMP END PARALLEL

		totalStrength = zero
		do i=1,nActiveCells
			iLat = activeCells(i,1)
			iPhi = activeCells(i,2)
			totalStrength = totalStrength+QOnGrid(iLat,iPhi)*sphericalGrid%dA(iLat)
		end do
		write(*,*)'totalVorStrength after interpolating from particles to grid =',totalStrength

		if(allocated(activeCells))deallocate(activeCells)

		return
	end subroutine interpolate_Q_FromParticlesToGrid

	subroutine setQonGridToZero()
		use sphericalGridDataModule
		implicit none
		QOnGrid=0.0D0
		return
	end subroutine setQonGridToZero

	subroutine interpolate_D_FromParticlesToGrid()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		use QuickSearchVariables
		!$ use omp_lib
		implicit none
		integer :: iLat, iphi, i, cutoff, sigma !,iLat0,iphi0,iPhiG
		real (kind=QR_K) :: cosLat,cosLat0,sinLat,sinLat0,lat0,phi0
		real (kind=QR_K) :: cosPhi,cosPhi0,sinPhi,sinPhi0,cosDPhi
		real (kind=QR_K) :: sigma2,d2,d,cosd,coef
		real (kind=QR_K) :: dumDiv, term, totalDivStrength
		integer :: CHUNK, i1, i2, threadID, nThreads
		integer :: iL, iP, iP1, iS, j1, j2, jj, izone

		write(*,*)'       markActiveGridCells()';call markActiveGridCells()
		write(*,*)'       reorderOnGrid()';call reorderOnGrid()
		if(.NOT.allocated(divOnGrid))allocate(divOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		divOnGrid=zero

		cutoff = 4 !stop at cutoff*sigma
		sigma = 2! sigma = 2.0*spacing, spacing = dtheta=dphi
		sigma2=(2.0*sphericalGrid%dLat)**2;
		coef = 1.0/(pi*sigma2)

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nActiveCells/nThreads
		i1=1; i2=nActiveCells

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &coef,nParticles,cutoff,sigma,sigma2,nActiveCells,activeCells,&
		!$OMP &zones, nLatZones,nZones,divOnGrid)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nActiveCells
		do i=i1,i2
			iLat = activeCells(i,1)
			iPhi = activeCells(i,2)
			cosLat=sphericalGrid%cosLatC(iLat)
			sinLat=sphericalGrid%sinLatC(iLat)
			cosPhi=sphericalGrid%cosPhiC(iPhi)
			sinPhi=sphericalGrid%sinPhiC(iPhi)
			dumDiv = zero
			do iL=max(iLat-sigma*cutoff,1),min(iLat+sigma*cutoff,sphericalGrid%nLat)
				do iP1=iPhi-sigma*cutoff,iPhi+sigma*cutoff
					iP=iP1
					if(iP1.lt.1)iP=iP1+sphericalGrid%nPhi
					if(iP1.gt.sphericalGrid%nPhi)iP=iP1-sphericalGrid%nPhi
					izone=(iP-1)*nLatZones + iL
					if(zones(izone)%nEls.gt.0) then
						j1 = zones(izone)%startLoc
						j2 = j1+zones(izone)%nEls-1
						do jj = j1, j2 !sources
							iS = particles(jj)%el(1)%index
							Lat0 = particles(iS)%el(1)%Lat
							cosLat0=particles(iS)%el(1)%cosLat
							sinLat0=particles(iS)%el(1)%sinLat
							phi0 = particles(iS)%el(1)%Phi
							cosPhi0=particles(iS)%el(1)%cosPhi
							sinPhi0=particles(iS)%el(1)%sinPhi
							cosDPhi = cosPhi*cosPhi0+sinPhi*sinPhi0 !
							cosd=sinLat*sinLat0+cosLat*cosLat0*cosDPhi
							if(cosd.gt.1.0)then
								cosd=1.0
							else if(cosd.lt.-1.0) then
								cosd=-1.0
							end if
							d = acos(cosd)
							if(d.lt.0.0)d=d+two*pi
							d2=d*d/sigma2
							term = coef*exp(-d2)
							dumDiv = dumDiv + term*particles(iS)%el(1)%strengthDiv
						end do
					end if
				end do
			end do
			divOnGrid(iLat,iPhi) = dumDiv
		end do
		!$OMP END PARALLEL

		totalDivStrength = zero
		do i=1,nActiveCells
			iLat = activeCells(i,1)
			iPhi = activeCells(i,2)
			totalDivStrength = totalDivStrength+divOnGrid(iLat,iPhi)*sphericalGrid%dA(iLat)
		end do
		write(*,*)'totalDivStrength after interpolating from particles to grid =',totalDivStrength

		!if(allocated(activeCells))deallocate(activeCells)

		return
	end subroutine interpolate_D_FromParticlesToGrid

	subroutine regrid()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		!$ use omp_lib
		implicit none
		integer :: CHUNK, i1, i2, threadID, nThreads
		real (kind=QR_K) :: Lat,Phi,dA,colatitude,totalStrength
		integer :: i, iPhi,iLat

		totalStrength = 0.0
		do i=1,nParticles
			if(allocated(particles(i)%el))deallocate(particles(i)%el)
		end do
		nParticles = 0

		i1=1; i2=sphericalGrid%nPhi

		!		nThreads = OMP_get_max_threads()
		!		call OMP_SET_NUM_THREADS(nThreads)
		!		CHUNK = sphericalGrid%nPhi/nThreads
		!		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!		!$OMP &nParticles,QOnGrid,velocityVorOnGrid)
		!		threadID = omp_get_thread_num()
		!		i1 = threadID*CHUNK+1
		!		i2=i1-1+CHUNK
		!		if(threadID+1.eq.nThreads)i2=sphericalGrid%nPhi

		do iPhi=i1,i2 !1,nParticles
			Phi = half*sphericalGrid%dphi+sphericalGrid%dphi*real(iPhi-1,QR_K)
			do iLat=1,sphericalGrid%nLat
				Lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
				if(ABS(QOnGrid(iLat,iPhi)/OMEGA).gt.1.0D-10.and.sphericalGrid%d(iLat,iPhi).gt.zero)then
					dA=sphericalGrid%dA(iLat)
					!					!$OMP CRITICAL
					nParticles=nParticles+1
					!					write(*,*)nParticles
					allocate(particles(nParticles)%el(1))
					particles(nParticles)%el(1)%lat = Lat
					particles(nParticles)%el(1)%Phi = phi
					particles(nParticles)%el(1)%strengthVor=QOnGrid(iLat,iPhi)*dA
					particles(nParticles)%el(1)%uLat = velocityOnGrid(1,iLat,iPhi)
					particles(nParticles)%el(1)%uLon = velocityOnGrid(2,iLat,iPhi)
					particles(nParticles)%el(1)%sinLat=sphericalGrid%sinLatC(iLat)
					particles(nParticles)%el(1)%cosLat=sphericalGrid%cosLatC(iLat)
					particles(nParticles)%el(1)%sinPhi=sphericalGrid%sinPhiC(iPhi)
					particles(nParticles)%el(1)%cosPhi=sphericalGrid%cosPhiC(iPhi)
					totalStrength = totalStrength+particles(nParticles)%el(1)%strengthVor
					particles(nParticles)%el(1)%d = sphericalGrid%d(iLat,iPhi)
					particles(nParticles)%el(1)%hPrev = hsPrevOnGrid(iLat,iPhi)+particles(nParticles)%el(1)%d
				end if
			end do
		end do
		!		!$OMP END PARALLEL

		write(*,*)'totalStrength after regridding = ',totalStrength

		if(allocated(QOnGrid))deallocate(QOnGrid)
		if(allocated(hsPrevOnGrid))deallocate(hsPrevOnGrid)

		return
	end subroutine regrid

	subroutine interpolate_V_FromGridToParticles()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		!$ use omp_lib
		implicit none
		integer :: CHUNK, i1, i2, threadID, nThreads
		real (kind=QR_K) :: Lat0,phi0,LatB,LatT,phiL,phiR, beta
		real (kind=QR_K) :: alpha,V1phiL,V1phiR,V2phiR,V2phiL
		integer :: i,iLatB,iLatT,iphiL,iphiR
		real (kind=QR_K) :: term, VBL, VTL, VBR, VTR

		i1=1; i2=nParticles

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nParticles/nThreads
		write(*,*)CHUNK
		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &nParticles,velocityOnGrid,psiOnGrid)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nParticles

		do i=i1,i2 !1,nParticles
			Lat0 = particles(i)%el(1)%Lat
			phi0 = particles(i)%el(1)%Phi
			iLatB = 1+int((Lat0+half*pi-half*sphericalGrid%dLat)/sphericalGrid%dLat)
			iLatT = iLatB+1

			iphiL = 1+int((phi0-half*sphericalGrid%dphi)/sphericalGrid%dphi)
			iphiR = iphiL+1; if(iphiR.gt.sphericalGrid%nPhi)iphiR=1

			LatB=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatB-1,QR_K)
			LatT=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatT-1,QR_K)
			phiL= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiL-1,QR_K)
			phiR= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiR-1,QR_K)

			alpha = (Lat0-LatB)/(LatT-LatB)
			V1phiL = velocityOnGrid(1,iLatB,iphiL) + alpha*(velocityOnGrid(1,iLatT,iphiL)-velocityOnGrid(1,iLatB,iphiL))
			V2phiL = velocityOnGrid(2,iLatB,iphiL) + alpha*(velocityOnGrid(2,iLatT,iphiL)-velocityOnGrid(2,iLatB,iphiL))
			V1phiR = velocityOnGrid(1,iLatB,iphiR) + alpha*(velocityOnGrid(1,iLatT,iphiR)-velocityOnGrid(1,iLatB,iphiR))
			V2phiR = velocityOnGrid(2,iLatB,iphiR) + alpha*(velocityOnGrid(2,iLatT,iphiR)-velocityOnGrid(2,iLatB,iphiR))
			beta = (phi0-phiL)/(phiR-phiL)
			particles(i)%el(1)%uLat = V1phiL + beta*(V1phiR-V1phiL)
			particles(i)%el(1)%uLon = V2phiL + beta*(V2phiR-V2phiL)

			!term = (two*OMEGA*particles(i)%el(1)%sinLat/GRAVITY)
			!VBL=sphericalGrid%d(iLatB,iphiL)+term*psiOnGrid(iLatB,iphiL)
			!VTL=sphericalGrid%d(iLatT,iphiL)+term*psiOnGrid(iLatT,iphiL)
			!VBR=sphericalGrid%d(iLatB,iphiR)+term*psiOnGrid(iLatB,iphiR)
			!VTR=sphericalGrid%d(iLatT,iphiR)+term*psiOnGrid(iLatT,iphiR)

			!V1phiL = VBL + alpha*(VTL-VBL)
			!V1phiR = VBR + alpha*(VTR-VBR)
			!particles(i)%el(1)%h = V1phiL + beta*(V1phiR-V1phiL)
		end do
		!$OMP END PARALLEL

		return
	end subroutine interpolate_V_FromGridToParticles

	subroutine initialize_h_d_hPrev_AtParticlesLocations()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		!$ use omp_lib
		implicit none
		integer :: CHUNK, i1, i2, threadID, nThreads
		real (kind=QR_K) :: Lat0,phi0,LatB,LatT,phiL,phiR
		real (kind=QR_K) :: alpha,V1phiL,V1phiR,V2phiR,V2phiL,beta
		integer :: i,iLatB,iLatT,iphiL,iphiR
		real (kind=QR_K) :: term, VBL, VTL, VBR, VTR, hPrev, dt

		i1=1; i2=nParticles

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nParticles/nThreads
		write(*,*)CHUNK
		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &nParticles,psiOnGrid)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nParticles

		do i=i1,i2 !1,nParticles
			Lat0 = particles(i)%el(1)%Lat
			phi0 = particles(i)%el(1)%Phi
			iLatB = 1+int((Lat0+half*pi-half*sphericalGrid%dLat)/sphericalGrid%dLat)
			iLatT = iLatB+1

			iphiL = 1+int((phi0-half*sphericalGrid%dphi)/sphericalGrid%dphi)
			iphiR = iphiL+1; if(iphiR.gt.sphericalGrid%nPhi)iphiR=1

			LatB=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatB-1,QR_K)
			LatT=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatT-1,QR_K)
			phiL= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiL-1,QR_K)
			phiR= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiR-1,QR_K)

			term = (two*OMEGA*particles(i)%el(1)%sinLat/GRAVITY)
			VBL=sphericalGrid%d(iLatB,iphiL)+term*psiOnGrid(iLatB,iphiL)
			VTL=sphericalGrid%d(iLatT,iphiL)+term*psiOnGrid(iLatT,iphiL)
			VBR=sphericalGrid%d(iLatB,iphiR)+term*psiOnGrid(iLatB,iphiR)
			VTR=sphericalGrid%d(iLatT,iphiR)+term*psiOnGrid(iLatT,iphiR)
			alpha = (Lat0-LatB)/(LatT-LatB)
			beta = (phi0-phiL)/(phiR-phiL)
			V1phiL = VBL + alpha*(VTL-VBL)
			V1phiR = VBR + alpha*(VTR-VBR)
			particles(i)%el(1)%h = V1phiL + beta*(V1phiR-V1phiL)
			particles(i)%el(1)%hPrev=particles(i)%el(1)%h
			particles(i)%el(1)%strengthDiv = zero

			term=zero
			VBL=sphericalGrid%d(iLatB,iphiL)+term*psiOnGrid(iLatB,iphiL)
			VTL=sphericalGrid%d(iLatT,iphiL)+term*psiOnGrid(iLatT,iphiL)
			VBR=sphericalGrid%d(iLatB,iphiR)+term*psiOnGrid(iLatB,iphiR)
			VTR=sphericalGrid%d(iLatT,iphiR)+term*psiOnGrid(iLatT,iphiR)
			alpha = (Lat0-LatB)/(LatT-LatB)
			beta = (phi0-phiL)/(phiR-phiL)
			V1phiL = VBL + alpha*(VTL-VBL)
			V1phiR = VBR + alpha*(VTR-VBR)
			particles(i)%el(1)%d = V1phiL + beta*(V1phiR-V1phiL)
			!particles(i)%el(1)%h = particles(i)%el(1)%d
		end do
		!$OMP END PARALLEL

		!		if(.NOT.allocated(velocityDivOnGrid))allocate(velocityDivOnGrid(2,sphericalGrid%nLat,sphericalGrid%nPhi))
		!		velocityDivOnGrid=zero
		!velocityOnGrid = velocityVorOnGrid

		return
	end subroutine initialize_h_d_hPrev_AtParticlesLocations 

	subroutine compute_d_AtParticlesLocations()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		!$ use omp_lib
		implicit none
		integer :: CHUNK, i1, i2, threadID, nThreads
		real (kind=QR_K) :: Lat0,phi0,LatB,LatT,phiL,phiR
		real (kind=QR_K) :: alpha,V1phiL,V1phiR,V2phiR,V2phiL,beta
		integer :: i,iLatB,iLatT,iphiL,iphiR
		real (kind=QR_K) :: term, VBL, VTL, VBR, VTR, hPrev, dt

		i1=1; i2=nParticles

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nParticles/nThreads
		write(*,*)CHUNK
		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &nParticles,psiOnGrid)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nParticles

		do i=i1,i2 !1,nParticles
			Lat0 = particles(i)%el(1)%Lat
			phi0 = particles(i)%el(1)%Phi
			iLatB = 1+int((Lat0+half*pi-half*sphericalGrid%dLat)/sphericalGrid%dLat)
			iLatT = iLatB+1

			iphiL = 1+int((phi0-half*sphericalGrid%dphi)/sphericalGrid%dphi)
			iphiR = iphiL+1; if(iphiR.gt.sphericalGrid%nPhi)iphiR=1

			LatB=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatB-1,QR_K)
			LatT=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatT-1,QR_K)
			phiL= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiL-1,QR_K)
			phiR= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiR-1,QR_K)

			term = (two*OMEGA*particles(i)%el(1)%sinLat/GRAVITY)
			VBL=sphericalGrid%d(iLatB,iphiL)+term*psiOnGrid(iLatB,iphiL)
			VTL=sphericalGrid%d(iLatT,iphiL)+term*psiOnGrid(iLatT,iphiL)
			VBR=sphericalGrid%d(iLatB,iphiR)+term*psiOnGrid(iLatB,iphiR)
			VTR=sphericalGrid%d(iLatT,iphiR)+term*psiOnGrid(iLatT,iphiR)
			alpha = (Lat0-LatB)/(LatT-LatB)
			beta = (phi0-phiL)/(phiR-phiL)
			V1phiL = VBL + alpha*(VTL-VBL)
			V1phiR = VBR + alpha*(VTR-VBR)
			particles(i)%el(1)%h = V1phiL + beta*(V1phiR-V1phiL)
			particles(i)%el(1)%hPrev=particles(i)%el(1)%h
			particles(i)%el(1)%strengthDiv = zero

			term=zero
			VBL=sphericalGrid%d(iLatB,iphiL)+term*psiOnGrid(iLatB,iphiL)
			VTL=sphericalGrid%d(iLatT,iphiL)+term*psiOnGrid(iLatT,iphiL)
			VBR=sphericalGrid%d(iLatB,iphiR)+term*psiOnGrid(iLatB,iphiR)
			VTR=sphericalGrid%d(iLatT,iphiR)+term*psiOnGrid(iLatT,iphiR)
			alpha = (Lat0-LatB)/(LatT-LatB)
			beta = (phi0-phiL)/(phiR-phiL)
			V1phiL = VBL + alpha*(VTL-VBL)
			V1phiR = VBR + alpha*(VTR-VBR)
			particles(i)%el(1)%d = V1phiL + beta*(V1phiR-V1phiL)
		end do
		!$OMP END PARALLEL

		return
	end subroutine compute_d_AtParticlesLocations

	! subroutine computeh0hAtParticlesLocations()
		! 	use sphericalGridDataModule
		! 	use basicDataStructures
		! 	use simulationParameters
		! 	use particlesDataModule
		! 	!$ use omp_lib
		! 	implicit none
		! 	integer :: CHUNK, i1, i2, threadID, nThreads
		! 	real (kind=QR_K) :: Lat0,phi0,LatB,LatT,phiL,phiR
		! 	real (kind=QR_K) :: alpha,V1phiL,V1phiR,V2phiR,V2phiL,beta
		! 	integer :: i,iLatB,iLatT,iphiL,iphiR
		! 	real (kind=QR_K) :: term, VBL, VTL, VBR, VTR, hPrev, dt,h,h0,hs
		!
		! 	i1=1; i2=nParticles
		!
		! 	nThreads = OMP_get_max_threads()
		! 	call OMP_SET_NUM_THREADS(nThreads)
		! 	CHUNK = nParticles/nThreads
		! 	write(*,*)CHUNK
		! 	!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		! 	!$OMP &nParticles)
		! 	threadID = omp_get_thread_num()
		! 	i1 = threadID*CHUNK+1
		! 	i2=i1-1+CHUNK
		! 	if(threadID+1.eq.nThreads)i2=nParticles
		!
		! 	do i=i1,i2 !1,nParticles
		! 		Lat0 = particles(i)%el(1)%Lat
		! 		phi0 = particles(i)%el(1)%Phi
		! 		iLatB = 1+int((Lat0+half*pi-half*sphericalGrid%dLat)/sphericalGrid%dLat)
		! 		iLatT = iLatB+1
		!
		! 		iphiL = 1+int((phi0-half*sphericalGrid%dphi)/sphericalGrid%dphi)
		! 		iphiR = iphiL+1; if(iphiR.gt.sphericalGrid%nPhi)iphiR=1
		!
		! 		LatB=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatB-1,QR_K)
		! 		LatT=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatT-1,QR_K)
		! 		phiL= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiL-1,QR_K)
		! 		phiR= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiR-1,QR_K)
		!
		! 		alpha = (Lat0-LatB)/(LatT-LatB)
		! 		beta = (phi0-phiL)/(phiR-phiL)
		!
		! 		VBL=sphericalGrid%d(iLatB,iphiL)
		! 		VTL=sphericalGrid%d(iLatT,iphiL)
		! 		VBR=sphericalGrid%d(iLatB,iphiR)
		! 		VTR=sphericalGrid%d(iLatT,iphiR)
		! 		V1phiL = VBL + alpha*(VTL-VBL)
		! 		V1phiR = VBR + alpha*(VTR-VBR)
		! 		h0 = V1phiL + beta*(V1phiR-V1phiL)
		!
		! 		term = (two*OMEGA*particles(i)%el(1)%sinLat/GRAVITY)
		! 		VBL=term*psiOnGrid(iLatB,iphiL)
		! 		VTL=term*psiOnGrid(iLatT,iphiL)
		! 		VBR=term*psiOnGrid(iLatB,iphiR)
		! 		VTR=term*psiOnGrid(iLatT,iphiR)
		! 		V1phiL = VBL + alpha*(VTL-VBL)
		! 		V1phiR = VBR + alpha*(VTR-VBR)
		! 		hs = V1phiL + beta*(V1phiR-V1phiL)
		!
		! 		particles(i)%el(1)%d = h0
		! 		particles(i)%el(1)%h = h0+hs
		! 	end do
		! 	!$OMP END PARALLEL
		!
		! 	return
	! end subroutine computeh0hAtParticlesLocations

	subroutine compute_d_h_h_s_div_AtParticlesLocations(step,iTime)
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use particlesDataModule
		!$ use omp_lib
		implicit none
		integer, intent(in) :: step, iTime
		integer :: CHUNK, i1, i2, threadID, nThreads
		real (kind=QR_K) :: Lat0,phi0,LatB,LatT,phiL,phiR
		real (kind=QR_K) :: alpha,V1phiL,V1phiR,V2phiR,V2phiL,beta
		integer :: i,iLatB,iLatT,iphiL,iphiR
		real (kind=QR_K) :: term, VBL, VTL, VBR, VTR, hPrev, dt,h,d,hs

		i1=1; i2=nParticles

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = nParticles/nThreads
		write(*,*)CHUNK
		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,particles,sphericalGrid,&
		!$OMP &nParticles,psiOnGrid,step,iTime)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=nParticles

		do i=i1,i2 !1,nParticles
			Lat0 = particles(i)%el(1)%Lat
			phi0 = particles(i)%el(1)%Phi
			iLatB = 1+int((Lat0+half*pi-half*sphericalGrid%dLat)/sphericalGrid%dLat)
			iLatT = iLatB+1

			iphiL = 1+int((phi0-half*sphericalGrid%dphi)/sphericalGrid%dphi)
			iphiR = iphiL+1; if(iphiR.gt.sphericalGrid%nPhi)iphiR=1

			LatB=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatB-1,QR_K)
			LatT=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLatT-1,QR_K)
			phiL= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiL-1,QR_K)
			phiR= half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphiR-1,QR_K)

			alpha = (Lat0-LatB)/(LatT-LatB)
			beta = (phi0-phiL)/(phiR-phiL)

			VBL=sphericalGrid%d(iLatB,iphiL)
			VTL=sphericalGrid%d(iLatT,iphiL)
			VBR=sphericalGrid%d(iLatB,iphiR)
			VTR=sphericalGrid%d(iLatT,iphiR)
			V1phiL = VBL + alpha*(VTL-VBL)
			V1phiR = VBR + alpha*(VTR-VBR)
			d = V1phiL + beta*(V1phiR-V1phiL)

			if(include_hs)then
				term = (two*OMEGA*particles(i)%el(1)%sinLat/GRAVITY)
				VBL=term*psiOnGrid(iLatB,iphiL)
				VTL=term*psiOnGrid(iLatT,iphiL)
				VBR=term*psiOnGrid(iLatB,iphiR)
				VTR=term*psiOnGrid(iLatT,iphiR)
				V1phiL = VBL + alpha*(VTL-VBL)
				V1phiR = VBR + alpha*(VTR-VBR)
				hs = V1phiL + beta*(V1phiR-V1phiL)
			else
				hs = zero
			end if

			if(step.eq.1)then
				hPrev = particles(i)%el(1)%hPrev
				dt = deltat
			else
				hPrev = particles(i)%el(1)%h
				dt = half*deltat
			end if
			particles(i)%el(1)%d = d
			particles(i)%el(1)%h = d+hs
			particles(i)%el(1)%hs = hs

			if(iTime.eq.1.and.step.eq.1)then
				particles(i)%el(1)%strengthDiv = zero
			else
				particles(i)%el(1)%strengthDiv = -(log(max(ten,particles(i)%el(1)%h)) - log(max(ten,hPrev)))/dt
			end if
			!if(ISNAN(particles(i)%el(1)%strengthDiv))then
			!	write(*,*)i,particles(i)%el(1)%h,hPrev
			!	stop
			!end if

			if(step.eq.1)particles(i)%el(1)%hPrev=particles(i)%el(1)%h
		end do
		!$OMP END PARALLEL

		return
	end subroutine compute_d_h_h_s_div_AtParticlesLocations

	subroutine interpolateBathymetryToGrid()
		use sphericalGridDataModule
		use basicDataStructures
		use simulationParameters
		use bathymetryDataModule
		!$ use omp_lib
		implicit none
		integer :: CHUNK, i1, i2, threadID, nThreads
		real (kind=QR_K) :: Lat0,phi0,LatB,LatT,phiL,phiR
		real (kind=QR_K) :: alpha,zphiL,zphiR
		integer :: i,iLatB,iLatT,iphiL,iphiR,iLat,iPhi,nnn
		real (kind=QR_K) :: dLatB, dPhiB, zT, zB, aveDepth

		if(allocated(sphericalGrid%d))deallocate(sphericalGrid%d)
		if(.NOT.allocated(sphericalGrid%d))allocate(sphericalGrid%d(sphericalGrid%nLat,sphericalGrid%nPhi))
		sphericalGrid%d=zero

		dLatB = pi/real(bathymetry%nLat,QR_K)
		dPhiB = two*pi/real(bathymetry%nLong,QR_K)

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = sphericalGrid%nLat/nThreads
		write(*,*)CHUNK
		i1=1; i2=sphericalGrid%nLat

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,bathymetry,sphericalGrid,&
		!$OMP &dLatB,dPhiB)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=sphericalGrid%nLat

		do iLat=i1,i2
			Lat0=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(ilat-1,QR_K)
			iLatB = 1+int((Lat0+half*pi-half*dLatB)/dLatB)
			iLatT = iLatB+1
			LatB=-half*pi+half*dLatB+dLatB*real(iLatB-1,QR_K)
			LatT=-half*pi+half*dLatB+dLatB*real(iLatT-1,QR_K)

			do iphi = 1,sphericalGrid%nPhi
				phi0=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi-1,QR_K)

				iphiL = 1+int((phi0-half*dphiB)/dphiB)
				iphiR = iphiL+1; if(iphiR.gt.bathymetry%nLong)iphiR=1
				phiL= half*dphiB+dphiB*real(iphiL-1,QR_K)
				phiR= half*dphiB+dphiB*real(iphiR-1,QR_K)

				alpha = (Lat0-LatB)/(LatT-LatB)
				zT = real(bathymetry%z(iphiL,iLatT),QR_K)
				zB = real(bathymetry%z(iphiL,iLatB),QR_K)
				zphiL = zB + alpha*(zT-zB)
				zT = real(bathymetry%z(iphiR,iLatT),QR_K)
				zB = real(bathymetry%z(iphiR,iLatB),QR_K)
				zphiR = zB + alpha*(zT-zB)
				alpha = (phi0-phiL)/(phiR-phiL)
				sphericalGrid%d(iLat,iPhi) = -(zphiL + alpha*(zphiR-zphiL))
			end do
		end do
		!$OMP END PARALLEL

		nnn=0
		aveDepth=zero
		do iLat=1,sphericalGrid%nLat
			do iphi = 1,sphericalGrid%nPhi
				if(sphericalGrid%d(iLat,iPhi).gt.zero) then
					aveDepth=(real(nnn,QR_K)*aveDepth+sphericalGrid%d(iLat,iPhi))/real(nnn+1,QR_K)
					nnn=nnn+1
				end if
			end do
		end do
		sphericalGrid%d0 = aveDepth

		Ld0 = sqrt(GRAVITY*sphericalGrid%d0)/(two*OMEGA)

		write(*,*)'*** Mean depth is calculated ***',sphericalGrid%d0

		write(*,*)'*** The Rossby Radius of Deformation is computed using, average depth: Ld0 = ***',Ld0

		return
	end subroutine interpolateBathymetryToGrid

!ADDED SUUBROUTINE: Interpolating the Bathymetry to the OceanGrid.
	subroutine interpolateBathymetryToOceanGrid()
		use OceanGridDataModule
		use basicDataStructures
		use simulationParameters
		use bathymetryDataModule
		!$ use omp_lib
		implicit none
		integer :: CHUNK, i1, i2, threadID, nThreads
		real (kind=QR_K) :: Lat0,phi0,LatB,LatT,phiL,phiR
		real (kind=QR_K) :: alpha,zphiL,zphiR,lat,sinLat,beta,aalpha
		integer :: i,iLatB,iLatT,iphiL,iphiR,iLat,iPhi,nnn
		real (kind=QR_K) :: dLatB, dPhiB, zT, zB, aveDepth
		! Mapping arrays for ocean points only
		integer, allocatable :: oceanPointIndex(:,:)  ! Maps (iLat,iLon) to ocean point index
		integer :: nOceanPoints,iLon

		if(allocated(OceanGrid%d))deallocate(OceanGrid%d)
		if(.NOT.allocated(OceanGrid%d))allocate(OceanGrid%d(OceanGrid%nLat,OceanGrid%nPhi))
		OceanGrid%d=zero

		dLatB = pi/real(bathymetry%nLat,QR_K)
		dPhiB = two*pi/real(bathymetry%nLong,QR_K)

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = OceanGrid%nLat/nThreads
		write(*,*) 'Chunks', CHUNK
		
		i1=1; i2=OceanGrid%nLat

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,bathymetry,OceanGrid,&
		!$OMP &dLatB,dPhiB)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=OceanGrid%nLat
		do iLat=i1,i2
			Lat0=-half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(ilat-1,QR_K)
			iLatB = 1+int((Lat0+half*pi-half*dLatB)/dLatB)
			iLatT = iLatB+1
			LatB=-half*pi+half*dLatB+dLatB*real(iLatB-1,QR_K)
			LatT=-half*pi+half*dLatB+dLatB*real(iLatT-1,QR_K)

			do iphi = 1,OceanGrid%nPhi
				phi0=half*OceanGrid%dphi+OceanGrid%dphi*real(iphi-1,QR_K)

				iphiL = 1+int((phi0-half*dphiB)/dphiB)
				iphiR = iphiL+1; if(iphiR.gt.bathymetry%nLong)iphiR=1
				phiL= half*dphiB+dphiB*real(iphiL-1,QR_K)
				phiR= half*dphiB+dphiB*real(iphiR-1,QR_K)

				alpha = (Lat0-LatB)/(LatT-LatB)
				zT = real(bathymetry%z(iphiL,iLatT),QR_K)
				zB = real(bathymetry%z(iphiL,iLatB),QR_K)
				zphiL = zB + alpha*(zT-zB)
				zT = real(bathymetry% z(iphiR,iLatT),QR_K)
				zB = real(bathymetry%z(iphiR,iLatB),QR_K)
				zphiR = zB + alpha*(zT-zB)
				alpha = (phi0-phiL)/(phiR-phiL)
				OceanGrid%d(iLat,iPhi) = -(zphiL + alpha*(zphiR-zphiL)) !!all points ocean and land
			end do
		end do
		!$OMP END PARALLEL
		!COMMENT: Check the values of positive depths
		nOceanPoints = 0
		do iLat=1,OceanGrid%nLat
			do iLon=1,OceanGrid%nPhi
				if (OceanGrid%d(iLat,iLon) > zero) then
					nOceanPoints = nOceanPoints + 1
				endif
			end do
		end do

		write(*,*)'Number of Only Ocean Particles = ', nOceanPoints
		write(*,*)'OceanGrid%d: Land and Ocean Particles =', size(OceanGrid%d)
		! write(*,*)'OceanGrid%d land and ocean=',OceanGrid%d	
		
		!COMMENT: Average Depth d₀ Calculations and Ld₀ 
		nnn=0
		aveDepth=zero
		do iLat=1,OceanGrid%nLat
			do iphi = 1,OceanGrid%nPhi
				if(OceanGrid%d(iLat,iPhi).gt.zero) then
					aveDepth=(real(nnn,QR_K)*aveDepth+OceanGrid%d(iLat,iPhi))/real(nnn+1,QR_K)
					nnn=nnn+1
				end if
			end do
		end do
		OceanGrid%d0 = aveDepth
		Ld0 = sqrt(GRAVITY*OceanGrid%d0)/(two*OMEGA)

		write(*,*)'⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯'
		write(*,*)'d Max , d Min =', maxval(OceanGrid%d) , ',', minval(OceanGrid%d)
		write(*,*)'Mean depth is computed, d₀ = ', OceanGrid%d0
		write(*,*)'Ocean particles, d>0 = ', nnn
		write(*,*)'The Rossby Radius of Deformation is computed Ld₀ = sqrt(g*d₀)/(2Ω) =  ', Ld0
		write(*,*)'⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯'
		return
	end subroutine interpolateBathymetryToOceanGrid

	!ADDED SUBROUTINE: Constructing the OceanGrid (θ,ƛ).
	subroutine constructOceanGrid()
		use OceanGridDataModule 
		use simulationParameters
		implicit none
		integer :: iphi,ilat
		real (kind=QR_K) :: lat, phi, colatitude

		OceanGrid%nPhi=2*OceanGrid%nLat 

		OceanGrid%dLat=(pi)/real(OceanGrid%nLat,QR_K)
		OceanGrid%dphi=(two*pi)/real(OceanGrid%nPhi,QR_K)

		write(*,*)' ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯ '
		write(*,*) 'Grid Spacing: dθ = ', OceanGrid%dLat , 'dƛ= ', OceanGrid%dphi
		write(*,*) 'Number of Grids: nθ , nƛ =',OceanGrid%nLat, OceanGrid%nPhi
		write(*,*) 'Total number of Grids: nθ x nƛ = ',OceanGrid%nLat * OceanGrid%nPhi
		write(*,*)' ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯ '
		write(*,*)
	
		if(.NOT.allocated(OceanGrid%sinLat))allocate(OceanGrid%sinLat(OceanGrid%nLat+1))
		if(.NOT.allocated(OceanGrid%cosLat))allocate(OceanGrid%cosLat(OceanGrid%nLat+1))
		if(.NOT.allocated(OceanGrid%sinLatC))allocate(OceanGrid%sinLatC(OceanGrid%nLat))
		if(.NOT.allocated(OceanGrid%cosLatC))allocate(OceanGrid%cosLatC(OceanGrid%nLat))

		if(.NOT.allocated(OceanGrid%sinPhi))allocate(OceanGrid%sinPhi(OceanGrid%nPhi))
		if(.NOT.allocated(OceanGrid%cosPhi))allocate(OceanGrid%cosPhi(OceanGrid%nPhi))
		if(.NOT.allocated(OceanGrid%sinPhiC))allocate(OceanGrid%sinPhiC(OceanGrid%nPhi))
		if(.NOT.allocated(OceanGrid%cosPhiC))allocate(OceanGrid%cosPhiC(OceanGrid%nPhi))

		!COMMENT: Cell Vertices
		do ilat=1,OceanGrid%nLat+1
			lat=-half*pi+OceanGrid%dLat*real(ilat-1,QR_K)
			OceanGrid%sinLat(ilat)=sin(lat)
			OceanGrid%cosLat(ilat)=cos(lat)
		end do
		do iphi=1,OceanGrid%nPhi
			phi=OceanGrid%dphi*real(iphi-1,QR_K)
			OceanGrid%sinPhi(iphi)=sin(phi)
			OceanGrid%cosPhi(iphi)=cos(phi)
		end do

		!COMMENT: Cell Centers
		do ilat=1,OceanGrid%nLat
			lat= -half*pi + half*OceanGrid%dLat + OceanGrid%dLat*real(ilat-1,QR_K)
			OceanGrid%sinLatC(ilat)=sin(lat)
			OceanGrid%cosLatC(ilat)=cos(lat)
			! write(*,*) 'Center: OceanGrid%dLat*real(ilat-1,QR_K) ,lat , sinLat = ',  OceanGrid%dLat*real(ilat-1,QR_K), lat ,OceanGrid%sinLat(ilat) 
		end do

		do iphi=1,OceanGrid%nPhi
			phi=half*OceanGrid%dphi+OceanGrid%dphi*real(iphi-1,QR_K)
			OceanGrid%sinPhiC(iphi)=sin(phi)
			OceanGrid%cosPhiC(iphi)=cos(phi)
		end do

		!COMMENT: If Bathymetry file is not read, then appending depth of our choice
		if(.not.include_bathymetry)then
			if(allocated(OceanGrid%d))deallocate(OceanGrid%d)
			if(.NOT.allocated(OceanGrid%d))allocate(OceanGrid%d(OceanGrid%nLat,OceanGrid%nPhi))
			!ALTERED: d0Constant value is set in the data structure. !3423.56 !((two*OMEGA*Ld0)**2)/gravity
			OceanGrid%d0 = d0Constant ; 
			do ilat=1,OceanGrid%nLat
				do iphi=1,OceanGrid%nPhi
					OceanGrid%d(ilat,iPhi) = OceanGrid%d0
					! write(*,*) 'OceanGrid%d(ilat,iPhi) = ', OceanGrid%d(ilat,iPhi)  
				end do
			end do
		write(*,*)'Bathymetry is not included and the depth is set to d₀ =', OceanGrid%d0
	end if
	write(*,*)'Ocean Grid is Constructed ***' 
	end subroutine constructOceanGrid

	!TOFIX: Subroutine (in progress — not yet ready for use or aligned with current workflow). to Construct the AdaptiveOceanGrid 
	subroutine interpolateBathymetryToAdaptiveGrid()
		use OceanAdaptiveGridDataModule
		use basicDataStructures
		use simulationParameters
		use bathymetryDataModule
		!$ use omp_lib
		implicit none
		integer :: CHUNK, i1, i2, threadID, nThreads
		real (kind=QR_K) :: Lat0,phi0,LatB,LatT,phiL,phiR
		real (kind=QR_K) :: alpha,zphiL,zphiR,lat,sinLat,beta,aalpha
		integer :: i,iLatB,iLatT,iphiL,iphiR,iLat,iPhi,nnn
		real (kind=QR_K) :: dLatB, dPhiB, zT, zB, aveDepth
		! Mapping arrays for ocean points only
		integer, allocatable :: oceanPointIndex(:,:)  ! Maps (iLat,iLon) to ocean point index
		integer :: nOceanPoints,iLon

		if(allocated(OceanAdaptiveGrid%d))deallocate(OceanAdaptiveGrid%d)
		if(.NOT.allocated(OceanAdaptiveGrid%d))allocate(OceanAdaptiveGrid%d(OceanAdaptiveGrid%nLat,OceanAdaptiveGrid%nPhi))
		OceanAdaptiveGrid%d=zero

		dLatB = pi/real(bathymetry%nLat,QR_K)
		dPhiB = two*pi/real(bathymetry%nLong,QR_K)

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = OceanAdaptiveGrid%nLat/nThreads
		write(*,*) 'Chunks', CHUNK
		
		i1=1; i2=OceanAdaptiveGrid%nLat

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,nThreads,bathymetry,OceanAdaptiveGrid,&
		!$OMP &dLatB,dPhiB)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		if(threadID+1.eq.nThreads)i2=OceanAdaptiveGrid%nLat
		do iLat=i1,i2
			Lat0=-half*pi+half*OceanAdaptiveGrid%dLat+OceanAdaptiveGrid%dLat*real(ilat-1,QR_K)
			iLatB = 1+int((Lat0+half*pi-half*dLatB)/dLatB)
			iLatT = iLatB+1
			LatB=-half*pi+half*dLatB+dLatB*real(iLatB-1,QR_K)
			LatT=-half*pi+half*dLatB+dLatB*real(iLatT-1,QR_K)

			do iphi = 1,OceanAdaptiveGrid%nPhi
				phi0=half*OceanAdaptiveGrid%dphi+OceanAdaptiveGrid%dphi*real(iphi-1,QR_K)

				iphiL = 1+int((phi0-half*dphiB)/dphiB)
				iphiR = iphiL+1; if(iphiR.gt.bathymetry%nLong)iphiR=1
				phiL= half*dphiB+dphiB*real(iphiL-1,QR_K)
				phiR= half*dphiB+dphiB*real(iphiR-1,QR_K)

				alpha = (Lat0-LatB)/(LatT-LatB)
				zT = real(bathymetry%z(iphiL,iLatT),QR_K)
				zB = real(bathymetry%z(iphiL,iLatB),QR_K)
				zphiL = zB + alpha*(zT-zB)
				zT = real(bathymetry% z(iphiR,iLatT),QR_K)
				zB = real(bathymetry%z(iphiR,iLatB),QR_K)
				zphiR = zB + alpha*(zT-zB)
				alpha = (phi0-phiL)/(phiR-phiL)
				OceanAdaptiveGrid%d(iLat,iPhi) = -(zphiL + alpha*(zphiR-zphiL)) !!all points ocean and land
			end do
		end do
		!$OMP END PARALLEL
		nOceanPoints = 0
		do iLat=1,OceanAdaptiveGrid%nLat
			do iLon=1,OceanAdaptiveGrid%nPhi
				if (OceanAdaptiveGrid%d(iLat,iLon) > zero) then
					nOceanPoints = nOceanPoints + 1
				endif
			end do
		end do
		
		write(*,*)'Number of only ocean points: ', nOceanPoints
		write(*,*)'OceanAdaptiveGrid%d land and ocean=', size(OceanAdaptiveGrid%d)
		! write(*,*)'OceanAdaptiveGrid%d land and ocean=',OceanAdaptiveGrid%d	
		
		nnn=0
		aveDepth=zero
		do iLat=1,OceanAdaptiveGrid%nLat
			do iphi = 1,OceanAdaptiveGrid%nPhi
				if(OceanAdaptiveGrid%d(iLat,iPhi).gt.zero) then
					aveDepth=(real(nnn,QR_K)*aveDepth+OceanAdaptiveGrid%d(iLat,iPhi))/real(nnn+1,QR_K)
					nnn=nnn+1
				end if
			end do
		end do

		OceanAdaptiveGrid%d0 = aveDepth
		Ld0 = sqrt(GRAVITY*OceanAdaptiveGrid%d0)/(two*OMEGA)

		write(*,*)'*** Mean depth is computed, d0 = ***',OceanAdaptiveGrid%d0

		write(*,*)'*** The Rossby Radius of Deformation is computed @ Ld0 = sqrt(GRAVITY*OceanAdaptiveGrid%d0)/(two*OMEGA) =  ***',Ld0

		return
	end subroutine interpolateBathymetryToAdaptiveGrid

	!TOFIX: Subroutine (in progress — not yet ready for use or aligned with current workflow).
	subroutine constructAdaptiveGrid()
		use OceanAdaptiveGridDataModule
		use simulationParameters
		implicit none
		integer :: iphi,ilat
		real (kind=QR_K) :: lat, phi, colatitude

		OceanAdaptiveGrid%nPhi=2*OceanAdaptiveGrid%nLat !!Should be twice n theta 

		OceanAdaptiveGrid%dLat=pi/real(OceanAdaptiveGrid%nLat,QR_K)
		OceanAdaptiveGrid%dphi=(two*pi)/real(OceanAdaptiveGrid%nPhi,QR_K)

		write(*,*) '*** CONSTRUCT GRID ADAPTIVE GRID ***'
		
		write(*,*)' '
		write(*,*) 'nTheta, nPhi = ', OceanAdaptiveGrid%nLat, OceanAdaptiveGrid%nPhi
		
		write(*,*) '(dLat, dPhi) = (', OceanAdaptiveGrid%dLat, ' rad,', OceanAdaptiveGrid%dPhi , ' rad )'
		write(*,*)  '(dLat, dPhi) = (',OceanAdaptiveGrid%dLat * 180.0/pi, ' deg, ',OceanAdaptiveGrid%dPhi * 180.0/pi, ' deg )'
		write(*,*)  '(dLat, dPhi) = (',OceanAdaptiveGrid%dLat* R*0.0001, ' km, ', OceanAdaptiveGrid%dPhi*0.0001* R, ' km)'
	 	write(*,*)' '

		! 'dLat = ', OceanAdaptiveGrid%dLat , 'dphi = ', OceanAdaptiveGrid%dphi
		! write(*,*)'ntheta , nPhi=',OceanAdaptiveGrid%nLat, OceanAdaptiveGrid%nPhi
		!! Allocate arrays for sine and cosine of latitudes and longitudes
		if(.NOT.allocated(OceanAdaptiveGrid%sinLat))allocate(OceanAdaptiveGrid%sinLat(OceanAdaptiveGrid%nLat+1))
		if(.NOT.allocated(OceanAdaptiveGrid%cosLat))allocate(OceanAdaptiveGrid%cosLat(OceanAdaptiveGrid%nLat+1))
		if(.NOT.allocated(OceanAdaptiveGrid%sinLatC))allocate(OceanAdaptiveGrid%sinLatC(OceanAdaptiveGrid%nLat))
		if(.NOT.allocated(OceanAdaptiveGrid%cosLatC))allocate(OceanAdaptiveGrid%cosLatC(OceanAdaptiveGrid%nLat))

		if(.NOT.allocated(OceanAdaptiveGrid%sinPhi))allocate(OceanAdaptiveGrid%sinPhi(OceanAdaptiveGrid%nPhi))
		if(.NOT.allocated(OceanAdaptiveGrid%cosPhi))allocate(OceanAdaptiveGrid%cosPhi(OceanAdaptiveGrid%nPhi))
		if(.NOT.allocated(OceanAdaptiveGrid%sinPhiC))allocate(OceanAdaptiveGrid%sinPhiC(OceanAdaptiveGrid%nPhi))
		if(.NOT.allocated(OceanAdaptiveGrid%cosPhiC))allocate(OceanAdaptiveGrid%cosPhiC(OceanAdaptiveGrid%nPhi))

		!!cell vertices: -pi/2 -> pi/2
		do ilat=1,OceanAdaptiveGrid%nLat+1
			lat=-half*pi+OceanAdaptiveGrid%dLat*real(ilat-1,QR_K)
			OceanAdaptiveGrid%sinLat(ilat)=sin(lat)
			OceanAdaptiveGrid%cosLat(ilat)=cos(lat)
		end do

		do iphi=1,OceanAdaptiveGrid%nPhi
			phi=OceanAdaptiveGrid%dphi*real(iphi-1,QR_K)
			OceanAdaptiveGrid%sinPhi(iphi)=sin(phi)
			OceanAdaptiveGrid%cosPhi(iphi)=cos(phi)
		end do

		!!cell centers
		do ilat=1,OceanAdaptiveGrid%nLat
			lat=-half*pi+half*OceanAdaptiveGrid%dLat+OceanAdaptiveGrid%dLat*real(ilat-1,QR_K)
			OceanAdaptiveGrid%sinLatC(ilat)=sin(lat)
			OceanAdaptiveGrid%cosLatC(ilat)=cos(lat)
		end do

		do iphi=1,OceanAdaptiveGrid%nPhi
			phi=half*OceanAdaptiveGrid%dphi+OceanAdaptiveGrid%dphi*real(iphi-1,QR_K)
			OceanAdaptiveGrid%sinPhiC(iphi)=sin(phi)
			OceanAdaptiveGrid%cosPhiC(iphi)=cos(phi)
		end do

		if(allocated(OceanAdaptiveGrid%dA))deallocate(OceanAdaptiveGrid%dA)
		if(.NOT.allocated(OceanAdaptiveGrid%dA))allocate(OceanAdaptiveGrid%dA(OceanAdaptiveGrid%nLat))

		do iLat = 1, OceanAdaptiveGrid%nLat
			Lat = -half*pi+half*OceanAdaptiveGrid%dLat+OceanAdaptiveGrid%dLat*real(ilat-1,QR_K)
			colatitude = half*pi-Lat;
			OceanAdaptiveGrid%dA(iLat) =  (cos(colatitude-half*OceanAdaptiveGrid%dLat)-cos(colatitude + half*OceanAdaptiveGrid%dLat))*OceanAdaptiveGrid%dphi
			! write(*,*) 'dA = ', OceanAdaptiveGrid%dA(iLat)
		end do

		!if we do not have bathymetry then we assign an averaged depth on our own. from bathymetry the averaged depth is d = 
		if(.not.include_bathymetry)then
			if(allocated(OceanAdaptiveGrid%d))deallocate(OceanAdaptiveGrid%d)
			if(.NOT.allocated(OceanAdaptiveGrid%d))allocate(OceanAdaptiveGrid%d(OceanAdaptiveGrid%nLat,OceanAdaptiveGrid%nPhi))
			OceanAdaptiveGrid%d0 = 4000 !((two*OMEGA*Ld0)**2)/gravity
			do ilat=1,OceanAdaptiveGrid%nLat
				do iphi=1,OceanAdaptiveGrid%nPhi
					OceanAdaptiveGrid%d(ilat,iPhi) = OceanAdaptiveGrid%d0
				end do
			end do
			write(*,*)'Bathymetry is not included, OceanAdaptiveGrid%d0=', OceanAdaptiveGrid%d0
		end if

		write(*,*)'*** Ocean: Adaptive Grid constructed ***' 
		return
	end subroutine constructAdaptiveGrid

end module gridOperationsModule

