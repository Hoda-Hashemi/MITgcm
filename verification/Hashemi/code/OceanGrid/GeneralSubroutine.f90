module GeneralSubroutine
use CoreFunctions
implicit none

public :: ReorderNew

contains
!############################################################
subroutine ReorderNew()
	use QuickSearchVariables
	implicit none

	integer :: i, j
	integer :: n, m, l, k
	integer :: iPhi, iLat, iz, izone
	Double precision :: phiC, latC, zc
	integer, allocatable, dimension(:) :: counter

! step 1 find phiMin, phiMax, phiMin, latMax

    phiMin = MINVAL(particles(1:nParticles)%Lon)
    phiMax = MAXVAL(particles(1:nParticles)%Lon)
    latMin = MINVAL(particles(1:nParticles)%lat)
    latMax = MAXVAL(particles(1:nParticles)%lat)

	phiC = 0.50D0*(phiMax+phiMin)
	latC = 0.50D0*(latMax+latMin)

! calculate the number of zones based on max, min and zone size
	nPhiZones = INT((phiMax-phiMin)/zoneSize+Epsilons) + 3
	nLatZones = INT((latMax-latMin)/zoneSize+Epsilons) + 3
	!nBlck = nyZones*nzZones
	!nBlck =nLatZones
	nZones=nPhiZones*nLatZones

! All the frontiere zones does not contain any element

   phiMin = phiC - 0.5D0*DBLE(nPhiZones)*zoneSize
   latMin = latC - 0.5D0*DBLE(nLatZones)*zoneSize

   phiMax = phiC + 0.5D0*DBLE(nPhiZones)*zoneSize
   latMax = latC + 0.5D0*DBLE(nLatZones)*zoneSize

	if(allocated(zones))deallocate(zones)
	allocate(zones(nZones))
	allocate(counter(nZones))

! initialize number of elements of zone i to zero
	zones(1:nZones)%nEls=0
	counter(1:nZones)=0

! The zones are indexed from 1 to nZones marching first in the y direction then z and then x
! march over the elements to determine:
! (a) number of elements in each zone
! (b) the zone to which element i belongs to
	do i=1,nParticles
		iPhi=int((particles(i)%Lon-phiMin)/zoneSize)+1
		iLat=int((particles(i)%lat-latMin)/zoneSize)+1
		!izone=(ix-1)*nBlck + (iz-1)*nyZones + iy
	    izone=(iPhi-1)*nLatZones + iLat
	    zones(izone)%nEls=zones(izone)%nEls+1
	    particles(i)%indexZone=izone
	end do

! store the index of the fist element in each zone in startLoc
	zones(1)%startLoc=1
	do i=2,nZones
		zones(i)%startLoc=zones(i-1)%startLoc+zones(i-1)%nEls
	end do

! reord all elements in zones
	do i=1,nParticles
	  izone=particles(i)%indexZone
	  counter(izone)=counter(izone)+1
	  j=zones(izone)%startLoc+counter(izone)-1
	  particles(j)%index = i
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
!############################################################

end module GeneralSubroutine

