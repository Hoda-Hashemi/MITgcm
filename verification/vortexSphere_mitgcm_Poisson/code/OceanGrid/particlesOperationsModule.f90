module particlesOperationsModule
	implicit none
	public :: moveParticles


contains

    !COMMENT: Inside this subroutine Boundary conditions were revised and added.
	subroutine moveParticles(fraction)
		use particlesDataModule
		use simulationParameters
		use basicDataStructures
		implicit none
		integer :: i
		real (kind=QR_K) :: fraction
		do i=1,nParticles
			particles(i)%el(1)%lat = particles(i)%el(1)%lat + fraction*(deltat/R)*particles(i)%el(1)%uLat
			particles(i)%el(1)%Phi = particles(i)%el(1)%Phi + fraction*(deltat/R)*particles(i)%el(1)%uLon

			! if(particles(i)%el(1)%Phi.gt.2.0*pi)particles(i)%el(1)%Phi=particles(i)%el(1)%Phi-2.0*pi
			! if(particles(i)%el(1)%Phi.lt.0.0)particles(i)%el(1)%Phi=particles(i)%el(1)%Phi+2.0*pi
			
			!COMMENT: Added Neumann Boundary conditions and Periodic BC
			if (particles(i)%el(1)%lat > half*pi) then
				particles(i)%el(1)%lat = pi - particles(i)%el(1)%lat
				particles(i)%el(1)%Phi = particles(i)%el(1)%Phi + pi
			else if (particles(i)%el(1)%lat < -half*pi) then
				particles(i)%el(1)%lat = -pi - particles(i)%el(1)%lat
				particles(i)%el(1)%Phi = particles(i)%el(1)%Phi + pi
			end if
			if (particles(i)%el(1)%Phi > 2.0*pi) then
				particles(i)%el(1)%Phi = particles(i)%el(1)%Phi - 2.0*pi
			else if (particles(i)%el(1)%Phi < 0.0) then
				particles(i)%el(1)%Phi = particles(i)%el(1)%Phi + 2.0*pi
			end if 
			
			particles(i)%el(1)%sinLat=sin(particles(i)%el(1)%lat)
			particles(i)%el(1)%cosLat=cos(particles(i)%el(1)%lat)
			particles(i)%el(1)%sinPhi=sin(particles(i)%el(1)%Phi)
			particles(i)%el(1)%cosPhi=cos(particles(i)%el(1)%Phi)
		end do

		return
	end subroutine moveParticles




end module particlesOperationsModule
