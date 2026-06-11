module outputModule

	implicit none
	!COMMENT: Original Legacy Code
	public :: plotXY
	public :: saveParticlesVTU
	public :: saveMesh
	public :: saveBathymetry

	!COMMENT: Added subroutines.
	public :: saveOceanAdaptiveMesh
	public :: saveOceanMesh
	public :: saveOceanMesh2D
	public :: saveBathymetry2D
	public :: windVelocity

contains
	subroutine plotXY(fileName,xLabel,yLabel)
		use simulationParameters
		implicit none
		CHARACTER(len=256) :: xlabel,ylabel,fileName,command
		write(*,*)trim(adjustl(osName))
		if(trim(adjustl(osName)).eq.'Darwin')then
			call system('cp ../resources/plotXY0MAC.plt ../resources/plotXY.plt')
			command="sed -i '' "//"'s/"//"SEDxLabel/"//trim(adjustl(xlabel))//"/g' ../resources/plotXY.plt"
			call system(trim(adjustl(command)))
			command="sed -i '' "//"'s/"//"SEDyLabel/"//trim(adjustl(ylabel))//"/g' ../resources/plotXY.plt"
			call system(trim(adjustl(command)))
			command="sed -i '' "//"'s/"//"SEDfileName/"//trim(adjustl(fileNAme))//"/g' ../resources/plotXY.plt"
			call system(trim(adjustl(command)))
			call system('gnuplot ../resources/plotXY.plt')
		end if

		if(trim(adjustl(osName)).eq.'Linux')then
			call system('cp ../resources/plotXY0.plt ../resources/plotXY.plt')
			command="sed -i "//"'s/"//"SEDxLabel/"//trim(adjustl(xlabel))//"/g' ../resources/plotXY.plt"
			call system(trim(adjustl(command)))
			command="sed -i "//"'s/"//"SEDyLabel/"//trim(adjustl(ylabel))//"/g' ../resources/plotXY.plt"
			call system(trim(adjustl(command)))
			command="sed -i "//"'s/"//"SEDfileName/"//trim(adjustl(fileNAme))//"/g' ../resources/plotXY.plt"
			call system(trim(adjustl(command)))
			call system('gnuplot ../resources/plotXY.plt')
		end if

	end subroutine plotXY

	subroutine saveParticlesVTU(kkk)
		use basicDataStructures
		use simulationParameters
		use sphericalGridDataModule
		use particlesDataModule
		implicit none
		integer, intent(in) :: kkk
		CHARACTER  :: STRING*80
		CHARACTER(len=32) :: xrstring
		CHARACTER(len=32) :: yrstring
		CHARACTER(len=32) :: zrstring
		CHARACTER(len=10) :: xstring
		!	CHARACTER(len=10) :: ystring
		!	CHARACTER(len=10) :: zstring
		integer :: i !,stat
		real (kind=QR_K) :: cosTheta,sinTheta,cosPhi,sinPhi,Vr,Vtheta,Vphi,ux,uy,uz

		write(*,*)'saveParticlesVTU ',kkk
		WRITE(STRING, FMT='(I5)') kkk

		open(220,file='../results/particles'//trim(adjustl(STRING))//'.vtu')

		write(220,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
		write(220,*)'  <UnstructuredGrid>'

		write(xstring,'(i10)' )nParticles
		write(220,*)'    <Piece NumberOfPoints="'//trim(adjustl(xstring))//'" NumberOfCells="'//trim(adjustl(xstring))//'">'

		write(220,*)'      <PointData Scalars="scalars">'
		write(220,*)'        <DataArray type="Float32" Name="strength" Format="ascii">'
		do i=1,nParticles
		!			!write(220,*)particles(i)%el(1)%Q
			write(xrstring,'(f16.10)' )particles(i)%el(1)%strengthVor/(sphericalGrid%dLat**2)
			write(220,*)xrstring !DBLE(particles(i)%el(1)%strengthVor/(sphericalGrid%dLat**2))
		end do
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="Float32" Name="h" Format="ascii">'
		do i=1,nParticles
		!			write(xrstring,'(f16.10)' )particles(i)%el(1)%h
			write(220,*)real(particles(i)%el(1)%h) !xrstring
		end do
		write(220,*)'        </DataArray>'
		write(220,*)'        <DataArray type="Float32" Name="hs" Format="ascii">'
		do i=1,nParticles
			!write(xrstring,'(f16.10)' )particles(i)%el(1)%h-particles(i)%el(1)%d
			write(220,*)real(particles(i)%el(1)%h-particles(i)%el(1)%d)
		end do
		write(220,*)'        </DataArray>'
		write(220,*)'        <DataArray type="Float32" Name="d" Format="ascii">'
		do i=1,nParticles
			!write(xrstring,'(f16.10)' )particles(i)%el(1)%d
			write(220,*)real(particles(i)%el(1)%d)
		end do
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="Float32" Name="div" Format="ascii">'
		do i=1,nParticles
			!write(xrstring,'(f16.10)' )particles(i)%el(1)%strengthDiv
			write(220,*)real(particles(i)%el(1)%strengthDiv)
		end do
		write(220,*)'        </DataArray>'

		!		write(220,*)'        <DataArray type="Float32" Name="Qmf" Format="ascii">'
		!		do i=1,nParticles
		!			write(xrstring,'(f16.10)' )particles(i)%el(1)%Q-two*OMEGA*sin(particles(i)%el(1)%Lat)
		!			write(220,*)xrstring
		!			!			write(*,*)particles(i)%el(1)%Q-two*OMEGA*sin(particles(i)%el(1)%Lat)
		!		end do
		!		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="Float32" NumberOfComponents="3" Name="velocity" Format="ascii">'
		do i=1,nParticles
			! Vx = sin(theta)*cos(phi) Vr + cos(theta) cos(phi) Vtheta - sin(phi) Vphi
			! Vy = sin(theta)*sin(phi) Vr + cos(theta) sin(phi) Vtheta + cos(phi) Vphi
			! Vz = cos(theta)          Vr - sin(theta)          Vtheta + 0        Vphi
			cosTheta=particles(i)%el(1)%sinLat
			sinTheta=particles(i)%el(1)%cosLat
			cosPhi=particles(i)%el(1)%cosPhi
			sinPhi=particles(i)%el(1)%sinPhi
			Vr=zero;

			!if(kkk.eq.9999)then
			!	Vtheta=particles(i)%el(1)%uLatExact; Vphi=particles(i)%el(1)%uLonExact;
			!else
				Vtheta=particles(i)%el(1)%uLat; Vphi=particles(i)%el(1)%uLon;
			!end if
			ux=(sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi)
			uy=(sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi)
			uz=(cosTheta*Vr - sinTheta*Vtheta)

		if(isnan(ux).or.isnan(uy).or.isnan(uz))then
		write(*,*)'in save particles',Vtheta,Vphi,cosTheta,sinPhi,cosPhi
		stop
		end if

			write(xrstring,'(f16.10)' )ux
			write(yrstring,'(f16.10)' )uy
			write(zrstring,'(f16.10)' )uz
			!write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
			write(220,*)real(ux)
			write(220,*)real(uy)
			write(220,*)real(uz)
		end do
		write(220,*)'        </DataArray>'
		write(220,*)'      </PointData>'

		write(220,*)'      <Points>'
		write(220,*)'        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
		do i=1,nParticles
			! x = R*sin(colat)*cos(phi)=R*cos(lat)*cos(phi)
			! y = R*sin(colat)*sin(phi)=R*cos(lat)*sin(phi)
			! z = R*cos(colat) = R*sin(lat)
			write(xrstring,'(f16.10)' )particles(i)%el(1)%cosLat*cos(particles(i)%el(1)%Phi) !*R
			write(yrstring,'(f16.10)' )particles(i)%el(1)%cosLat*sin(particles(i)%el(1)%Phi) !*R
			write(zrstring,'(f16.10)' )particles(i)%el(1)%sinLat !*R
			write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
		end do
		write(220,*)'        </DataArray>'
		write(220,*)'      </Points>'

		write(220,*)'     <Cells>'
		write(220,*)'        <DataArray type="Int32" Name="connectivity" Format="ascii">'
		do i=1,nParticles
			write(xstring,'(i10)' )i-1
			write(220,*)trim(adjustl(xstring))
		end do
		write(220,*)'	        </DataArray>'
		write(220,*)'        <DataArray type="Int32" Name="offsets" Format="ascii">'
		do i=1,nParticles
			write(xstring,'(i10)' )i
			write(220,*)trim(adjustl(xstring))
		end do
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="UInt8" Name="types" Format="ascii">'
		do i=1,nParticles
			write(xstring,'(i10)' )1
			write(220,*)trim(adjustl(xstring))
		end do
		write(220,*)'        </DataArray>'
		write(220,*)'      </Cells>'

		write(220,*)'    </Piece>'

		write(220,*)'  </UnstructuredGrid>'

		write(220,*)'</VTKFile>'
		close(220)

	end subroutine saveParticlesVTU

	subroutine saveBathymetry()
		use OceanGridDataModule
		use basicDataStructures
		use simulationParameters
		use bathymetryDataModule
		
		implicit none
		!CHARACTER  :: STRING*80
		CHARACTER(len=32) :: xrstring
		CHARACTER(len=32) :: yrstring
		CHARACTER(len=32) :: zrstring
		CHARACTER(len=10) :: xstring
		CHARACTER(len=10) :: ystring
		!	CHARACTER(len=10) :: zstring
		integer :: i,indxLL,indxLR,indxUR,indxUL, iLat, iphi,iphiP1,Rb, HodaRb,LdPM2,f,d,meow
		real (kind=QR_K) :: dtheta
		real (kind=QR_K) :: phi, sinLat, cosLat
		real (kind=QR_K) :: sinPhi, cosPhi, lat, zzz, scale !, Vr, Vtheta, Vphi

		scale = real(30,QR_K)
		dtheta = pi/(bathymetry%nLat-1)

		open(220,file='../results/bathymetry.vtu')

		write(220,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
		write(220,*)'  <UnstructuredGrid>'

		write(xstring,'(i10)' )(bathymetry%nLat-1)*bathymetry%nLong
		write(ystring,'(i10)' )bathymetry%nLat*bathymetry%nLong
		write(220,*)'    <Piece NumberOfPoints="'//trim(adjustl(ystring))//'" NumberOfCells="'//trim(adjustl(xstring))//'">'

		 write(220,*)'      <PointData Scalars="scalars">'

		 write(220,*)'        <DataArray type="Int32" Name="Z" Format="ascii">'
		 do iLat=1,bathymetry%nLat
		 	do iphi=1,bathymetry%nLong
		 		write(xrstring,'(I7)' )bathymetry%z(iPhi,iLat)
		 		write(220,*)xrstring
		 	end do
		 end do
		 write(220,*)'        </DataArray>'

		 write(220,*)'      </PointData>'

		write(220,*)'      <Points>'
		write(220,*)'        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

		do iLat=1,bathymetry%nLat
			lat=-half*pi+dtheta*real(iLat-1,QR_K)
			cosLat=cos(lat)
			sinLat=sin(lat)
			do iphi=1,bathymetry%nLong
				phi=dtheta*real(iphi-1,QR_K)
				sinPhi=sin(phi)
				cosPhi=cos(phi)
				zzz=real(bathymetry%z(iPhi,iLat),QR_K)
				write(xrstring,'(f16.10)' )cosLat*cosPhi*(one+scale*zzz/R)
				write(yrstring,'(f16.10)' )cosLat*sinPhi*(one+scale*zzz/R)
				write(zrstring,'(f16.10)' )sinLat*(one+scale*zzz/R)
				write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
			end do
		end do

		write(220,*)'        </DataArray>'
		write(220,*)'      </Points>'

		write(220,*)'     <Cells>'
		write(220,*)'        <DataArray type="Int32" Name="connectivity" Format="ascii">'

		do iLat=1,bathymetry%nLat-1
			do iphi=1,bathymetry%nLong
				iphiP1=iphi+1
				if(iphi.eq.bathymetry%nLong)iphiP1=1
				indxLL=(iLat-1)*bathymetry%nLong+iphi
				indxLR=(iLat-1)*bathymetry%nLong+iphiP1
				indxUR=iLat*bathymetry%nLong+iphiP1
				indxUL=iLat*bathymetry%nLong+iphi
				write(220,*)indxLL-1,' ',indxLR-1,' ',indxUL-1,' ',indxUR-1
			end do
		end do
		write(220,*)'	        </DataArray>'

		write(220,*)'        <DataArray type="Int32" Name="offsets" Format="ascii">'
		do i=1,bathymetry%nLong*(bathymetry%nLat-1)
			write(xstring,'(i10)' )i*4
			write(220,*)trim(adjustl(xstring))
		end do
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="UInt8" Name="types" Format="ascii">'
		do i=1,bathymetry%nLong*(bathymetry%nLat-1)
			write(xstring,'(i10)' )8
			write(220,*)trim(adjustl(xstring))
		end do
		write(220,*)'        </DataArray>'
		write(220,*)'      </Cells>'

		write(220,*)'    </Piece>'

		write(220,*)'  </UnstructuredGrid>'

		write(220,*)'</VTKFile>'
		close(220)
		
	end subroutine saveBathymetry

	!COMMENT: this is from the old code. commented out some quantities. maybe added new quantities too.
	subroutine saveMesh(kkk)
		use basicDataStructures
		use simulationParameters
		use sphericalGridDataModule
		use particlesDataModule
		implicit none 
		integer, intent(in) :: kkk
		CHARACTER  :: STRING*80
		CHARACTER(len=32) :: xrstring
		CHARACTER(len=32) :: yrstring
		CHARACTER(len=32) :: zrstring
		CHARACTER(len=10) :: xstring
		CHARACTER(len=10) :: ystring
		!	CHARACTER(len=10) :: zstring
		integer :: i,indxLL,indxLR,indxUR,indxUL, iLat, iphi,iphiP1
		real (kind=QR_K) , allocatable, dimension(:,:) :: x, y, z
		real (kind=QR_K) :: phi, sinTheta, cosTheta
		real (kind=QR_K) :: sinPhi, cosPhi, lat, Vr, Vtheta, Vphi

		! sphericalGrid%nLat = 6
		! sphericalGrid%nPhi = 12
		! sphericalGrid%dLat=pi/real(sphericalGrid%nLat,QR_K)
		! sphericalGrid%dphi=two*pi/real(sphericalGrid%nPhi,QR_K)
		
		if(.NOT.allocated(velocityOnGrid))then
			allocate(velocityOnGrid(2,sphericalGrid%nLat,sphericalGrid%nPhi))
			velocityOnGrid=zero
		end if

		! if(.NOT.allocated(phiOnGrid))then
		! 	allocate(phiOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		! 	phiOnGrid=zero
		! end if
		if(.NOT.allocated(psiOnGrid))then
			allocate(psiOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
			psiOnGrid=zero
		end if
		! if(.NOT.allocated(QOnGrid))then
		! 	allocate(QOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		! 	QOnGrid=zero
		! end if
		!if(.NOT.allocated(velocityOnGridExact))then
		!	allocate(velocityOnGridExact(2,sphericalGrid%nLat,sphericalGrid%nPhi))
		!	velocityOnGridExact=zero
		!end if
		!if(.NOT.allocated(psiOnGridExact))then
		!	allocate(psiOnGridExact(sphericalGrid%nLat,sphericalGrid%nPhi))
		!	psiOnGridExact=zero
		!end if

		allocate(x(sphericalGrid%nLat+1,sphericalGrid%nPhi))
		allocate(y(sphericalGrid%nLat+1,sphericalGrid%nPhi))
		allocate(z(sphericalGrid%nLat+1,sphericalGrid%nPhi))

		do iLat=1,sphericalGrid%nLat+1
			lat=-half*pi+sphericalGrid%dLat*real(iLat-1,QR_K)
			do iphi=1,sphericalGrid%nPhi
				phi=sphericalGrid%dphi*real(iphi-1,QR_K)
				x(iLat,iphi)=sphericalGrid%cosLat(ilat)*sphericalGrid%cosPhi(iphi) !*R
				y(iLat,iphi)=sphericalGrid%cosLat(ilat)*sphericalGrid%sinPhi(iphi) !*R
				z(iLat,iphi)=sphericalGrid%sinLat(ilat) !*R
			end do
		end do

		WRITE(STRING, FMT='(I5)') kkk
		open(220,file='../results/mesh'//trim(adjustl(STRING))//'.vtu')

		write(220,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
		write(220,*)'  <UnstructuredGrid>'

		write(xstring,'(i10)' )sphericalGrid%nLat*sphericalGrid%nPhi
		write(ystring,'(i10)' )(sphericalGrid%nLat+1)*sphericalGrid%nPhi
		write(220,*)'    <Piece NumberOfPoints="'//trim(adjustl(ystring))//'" NumberOfCells="'//trim(adjustl(xstring))//'">'

		! write(220,*)'      <PointData Scalars="scalars">'
		! write(220,*)'        <DataArray type="Float32" Name="Psi/R Orig" Format="ascii">'
		! m=4; a = real(4.1D7,QR_K)
		! do iLat=1,sphericalGrid%nLat+1
		! 	lat=-half*pi+sphericalGrid%dLat*real(iLat-1,QR_K)
		! 	cosLatPowerm=sphericalGrid%cosLat(ilat)**m
		! 	do iphi=1,sphericalGrid%nPhi
		! 		phi=sphericalGrid%dphi*real(iphi-1,QR_K)
		! 		psiOverR=-(a/R)*sphericalGrid%sinLat(ilat)*cosLatPowerm*cos(real(m,QR_K)*phi)
		! 		write(xrstring,'(f16.10)' )psiOverR
		! 		write(220,*)xrstring
		! 	end do
		! end do
		! write(220,*)'        </DataArray>'
		! write(220,*)'        <DataArray type="Float32" Name="Psi/R Computed" Format="ascii">'
		! do iLat=1,sphericalGrid%nLat+1
		! 	do iphi=1,sphericalGrid%nPhi
		! 		write(xrstring,'(f16.10)' )psiOnGrid(iLat,iphi)/R
		! 		write(220,*)xrstring
		! 	end do
		! end do
		! write(220,*)'        </DataArray>'
		! write(220,*)'      </PointData>'

		write(220,*)'      <CellData Scalars="scalars">'

		write(220,*)'        <DataArray type="Float32" Name="Psi Poisson" Format="ascii">'
		do iLat=1,sphericalGrid%nLat
			do iphi=1,sphericalGrid%nPhi
				! write(xrstring,'(f16.10)' ) psiOnGrid(iLat,iphi)
				write(xrstring,'(E20.10)') psiOnGrid(iLat,iphi)
				write(220,*)xrstring
			end do
		end do
		write(220,*)'        </DataArray>'

		! write(220,*)'        <DataArray type="Float32" Name="Phi Poisson" Format="ascii">'
		! do iLat=1,sphericalGrid%nLat
		! 	do iphi=1,sphericalGrid%nPhi
		! 		write(xrstring,'(f16.10)' )phiOnGrid(iLat,iphi)

		! 		write(220,*)phiOnGrid(iLat,iphi) !xrstring
		! 	end do
		! end do
		! write(220,*)'        </DataArray>'

		!if(1.eq.2)then
		!write(220,*)'        <DataArray type="Float32" Name="Psi/R Exact" Format="ascii">'
		!do iLat=1,sphericalGrid%nLat
		!	do iphi=1,sphericalGrid%nPhi
		!		write(xrstring,'(f16.10)' )psiOnGridExact(iLat,iphi)/R
		!		write(220,*)xrstring
		!	end do
		!end do
		!write(220,*)'        </DataArray>'
		!end if

		! write(220,*)'        <DataArray type="Float32" Name="Q" Format="ascii">'
		! do iLat=1,sphericalGrid%nLat
		! lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
		! do iphi=1,sphericalGrid%nPhi
		! 	phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi-1,QR_K)
		! 	write(xrstring,'(f16.10)' )QOnGrid(iLat,iPhi)
		! 	write(220,*)QOnGrid(iLat,iPhi) !xrstring
		! 	end do
		! end do
		! write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="Float32" Name="depth" Format="ascii">'
		do iLat=1,sphericalGrid%nLat
			do iphi=1,sphericalGrid%nPhi
				!write(xrstring,'(f16.10)' )sphericalGrid%d(iLat,iphi)
				write(220,*)sphericalGrid%d(iLat,iphi)
				end do
			end do
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="Float32" Name=" - f" Format="ascii">'
		do iLat=1,sphericalGrid%nLat
			lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
			do iphi=1,sphericalGrid%nPhi
				phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi-1,QR_K)
				write(xrstring,'(f16.10)' )-two*OMEGA*sphericalGrid%sinLatC(ilat)
				write(220,*)xrstring
			end do
		end do
		write(220,*)'        </DataArray>'

		! if(1.eq.2)then
		! 		write(220,*)'        <DataArray type="Float32" NumberOfComponents="3" Name="velocity" Format="ascii">'
		! 		do iLat=1,sphericalGrid%nLat
		! 			lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
		! 			sinTheta=sphericalGrid%cosLatC(ilat)
		! 			cosTheta=sphericalGrid%sinLatC(ilat)
		! 			do iphi=1,sphericalGrid%nPhi
		! 				phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi-1,QR_K)
		! 				sinPhi=sphericalGrid%sinPhiC(iphi)
		! 				cosPhi=sphericalGrid%cosPhiC(iphi)
		! !				if(kkk.eq.9999)then
		! !					Vr=zero;Vtheta=velocityOnGridExact(1,iLat,iphi); Vphi=velocityOnGridExact(2,iLat,iphi);
		! !				else
		! 					Vr=zero;Vtheta=velocityOnGrid(1,iLat,iphi); Vphi=velocityOnGrid(2,iLat,iphi);
		! !				end if
		! 				write(xrstring,'(f16.10)' )sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
		! 				write(yrstring,'(f16.10)' )sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
		! 				write(zrstring,'(f16.10)' )cosTheta*Vr - sinTheta*Vtheta
		! 				!write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
		! 				write(220,*)sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
		! 				write(220,*)sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
		! 				write(220,*)cosTheta*Vr - sinTheta*Vtheta
		! 			end do
		! 		end do
		! 		write(220,*)'        </DataArray>'
		! end if

				write(220,*)'        <DataArray type="Float32" NumberOfComponents="3" Name="velocity Psi Poisson" Format="ascii">'
				do iLat=1,sphericalGrid%nLat
					lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
					sinTheta=sphericalGrid%cosLatC(ilat)
					cosTheta=sphericalGrid%sinLatC(ilat)
					do iphi=1,sphericalGrid%nPhi
						phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi-1,QR_K)
						sinPhi=sphericalGrid%sinPhiC(iphi)
						cosPhi=sphericalGrid%cosPhiC(iphi)
		!				if(kkk.eq.9999)then
		!					Vr=zero;Vtheta=velocityOnGridExact(1,iLat,iphi); Vphi=velocityOnGridExact(2,iLat,iphi);
		!				else
							Vr=zero;Vtheta=velocityOnGrid(1,iLat,iphi); Vphi=velocityOnGrid(2,iLat,iphi);
		!				end if
						write(xrstring,'(f16.10)' )sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
						write(yrstring,'(f16.10)' )sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
						write(zrstring,'(f16.10)' )cosTheta*Vr - sinTheta*Vtheta
						!write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
						write(220,*)sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
						write(220,*)sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
						write(220,*)cosTheta*Vr - sinTheta*Vtheta
					end do
				end do
				write(220,*)'        </DataArray>'

		! if(1.eq.2)then
		! 		write(220,*)'        <DataArray type="Float32" NumberOfComponents="3" Name="velocityExact" Format="ascii">'
		! 		do iLat=1,sphericalGrid%nLat
		! 			lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
		! 			sinTheta=sphericalGrid%cosLatC(ilat)
		! 			cosTheta=sphericalGrid%sinLatC(ilat)
		! 			do iphi=1,sphericalGrid%nPhi
		! 				phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iphi-1,QR_K)
		! 				sinPhi=sphericalGrid%sinPhiC(iphi)
		! 				cosPhi=sphericalGrid%cosPhiC(iphi)
		! !				if(kkk.eq.9999)then
		! 					Vr=zero;Vtheta=velocityOnGridExact(1,iLat,iphi); Vphi=velocityOnGridExact(2,iLat,iphi);
		! !				else
		! !					Vr=zero;Vtheta=velocityOnGrid(1,iLat,iphi); Vphi=velocityOnGrid(2,iLat,iphi);
		! !				end if
		! 				write(xrstring,'(f16.10)' )sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
		! 				write(yrstring,'(f16.10)' )sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
		! 				write(zrstring,'(f16.10)' )cosTheta*Vr - sinTheta*Vtheta
		! 				!write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
		! 				write(220,*)sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
		! 				write(220,*)sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
		! 				write(220,*)cosTheta*Vr - sinTheta*Vtheta
		! 			end do
		! 		end do
		! 		write(220,*)'        </DataArray>'
		! end if
		write(220,*)'      </CellData>'

		write(220,*)'      <Points>'
		write(220,*)'        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
		do iLat=1,sphericalGrid%nLat+1
			do iphi=1,sphericalGrid%nPhi
				write(xrstring,'(f16.10)' )x(iLat,iphi)
				write(yrstring,'(f16.10)' )y(iLat,iphi)
				write(zrstring,'(f16.10)' )z(iLat,iphi)
				write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
			end do
		end do
		write(220,*)'        </DataArray>'
		write(220,*)'      </Points>'

		write(220,*)'     <Cells>'
		write(220,*)'        <DataArray type="Int32" Name="connectivity" Format="ascii">'

		do iLat=1,sphericalGrid%nLat
			do iphi=1,sphericalGrid%nPhi
				iphiP1=iphi+1
				if(iphi.eq.sphericalGrid%nPhi)iphiP1=1
				indxLL=(iLat-1)*sphericalGrid%nPhi+iphi
				indxLR=(iLat-1)*sphericalGrid%nPhi+iphiP1
				indxUR=iLat*sphericalGrid%nPhi+iphiP1
				indxUL=iLat*sphericalGrid%nPhi+iphi
				write(220,*)indxLL-1,' ',indxLR-1,' ',indxUL-1,' ',indxUR-1
			end do
		end do
		write(220,*)'	        </DataArray>'
		write(220,*)'        <DataArray type="Int32" Name="offsets" Format="ascii">'
		do i=1,sphericalGrid%nLat*sphericalGrid%nPhi
			write(xstring,'(i10)' )i*4
			write(220,*)trim(adjustl(xstring))
		end do
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="UInt8" Name="types" Format="ascii">'
		do i=1,sphericalGrid%nLat*sphericalGrid%nPhi
			write(xstring,'(i10)' )8
			write(220,*)trim(adjustl(xstring))
		end do
		write(220,*)'        </DataArray>'
		write(220,*)'      </Cells>'

		write(220,*)'    </Piece>'

		write(220,*)'  </UnstructuredGrid>'

		write(220,*)'</VTKFile>'
		close(220)

		deallocate(x)
		deallocate(y)
		deallocate(z)

	end subroutine saveMesh

	!ADDED used for debugging a single value.
	! subroutine saveMesh()
	! 	use basicDataStructures
	! 	use simulationParameters
	! 	use sphericalGridDataModule
	! 	! use particlesDataModule
	! 	implicit none 
	! 	CHARACTER  :: STRING*80
	! 	CHARACTER(len=32) :: xrstring
	! 	CHARACTER(len=32) :: yrstring
	! 	CHARACTER(len=32) :: zrstring
	! 	CHARACTER(len=10) :: xstring
	! 	CHARACTER(len=10) :: ystring
	! 	integer :: i,indxLL,indxLR,indxUR,indxUL, iLat, iphi,iphiP1
	! 	real (kind=QR_K) , allocatable, dimension(:,:) :: x, y, z
	! 	real (kind=QR_K) :: phi, sinTheta, cosTheta
	! 	real (kind=QR_K) :: sinPhi, cosPhi, lat, Vr, Vtheta, Vphi

	! 	if(.NOT.allocated(psiOnGrid))then
	! 		allocate(psiOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
	! 		psiOnGrid=zero
	! 	end if

	! 	allocate(x(sphericalGrid%nLat+1,sphericalGrid%nPhi))
	! 	allocate(y(sphericalGrid%nLat+1,sphericalGrid%nPhi))
	! 	allocate(z(sphericalGrid%nLat+1,sphericalGrid%nPhi))

	! 	do iLat=1,sphericalGrid%nLat+1
	! 		lat=-half*pi+sphericalGrid%dLat*real(iLat-1,QR_K)
	! 		do iphi=1,sphericalGrid%nPhi
	! 			phi=sphericalGrid%dphi*real(iphi-1,QR_K)
	! 			x(iLat,iphi)=sphericalGrid%cosLat(ilat)*sphericalGrid%cosPhi(iphi)
	! 			y(iLat,iphi)=sphericalGrid%cosLat(ilat)*sphericalGrid%sinPhi(iphi) 
	! 			z(iLat,iphi)=sphericalGrid%sinLat(ilat)
	! 		end do
	! 	end do

	! 	open(220,file='../results/meshPsi_OLD.vtu')

	! 	write(220,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
	! 	write(220,*)'  <UnstructuredGrid>'
		
	! 	write(xstring,'(i10)' )sphericalGrid%nLat*sphericalGrid%nPhi
	! 	write(ystring,'(i10)' )(sphericalGrid%nLat+1)*sphericalGrid%nPhi
	! 	write(220,*)'    <Piece NumberOfPoints="'//trim(adjustl(ystring))//'" NumberOfCells="'//trim(adjustl(xstring))//'">'

	! 	write(220,*)'      <CellData Scalars="scalars">'

	! 	write(220,*)'        <DataArray type="Float32" Name="Psi OLD" Format="ascii">'
	! 	do iLat=1,sphericalGrid%nLat
	! 		do iphi=1,sphericalGrid%nPhi
	! 			write(xrstring,'(e20.10)' ) psiOnGrid(iLat,iphi) 
	! 			write(220,*) xrstring
	! 		end do
	! 	end do
	! 	write(220,*)'        </DataArray>'

	! 	write(220,*)'      </CellData>'

	! 	write(220,*)'      <Points>'
	! 	write(220,*)'        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
	! 	do iLat=1,sphericalGrid%nLat+1
	! 		do iphi=1,sphericalGrid%nPhi
	! 			write(xrstring,'(f16.10)' )x(iLat,iphi)
	! 			write(yrstring,'(f16.10)' )y(iLat,iphi)
	! 			write(zrstring,'(f16.10)' )z(iLat,iphi)
	! 			write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
	! 		end do
	! 	end do
	! 	write(220,*)'        </DataArray>'
		
	! 	write(220,*)'      </Points>'

	! 	write(220,*)'     <Cells>'
	! 	write(220,*)'        <DataArray type="Int32" Name="connectivity" Format="ascii">'

	! 	do iLat=1,sphericalGrid%nLat
	! 		do iphi=1,sphericalGrid%nPhi
	! 			iphiP1=iphi+1
	! 			if(iphi.eq.sphericalGrid%nPhi)iphiP1=1
	! 			indxLL=(iLat-1)*sphericalGrid%nPhi+iphi
	! 			indxLR=(iLat-1)*sphericalGrid%nPhi+iphiP1
	! 			indxUR=iLat*sphericalGrid%nPhi+iphiP1
	! 			indxUL=iLat*sphericalGrid%nPhi+iphi
	! 			write(220,*)indxLL-1,' ',indxLR-1,' ',indxUL-1,' ',indxUR-1
	! 		end do
	! 	end do
	! 	write(220,*)'	        </DataArray>'
	! 	write(220,*)'        <DataArray type="Int32" Name="offsets" Format="ascii">'
	! 	do i=1,sphericalGrid%nLat*sphericalGrid%nPhi
	! 		write(xstring,'(i10)' )i*4
	! 		write(220,*)trim(adjustl(xstring))
	! 	end do
	! 	write(220,*)'        </DataArray>'

	! 	write(220,*)'        <DataArray type="UInt8" Name="types" Format="ascii">'
	! 	do i=1,sphericalGrid%nLat*sphericalGrid%nPhi
	! 		write(xstring,'(i10)' )8
	! 		write(220,*)trim(adjustl(xstring))
	! 	end do
	! 	write(220,*)'        </DataArray>'
	! 	write(220,*)'      </Cells>'

	! 	write(220,*)'    </Piece>'

	! 	write(220,*)'  </UnstructuredGrid>'

	! 	write(220,*)'</VTKFile>'
	! 	close(220)

	! 	deallocate(x)
	! 	deallocate(y)
	! 	deallocate(z)

	! end subroutine saveMesh

	!ADDED SUBROUTINE: This subroutine is for saving the OceanGrid Data.
	subroutine saveOceanMesh()
		use basicDataStructures
		use simulationParameters
		use OceanGridDataModule
		implicit none 
		CHARACTER  :: STRING*80
		CHARACTER(len=32) :: xrstring
		CHARACTER(len=32) :: yrstring
		CHARACTER(len=32) :: zrstring
		CHARACTER(len=10) :: xstring
		CHARACTER(len=10) :: ystring
		!	CHARACTER(len=10) :: zstring
		integer :: i,indxLL,indxLR,indxUR,indxUL, iLat, iphi,iphiP1
		real (kind=QR_K) , allocatable, dimension(:,:) :: x, y, z
		real (kind=QR_K) :: phi, sinTheta, cosTheta
		real (kind=QR_K) :: sinPhi, cosPhi, lat, Vr, Vtheta, Vphi , OceanDepth 

		if(.NOT.allocated(psi1onGrid))then
			allocate(psi1onGrid(OceanGrid%nLat,OceanGrid%nPhi))
			psi1onGrid=zero
		end if

		! if(.NOT.allocated(phi1onGrid))then
		! 	allocate(phi1onGrid(OceanGrid%nLat,OceanGrid%nPhi))
		! 	phi1onGrid=zero
		! end if

		! if(.NOT.allocated(velocityVorOnGrid))then
		! 	allocate(velocityVorOnGrid(2,OceanGrid%nLat,OceanGrid%nPhi))
			! velocityVorOnGrid=zero
		! end if

		allocate(x(OceanGrid%nLat+1,OceanGrid%nPhi))
		allocate(y(OceanGrid%nLat+1,OceanGrid%nPhi))
		allocate(z(OceanGrid%nLat+1,OceanGrid%nPhi))

		do iLat=1,OceanGrid%nLat+1
			lat=-half*pi+OceanGrid%dLat*real(iLat-1,QR_K)
			do iphi=1,OceanGrid%nPhi
				phi=OceanGrid%dphi*real(iphi-1,QR_K)
				x(iLat,iphi)=OceanGrid%cosLat(ilat)*OceanGrid%cosPhi(iphi)
				y(iLat,iphi)=OceanGrid%cosLat(ilat)*OceanGrid%sinPhi(iphi) 
				z(iLat,iphi)=OceanGrid%sinLat(ilat)
			end do
		end do
		
		!COMMENT: OUTPUT FILE INTERPOLATED_DEPTH.VTU FOR SAVING THE INTERPOLATED DEPTH AND COMPARE WITH THE PYTHON INTERPOLATION TO SECURE THE PROCESS ON BOTH ENDS.
		! open(220,file='../results/INTERPOLATED_DEPTH.vtu')

		!COMMENT: OUTPUT FILE meshPsi.vtu 
		open(220,file='../results/meshPsi.vtu')

		write(220,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
		write(220,*)'  <UnstructuredGrid>'
		
		write(xstring,'(i10)' )OceanGrid%nLat*OceanGrid%nPhi
		write(ystring,'(i10)' )(OceanGrid%nLat+1)*OceanGrid%nPhi
		write(220,*)'    <Piece NumberOfPoints="'//trim(adjustl(ystring))//'" NumberOfCells="'//trim(adjustl(xstring))//'">'

		write(220,*)'      <CellData Scalars="scalars">'
		
		!COMMENT: ψ₁ for the eq: ∇ψ₁ − (1/Ld²) ψ₁ = − f(θ) 
		write(220,*)'        <DataArray type="Float32" Name="ψ₁" Format="ascii">'
		do iLat=1,OceanGrid%nLat
			do iphi=1,OceanGrid%nPhi
				write(xrstring,'(e20.10)') psi1onGrid(iLat,iphi) 
				write(220,*) xrstring
			end do
		end do
		write(220,*)'        </DataArray>'

		!COMMENT: depth d₀
		write(220,*)'        <DataArray type="Float32" Name="d₀" Format="ascii">'
		do iLat=1,OceanGrid%nLat
			do iphi=1,OceanGrid%nPhi
				write(220,*) OceanGrid%d(iLat,iphi)
				end do
			end do
		write(220,*)'        </DataArray>'

		!COMMENT: Latitudes θ
		write(220,*)'        <DataArray type="Float32" Name="lat" Format="ascii">'
		! θ
		do iLat=1,OceanGrid%nLat
			do iphi=1,OceanGrid%nPhi
				write(220,*) -half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)
				end do
			end do
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="Float32" Name="long" Format="ascii">'
		do iLat=1,OceanGrid%nLat
			do iphi=1,OceanGrid%nPhi
				write(220,*) half*OceanGrid%dphi+OceanGrid%dphi*real(iphi-1,QR_K)
				end do
			end do
		write(220,*)'        </DataArray>'
		
		!COMMENT: Coriolis Force f=2Ωsinθ
		write(220,*)'        <DataArray type="Float32" Name="f(θ)" Format="ascii">'
			do iLat=1,OceanGrid%nLat
				lat=-half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)
				do iphi=1,OceanGrid%nPhi
					phi=half*OceanGrid%dphi+OceanGrid%dphi*real(iphi-1,QR_K)
					write(xrstring,'(f16.10)' ) two*OMEGA*OceanGrid%sinLatC(ilat)
					write(220,*)xrstring
				end do
			end do
		write(220,*)'        </DataArray>'
		
		!COMMENT: u₁(θ,ƛ) : Velocity from ψ₁  
			! write(220,*)'       <DataArray type="Float64" NumberOfComponents="3" Name="u₁(θ,ƛ)" Format="ascii">'
			! do iLat=1,OceanGrid%nLat
			! 	lat=-half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)
			! 	sinTheta=OceanGrid%cosLatC(ilat)
			! 	cosTheta=OceanGrid%sinLatC(ilat)
			! 	do iphi=1,OceanGrid%nPhi
			! 		phi=half*OceanGrid%dphi+OceanGrid%dphi*real(iphi-1,QR_K)
			! 		sinPhi=OceanGrid%sinPhiC(iphi)
			! 		cosPhi=OceanGrid%cosPhiC(iphi)

			! 		Vr=zero; Vtheta=velocityPsiOneOnGrid(1,iLat,iphi); Vphi=velocityPsiOneOnGrid(2,iLat,iphi);

			! 		write(xrstring,'(e20.10)' )sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
			! 		write(yrstring,'(e20.10)' )sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
			! 		write(zrstring,'(e20.10)' )cosTheta*Vr - sinTheta*Vtheta
			! 		!write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
			! 		write(220,'(e20.10)')sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
			! 		write(220,'(e20.10)')sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
			! 		write(220,'(e20.10)')cosTheta*Vr - sinTheta*Vtheta
			! 	end do
			! end do
			! write(220,*)'        </DataArray>'

			write(220,*)'      </CellData>'

			write(220,*)'      <Points>'
			write(220,*)'        <DataArray type="Float64" NumberOfComponents="3" Format="ascii">'
			do iLat=1,OceanGrid%nLat+1
				do iphi=1,OceanGrid%nPhi
					write(xrstring,'(e20.10)' )x(iLat,iphi)
					write(yrstring,'(e20.10)' )y(iLat,iphi)
					write(zrstring,'(e20.10)' )z(iLat,iphi)
					write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
				end do
			end do
			write(220,*)'        </DataArray>'
			write(220,*)'      </Points>'

			write(220,*)'     <Cells>'
			write(220,*)'        <DataArray type="Int32" Name="connectivity" Format="ascii">'

			do iLat=1,OceanGrid%nLat
				do iphi=1,OceanGrid%nPhi
					iphiP1=iphi+1
					if(iphi.eq.OceanGrid%nPhi)iphiP1=1
					indxLL=(iLat-1)*OceanGrid%nPhi+iphi
					indxLR=(iLat-1)*OceanGrid%nPhi+iphiP1
					indxUR=iLat*OceanGrid%nPhi+iphiP1
					indxUL=iLat*OceanGrid%nPhi+iphi
					write(220,*)indxLL-1,' ',indxLR-1,' ',indxUL-1,' ',indxUR-1
				end do
			end do
			write(220,*)'	        </DataArray>'
			write(220,*)'        <DataArray type="Int32" Name="offsets" Format="ascii">'
			do i=1,OceanGrid%nLat*OceanGrid%nPhi
				write(xstring,'(i10)' )i*4
				write(220,*)trim(adjustl(xstring))
			end do
			write(220,*)'        </DataArray>'

			write(220,*)'        <DataArray type="UInt8" Name="types" Format="ascii">'
			do i=1,OceanGrid%nLat*OceanGrid%nPhi
				write(xstring,'(i10)' )8
				write(220,*)trim(adjustl(xstring))
			end do
			write(220,*)'        </DataArray>'
			write(220,*)'      </Cells>'
			write(220,*)'    </Piece>'
			write(220,*)'  </UnstructuredGrid>'
			write(220,*)'</VTKFile>'
			close(220)

		deallocate(x)
		deallocate(y)
		deallocate(z)

	end subroutine saveOceanMesh

	!TOFIX: (in progress... waiting to get Q first ...) this subroutine should save on AdaptiveOceanGrid separately from the OceanGrid. 
	subroutine saveOceanAdaptiveMesh(MeshMesh)
		use basicDataStructures
		use simulationParameters
		use OceanAdaptiveGridDataModule
	
		implicit none 
		integer, intent(in) :: MeshMesh

		CHARACTER  :: STRING*80
		CHARACTER(len=32) :: xrstring
		CHARACTER(len=32) :: yrstring
		CHARACTER(len=32) :: zrstring
		CHARACTER(len=10) :: xstring
		CHARACTER(len=10) :: ystring
		!	CHARACTER(len=10) :: zstring
		integer :: i,indxLL,indxLR,indxUR,indxUL, iLat, iphi,iphiP1, Qzero,Qnonzero
		real (kind=QR_K) , allocatable, dimension(:,:) :: x, y, z
		real (kind=QR_K) :: phi, sinTheta, cosTheta
		real (kind=QR_K) :: sinPhi, cosPhi, lat, Vr, Vtheta, Vphi , OceanDepth 
	
		allocate(x(OceanAdaptiveGrid%nLat+1,OceanAdaptiveGrid%nPhi))
		allocate(y(OceanAdaptiveGrid%nLat+1,OceanAdaptiveGrid%nPhi))
		allocate(z(OceanAdaptiveGrid%nLat+1,OceanAdaptiveGrid%nPhi))
	
		do iLat=1,OceanAdaptiveGrid%nLat+1
			lat=-half*pi+OceanAdaptiveGrid%dLat*real(iLat-1,QR_K)
			do iphi=1,OceanAdaptiveGrid%nPhi
				phi=OceanAdaptiveGrid%dphi*real(iphi-1,QR_K)
				x(iLat,iphi)=OceanAdaptiveGrid%cosLat(ilat)*OceanAdaptiveGrid%cosPhi(iphi)
				y(iLat,iphi)=OceanAdaptiveGrid%cosLat(ilat)*OceanAdaptiveGrid%sinPhi(iphi) 
				z(iLat,iphi)=OceanAdaptiveGrid%sinLat(ilat)
			end do
		end do
		
		WRITE(STRING, FMT='(I5)') MeshMesh
		open(220,file='../results/meshAdapt'//trim(adjustl(STRING))//'.vtu')
		! open(220,file='../results/meshAdapt.vtu')
	
		write(220,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
		write(220,*)'  <UnstructuredGrid>'
	
		write(xstring,'(i10)' )OceanAdaptiveGrid%nLat*OceanAdaptiveGrid%nPhi
		write(ystring,'(i10)' )(OceanAdaptiveGrid%nLat+1)*OceanAdaptiveGrid%nPhi
		write(220,*)'    <Piece NumberOfPoints="'//trim(adjustl(ystring))//'" NumberOfCells="'//trim(adjustl(xstring))//'">'
	
		write(220,*)'      <CellData Scalars="scalars">'
		write(220,*)'        <DataArray type="Float32" Name="d" Format="ascii">'
		do iLat=1,OceanAdaptiveGrid%nLat
			do iphi=1,OceanAdaptiveGrid%nPhi
				write(220,*) OceanAdaptiveGrid%d(iLat,iphi)
				end do
			end do
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="Float32" Name="f" Format="ascii">'
			do iLat=1,OceanAdaptiveGrid%nLat
				lat=-half*pi+half*OceanAdaptiveGrid%dLat+OceanAdaptiveGrid%dLat*real(iLat-1,QR_K)
				do iphi=1,OceanAdaptiveGrid%nPhi
					phi=half*OceanAdaptiveGrid%dphi+OceanAdaptiveGrid%dphi*real(iphi-1,QR_K)
					write(xrstring,'(f16.10)' ) two*OMEGA*OceanAdaptiveGrid%sinLatC(ilat)
					write(220,*)xrstring
				end do
			end do
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="Float32" Name="Q" Format="ascii">'
		Qnonzero=0
		Qzero=0 
		do iLat=1,OceanAdaptiveGrid%nLat
			do iphi=1,OceanAdaptiveGrid%nPhi
				! write(*,*) 'HERE = ', OceanAdaptiveGrid%Q(iLat,iphi)
				write(220,*) OceanAdaptiveGrid%Q(iLat,iphi) 
				if (OceanAdaptiveGrid%Q(iLat,iphi)== zero) then
					Qzero=Qzero+1
				else 
					Qnonzero=Qnonzero+1
				end if 
			end do
		end do
		write(*,*) 'Qnonzero, Qzero = ',Qnonzero,Qzero
		write(220,*)'        </DataArray>'

		write(220,*)'        <DataArray type="Float32" Name="Strength Q" Format="ascii">'
		do iLat=1,OceanAdaptiveGrid%nLat
			do iphi=1,OceanAdaptiveGrid%nPhi
				write(220,*) OceanAdaptiveGrid%strengthQ(iLat,iphi)
			end do
		end do
		write(220,*)'        </DataArray>'

		write(220,*)'      </CellData>'

			write(220,*)'      <Points>'
			write(220,*)'        <DataArray type="Float64" NumberOfComponents="3" Format="ascii">'
			do iLat=1,OceanAdaptiveGrid%nLat+1
				do iphi=1,OceanAdaptiveGrid%nPhi
					write(xrstring,'(e20.10)' )x(iLat,iphi)
					write(yrstring,'(e20.10)' )y(iLat,iphi)
					write(zrstring,'(e20.10)' )z(iLat,iphi)
					write(220,*)trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
				end do
			end do
			write(220,*)'        </DataArray>'
			write(220,*)'      </Points>'
	
			write(220,*)'     <Cells>'
			write(220,*)'        <DataArray type="Int32" Name="connectivity" Format="ascii">'
	
			do iLat=1,OceanAdaptiveGrid%nLat
				do iphi=1,OceanAdaptiveGrid%nPhi
					iphiP1=iphi+1
					if(iphi.eq.OceanAdaptiveGrid%nPhi)iphiP1=1
					indxLL=(iLat-1)*OceanAdaptiveGrid%nPhi+iphi
					indxLR=(iLat-1)*OceanAdaptiveGrid%nPhi+iphiP1
					indxUR=iLat*OceanAdaptiveGrid%nPhi+iphiP1
					indxUL=iLat*OceanAdaptiveGrid%nPhi+iphi
					write(220,*)indxLL-1,' ',indxLR-1,' ',indxUL-1,' ',indxUR-1
				end do
			end do
			write(220,*)'	        </DataArray>'
			
			write(220,*)'        <DataArray type="Int32" Name="offsets" Format="ascii">'
			do i=1,OceanAdaptiveGrid%nLat*OceanAdaptiveGrid%nPhi
				write(xstring,'(i10)' )i*4
				write(220,*)trim(adjustl(xstring))
			end do
			write(220,*)'        </DataArray>'
	
			write(220,*)'        <DataArray type="UInt8" Name="types" Format="ascii">'
			do i=1,OceanAdaptiveGrid%nLat*OceanAdaptiveGrid%nPhi
				write(xstring,'(i10)' )8
				write(220,*)trim(adjustl(xstring))
			end do
			write(220,*)'        </DataArray>'
			write(220,*)'      </Cells>'
	
			write(220,*)'    </Piece>'
	
			write(220,*)'  </UnstructuredGrid>'
	
			write(220,*)'</VTKFile>'
			close(220)
	
		deallocate(x)
		deallocate(y)
		deallocate(z)
	
	end subroutine saveOceanAdaptiveMesh

	!ADDED SUBROUTINE not used: viewing the wind Velocity.
	subroutine windVelocity()
		use  windDataModule
		implicit none
		integer :: iLon,iLat,ocean
		logical :: file_exists

		if(.NOT.allocated(wind%velocity))allocate(wind%velocity(wind%nLong,wind%nLat))

		INQUIRE(FILE='../windCopernicus/wind00/wind_components.csv', EXIST=file_exists)
		if(file_exists)then
			open(UNIT=2,FILE='../windCopernicus/wind00/wind_components.csv')
			read(2,*)
				do iLat=1,wind%nLat
					do iLon=1,wind%nLong
					read(2,*)wind%velocity(iLon,iLat)
				end do
			end do
			close(2)
		else
			write(*,*)'wind file does not exist'
			stop
		end if
		
	end subroutine windVelocity

	!COMMENT: 2d Bathymetry not needed
	subroutine saveBathymetry2D()
        use basicDataStructures
        use simulationParameters
        use bathymetryDataModule
        implicit none
        CHARACTER(len=32) :: xrstring, yrstring, zrstring
        CHARACTER(len=10) :: xstring
        integer :: iLat, iphi, Rb, HodaRb
        real (kind=QR_K) :: dtheta, phi, lat

        ! Initialize grid parameters
        dtheta = two*pi/real(bathymetry%nLong,QR_K)  ! Longitude step

        ! Open VTS file
        open(221, file='../results/bathymetry_2d.vts')

        ! Write VTK StructuredGrid header
        write(221,*) '<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian">'
        write(221,*) '  <StructuredGrid WholeExtent="0 ',bathymetry%nLong-1,' 0 ',bathymetry%nLat-1,' 0 0">'
        write(xstring,'(i10)') bathymetry%nLat*bathymetry%nLong
        write(221,*) '<Piece Extent="0 ',bathymetry%nLong-1,& 
		' 0 ',bathymetry%nLat-1,&
		' 0 0" NumberOfPoints="',trim(adjustl(xstring)),'">'

        ! Write PointData
        write(221,*) '      <PointData Scalars="scalars">'
        
        write(221,*) '        <DataArray type="Int32" Name="Z" Format="ascii">'
        do iLat=1,bathymetry%nLat
            do iphi=1,bathymetry%nLong
                write(xrstring,'(I7)') bathymetry%Z(iPhi,iLat)
                write(221,*) xrstring
            end do
        end do
        write(221,*) '        </DataArray>'

        write(221,*) '        <DataArray type="Int32" Name="RossbyRadius" Format="ascii">'
        do iLat=1,bathymetry%nLat
            lat = -half*pi + dtheta*real(iLat-1,QR_K)
            do iphi=1,bathymetry%nLong
                Rb = 0
                if (bathymetry%Z(iPhi,iLat).lt.0) then
                    Rb = int(sqrt(abs(real(bathymetry%Z(iPhi,iLat),QR_K))*real(9.81,QR_K))/(two*OMEGA*abs(sin(lat))*thousand))
                end if
                if (Rb.gt.10000) Rb = 10000
                write(xrstring,'(I7)') Rb
                write(221,*) xrstring
            end do
        end do
        write(221,*) '        </DataArray>'

        write(221,*) '        <DataArray type="Int32" Name="Hoda RossbyRadius" Format="ascii">'
        do iLat=1,bathymetry%nLat
            lat = -half*pi + dtheta*real(iLat-1,QR_K)
            do iphi=1,bathymetry%nLong
                HodaRb = 0
                if (bathymetry%Z(iPhi,iLat).lt.0) then
                    HodaRb = int(sqrt(abs(real(bathymetry%Z(iPhi,iLat),QR_K))*real(9.81,QR_K))/(two*OMEGA*abs(sin(lat))*1000))
                end if
                if (HodaRb.gt.10000) HodaRb = 10000
                write(xrstring,'(I7)') HodaRb
                write(221,*) xrstring
            end do
        end do
        write(221,*) '        </DataArray>'

        write(221,*) '      </PointData>'

        ! Write Points (2D: long, lat, 0)
        write(221,*) '      <Points>'
        write(221,*) '        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
        do iLat=1,bathymetry%nLat
            lat = -half*pi + dtheta*real(iLat-1,QR_K)  ! Latitude from -pi/2 to pi/2
            do iphi=1,bathymetry%nLong
                phi = -pi + dtheta*real(iphi-1,QR_K)  ! Longitude from -pi to pi
                write(xrstring,'(f16.10)') phi  ! Longitude as X
                write(yrstring,'(f16.10)') lat  ! Latitude as Y
                write(zrstring,'(f16.10)') 0.0
                write(221,*) trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
            end do
        end do
        write(221,*) '        </DataArray>'
        write(221,*) '      </Points>'

        ! Close VTK file
        write(221,*) '    </Piece>'
        write(221,*) '  </StructuredGrid>'
        write(221,*) '</VTKFile>'
        close(221)

    end subroutine saveBathymetry2D

	!COMMENT: 2d OceanMesh not needed
	subroutine saveOceanMesh2D()
        use basicDataStructures
        use simulationParameters
        use OceanGridDataModule
        implicit none 
        CHARACTER(len=32) :: xrstring, yrstring, zrstring
        CHARACTER(len=10) :: xstring
        integer :: iLat, iphi
        real (kind=QR_K) :: phi, sinTheta, cosTheta, lat, Vr, Vtheta, Vphi, OceanDepth
        real (kind=QR_K) :: dtheta, dphi, sinPhi, cosPhi

        ! Initialize grid parameters
        dtheta = pi/real(OceanGrid%nLat,QR_K)  ! Latitude step
        dphi = two*pi/real(OceanGrid%nPhi,QR_K)  ! Longitude step

        ! Ensure arrays are allocated
        if (.NOT. allocated(phi1onGrid)) then
            allocate(phi1onGrid(OceanGrid%nLat,OceanGrid%nPhi))
            phi1onGrid = zero
        end if
        if (.NOT. allocated(psi1onGrid)) then
            allocate(psi1onGrid(OceanGrid%nLat,OceanGrid%nPhi))
            psi1onGrid = zero
        end if
        if (.NOT. allocated(velocityVorOnGrid)) then
            allocate(velocityVorOnGrid(2,OceanGrid%nLat,OceanGrid%nPhi))
            velocityVorOnGrid = zero
        end if
        if (.NOT. allocated(velocityPsiOneOnGrid)) then
            allocate(velocityPsiOneOnGrid(2,OceanGrid%nLat,OceanGrid%nPhi))
            velocityPsiOneOnGrid = zero
        end if
        if (.NOT. allocated(velocityPhiOneOnGrid)) then
            allocate(velocityPhiOneOnGrid(2,OceanGrid%nLat,OceanGrid%nPhi))
            velocityPhiOneOnGrid = zero
        end if

        ! Open VTS file
        open(220, file='../results/oceanmesh_2d.vts')

        ! Write VTK StructuredGrid header
        write(220,*) '<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian">'
        write(220,*) '  <StructuredGrid WholeExtent="0 ',OceanGrid%nPhi-1,' 0 ',OceanGrid%nLat-1,' 0 0">'
        write(xstring,'(i10)') OceanGrid%nLat*OceanGrid%nPhi
        write(220,*) '    <Piece Extent="0 &
		',OceanGrid%nPhi-1,' 0 ',OceanGrid%nLat-1,&
		' 0 0" NumberOfPoints="',trim(adjustl(xstring)),'">'

        ! Write PointData
        write(220,*) '      <PointData Scalars="scalars">'

        write(220,*) '        <DataArray type="Float32" Name="Psi1 Poisson" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            do iphi=1,OceanGrid%nPhi
                write(xrstring,'(e20.10)') psi1onGrid(iLat,iphi)
                write(220,*) xrstring
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '        <DataArray type="Float32" Name="Phi1 Poisson" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            do iphi=1,OceanGrid%nPhi
                write(xrstring,'(e20.10)') phi1onGrid(iLat,iphi)
                write(220,*) xrstring
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '        <DataArray type="Float32" Name="depth" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            do iphi=1,OceanGrid%nPhi
                write(xrstring,'(f16.10)') OceanGrid%d(iLat,iphi)
                write(220,*) xrstring
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '        <DataArray type="Float32" Name="Oceandepth" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            do iphi=1,OceanGrid%nPhi
                if (OceanGrid%d(iLat,iphi) > 0.0) then
                    OceanDepth = OceanGrid%d(iLat,iphi)
                else
                    OceanDepth = 0.0
                end if
                write(xrstring,'(f16.10)') OceanDepth
                write(220,*) xrstring
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '        <DataArray type="Float32" Name="f" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            lat = -half*pi + half*OceanGrid%dLat + OceanGrid%dLat*real(iLat-1,QR_K)
            do iphi=1,OceanGrid%nPhi
                write(xrstring,'(f16.10)') two*OMEGA*OceanGrid%sinLatC(ilat)
                write(220,*) xrstring
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '        <DataArray type="Float32" Name="h1" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            do iphi=1,OceanGrid%nPhi
                write(xrstring,'(f16.10)') OceanGrid%h1(iLat,iphi)
                write(220,*) xrstring
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '        <DataArray type="Float32" Name="hs1" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            do iphi=1,OceanGrid%nPhi
                write(xrstring,'(f16.10)') OceanGrid%hs1(iLat,iphi)
                write(220,*) xrstring
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '        <DataArray type="Float32" NumberOfComponents="3" Name="u_psi_1" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            lat = -half*pi + half*OceanGrid%dLat + OceanGrid%dLat*real(iLat-1,QR_K)
            sinTheta = OceanGrid%cosLatC(ilat)
            cosTheta = OceanGrid%sinLatC(ilat)
            do iphi=1,OceanGrid%nPhi
                phi = -pi + dphi*real(iphi-1,QR_K)
                sinPhi = sin(phi)
                cosPhi = cos(phi)
                Vr = zero
                Vtheta = velocityPsiOneOnGrid(1,iLat,iphi)
                Vphi = velocityPsiOneOnGrid(2,iLat,iphi)
                write(xrstring,'(f16.10)') sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
                write(yrstring,'(f16.10)') sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
                write(zrstring,'(f16.10)') cosTheta*Vr - sinTheta*Vtheta
                write(220,*) trim(adjustl(xrstring))
                write(220,*) trim(adjustl(yrstring))
                write(220,*) trim(adjustl(zrstring))
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '        <DataArray type="Float32" NumberOfComponents="3" Name="u_phi_1" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            lat = -half*pi + half*OceanGrid%dLat + OceanGrid%dLat*real(iLat-1,QR_K)
            sinTheta = OceanGrid%cosLatC(ilat)
            cosTheta = OceanGrid%sinLatC(ilat)
            do iphi=1,OceanGrid%nPhi
                phi = -pi + dphi*real(iphi-1,QR_K)
                sinPhi = sin(phi)
                cosPhi = cos(phi)
                Vr = zero
                Vtheta = velocityPhiOneOnGrid(1,iLat,iphi)
                Vphi = velocityPhiOneOnGrid(2,iLat,iphi)
                write(xrstring,'(f16.10)') sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
                write(yrstring,'(f16.10)') sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
                write(zrstring,'(f16.10)') cosTheta*Vr - sinTheta*Vtheta
                write(220,*) trim(adjustl(xrstring))
                write(220,*) trim(adjustl(yrstring))
                write(220,*) trim(adjustl(zrstring))
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '        <DataArray type="Float32" NumberOfComponents="3" Name="u_1" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            lat = -half*pi + half*OceanGrid%dLat + OceanGrid%dLat*real(iLat-1,QR_K)
            sinTheta = OceanGrid%cosLatC(ilat)
            cosTheta = OceanGrid%sinLatC(ilat)
            do iphi=1,OceanGrid%nPhi
                phi = -pi + dphi*real(iphi-1,QR_K)
                sinPhi = sin(phi)
                cosPhi = cos(phi)
                Vr = zero
                Vtheta = velocityOnGrid(1,iLat,iphi)
                Vphi = velocityOnGrid(2,iLat,iphi)
                write(xrstring,'(f16.10)') sinTheta*cosPhi*Vr + cosTheta*cosPhi*Vtheta - sinPhi*Vphi
                write(yrstring,'(f16.10)') sinTheta*sinPhi*Vr + cosTheta*sinPhi*Vtheta + cosPhi*Vphi
                write(zrstring,'(f16.10)') cosTheta*Vr - sinTheta*Vtheta
                write(220,*) trim(adjustl(xrstring))
                write(220,*) trim(adjustl(yrstring))
                write(220,*) trim(adjustl(zrstring))
            end do
        end do
        write(220,*) '        </DataArray>'

        write(220,*) '      </PointData>'

        ! Write Points (2D: long, lat, 0)
        write(220,*) '      <Points>'
        write(220,*) '        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
        do iLat=1,OceanGrid%nLat
            lat = -half*pi + dtheta*real(iLat-1,QR_K)  ! Latitude from -pi/2 to pi/2
            do iphi=1,OceanGrid%nPhi
                phi = -pi + dphi*real(iphi-1,QR_K)  ! Longitude from -pi to pi
                write(xrstring,'(f16.10)') phi  ! Longitude as X
                write(yrstring,'(f16.10)') lat  ! Latitude as Y
                write(zrstring,'(f16.10)') 0.0
                write(220,*) trim(adjustl(xrstring))//' '//trim(adjustl(yrstring))//' '//trim(adjustl(zrstring))
            end do
        end do
        write(220,*) '        </DataArray>'
        write(220,*) '      </Points>'

        ! Close VTK file
        write(220,*) '    </Piece>'
        write(220,*) '  </StructuredGrid>'
        write(220,*) '</VTKFile>'
        close(220)
        
    end subroutine saveOceanMesh2D

end module outputModule

