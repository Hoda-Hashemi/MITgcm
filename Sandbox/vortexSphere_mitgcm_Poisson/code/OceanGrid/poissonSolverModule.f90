module poissonSolverModule
	!COMMENT: the subroutines are oreded the same order here.
	implicit none
	!COMMENT: Original Legacy Code — Some obsolete write statements or commented-out lines may remain.
	public :: buildCoefficientsMatrixAndSolve
	public :: buildCoefficientsMatrix
	public :: solvePoissonForPsi
	public :: solvePoissonForPhi
	public :: removeCoefficientsMatrix

	!COMMENT: The following subroutines are added and not all of them are finalized.
	public :: OceanbuildCoefficientsMatrixPsiOne !!FINALIZED
	public :: solvePoissonForPsiOne !!FINALIZED  
	public :: VelocityFieldsOneOnGrid !!FINALIZED
	! public :: SolvingLibrary 

	public :: OceanbuildCoefficientsMatrixPsiTwo
	public :: solvePoissonForPsiTwo

	public :: OceanbuildCoefficientsMatrixPhiOne 
	public :: solvePoissonForPhiOne

	public :: VelocityFieldsTwoOnGrid
	public :: Compute_RelativeVorticity_Zeta
contains 
	subroutine buildCoefficientsMatrixAndSolve()
		use sphericalGridDataModule
		use sparseMatricesDataModule
		use basicDataStructures
		use simulationParameters
		use lsqrModule
		use GMRESSettings
		use bathymetryDataModule
		!$ use omp_lib
		implicit none
		! the coefficient matrix of the operator Laplacian_s(Psi) - f(theta,phi) Psi/Ld^2
		! in spherical coordinates
		! R^{-2} [1/sin(theta) diff(sin(theta) diff(psi,theta),theta) + 1/sin^2(theta) diff(diff(psi,phi),phi)] - Ld^{-2} f(theta,phi) psi = RHS0
		!
		! [1/sin(theta) diff(sin(theta) diff(psi,theta),theta) + 1/sin^2(theta) diff(diff(psi,phi),phi)] - (R/Ld0)^2 f(theta,phi) psi = R^2 RHS0
		!
		! 1/sin(theta_{ij}) [sin(theta_{i,j+1/2}) [psi(i,j+1) - psi(i,j)]/Dtheta - sin(theta_{i,j-1/2}) [psi(i,j) - psi(i,j-1)]/Dtheta]/Dtheta
		! + 1/sin^2(theta_{ij}) [psi(i+1,j) - 2 psi(i,j) + psi(i-1,j)]/Dphi^2 - (R/Ld0)^2 f(i,j) psi(i,j) = R^2 RHS0
		!
		! coefficients
		! psi(i,j) : -1/sin(theta_{ij}) [sin(theta_{i,j+1/2}) + sin(theta_{i,j-1/2}]/Dtheta^2 - 2/sin^2(theta_{ij})/Dphi^2 - (R/Ld0)^2 f(i,j)
		!            -2/Dtheta^2 - 2/sin^2(theta_{ij})/Dphi^2 - (R/Ld0)^2 f(i,j)
		! psi(i-1,j): 1/sin^2(theta_{ij}) 1/Dphi^2
		! psi(i+1,j): 1/sin^2(theta_{ij}) 1/Dphi^2
		! psi(i,j-1): sin(theta_{i,j-1/2})/sin(theta_{ij}) 1/Dtheta^2
		! psi(i,j+1): sin(theta_{i,j+1/2})/sin(theta_{ij}) 1/Dtheta^2
		!
		! if dtheta=dphi, the new right hand side is dtheta^2 R^2 RHS0
		!
		! psi(i,j) : -2 - 2/sin^2(theta_{ij}) - (R/Ld0)^2 dtheta^2 f(i,j)
		! psi(i-1,j): 1/sin^2(theta_{ij})
		! psi(i+1,j): 1/sin^2(theta_{ij})
		! psi(i,j-1): sin(theta_{i,j-1/2})/sin(theta_{ij})
		! psi(i,j+1): sin(theta_{i,j+1/2})/sin(theta_{ij})
		! RHS = dtheta^2 R^2 RHS0
		!
		! If using the latitude
		!
		! R^{-2} [1/sin(theta) diff(sin(theta) diff(psi,theta),theta) + 1/sin^2(theta) diff(diff(psi,phi),phi)] - Ld^{-2} f(phi,theta) psi = RHS0
		! R^{-2} [1/cos(lat) diff(cos(lat) diff(psi,lat),lat) + 1/cos^2(theta) diff(diff(psi,phi),phi)] - Ld^{-2} f(phi,theta) psi = RHS0
		!
		! psi(i,j) : -2 - 2/cos^2(lat_{ij}) - (R/Ld0)^2 dtheta^2 f(ij)
		! psi(i-1,j): 1/cos^2(lat_{ij})
		! psi(i+1,j): 1/cos^2(lat_{ij})
		! psi(i,j-1): cos(lat_{i,j-1/2})/cos(lat_{ij}) = 1 + (dtheta/2) tan(lat_{ij})
		! psi(i,j+1): cos(lat_{i,j+1/2})/cos(lat_{ij}) = 1 - (dtheta/2) tan(lat_{ij})
		! RHS = dtheta^2 R^2 RHS0
		!
		! Multiplying everybody by cos^2(lat_{ij})
		!
		! psi(i,j) : -2 cos^2(lat_{ij}) - 2 - (R/Ld0)^2 dtheta^2 f(ij) cos^2(lat_{ij})
		! psi(i-1,j): 1
		! psi(i+1,j): 1
		! psi(i,j-1): cos(lat_{ij}) cos(lat_{i,j-1/2})
		! psi(i,j+1): cos(lat_{ij}) cos(lat_{i,j+1/2})
		! RHS = dtheta^2 R^2 cos^2(lat_{ij})	RHS0
		!

		integer :: i,k,iLat,iPhi,iM1,iP1,ki,m,nElements
		real (kind=QR_K) :: ROverLd02, dtheta, dtheta2,lat,cosLat,sinLat,tanLat,phi
		real (kind=QR_K) :: oneOverCosLat2,a,R2,aOverR2,aOverLd02,factor,psiOvera,Qminusf
		real (kind=QR_K) :: cosLatMHalf,cosLatPHalf,cosLat2,aOverR,sinLat2
		double precision, allocatable :: b(:),res(:),sol(:)
		logical :: signConstrained
		integer :: i1,i2,CHUNK,nThreads,threadID
		integer :: iLatP1,iLatM1,iPhiP1,iPhiM1
		double precision :: time1
		real (kind=QR_K) :: d2, lat0, phi0, sigma02, coriolisTerm,dOverd0

		time1=omp_get_wtime()

		write(*,*)'filling matrix',sphericalGrid%nLat,sphericalGrid%nPhi
		AC%length = (sphericalGrid%nLat-2)*sphericalGrid%nPhi*5+2*sphericalGrid%nPhi*4
		allocate(AC%j(AC%length))
		allocate(AC%value(AC%length))
		nElements=sphericalGrid%nLat*sphericalGrid%nPhi
		allocate(AC%nentries(nElements))
		allocate(b(nElements))
		AC%value=zero
		b = zero

		i=0
		k=0
		!		a = a_sim ! a is in m^2/s
		!		m=m_sim
		R2=R*R
		!		aOverR2=a/(R2)
		!		aOverR=a/R
		!		aOverLd02=a/(Ld0*Ld0)
		!		factor = -(real((1+m)*(2+m),QR_K)*aOverR2+aOverLd02) ! factor has the unit /sec

		!   R^2/Ld^2 = R^2/(g d/f^2)
		!            = R^2 f^2/gd
		!            = (2 Omega sinLat)^2 R^2/(gd)
		!            = (4 Omega^2 R^2)/(gd0) sinLat^2/(d/d0)
		!            = ROverLd02 sinLat^2/(d/d0)
		! if d=d0, we get
		!   R^2/Ld^2 = (4 Omega^2 R^2)/(gd0) sinLat^2
		!            = R^2/Ld0^2 sinLat^2
		!            = ROverLd02 sinLat^2
		ROverLd02=(R/Ld0)**2

		dtheta=sphericalGrid%dLat
		dtheta2=dtheta*dtheta
		do iLat=1,sphericalGrid%nLat
			lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
			cosLat=sphericalGrid%cosLatC(iLat); oneOverCosLat2=one/(cosLat*cosLat)
			cosLat2=cosLat*cosLat
			sinLat=sphericalGrid%sinLatC(iLat)
			sinLat2=sinLat*sinLat
			coriolisTerm = two*OMEGA*sinLat !coriolisTerm has the units of /sec
			cosLatMHalf=cos(Lat-half*dtheta)
			cosLatPHalf=cos(Lat+half*dtheta)
			if(iLat.eq.1)cosLatMHalf=zero
			if(iLat.eq.sphericalGrid%nLat)cosLatPHalf=zero

			tanLat=sinLat/cosLat
			do iPhi=1,sphericalGrid%nPhi
				phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iPhi-1,QR_K)

				if (sphericalGrid%d(iLat,iPhi).gt.zero)then
					dOverd0 = sphericalGrid%d(iLat,iPhi)/sphericalGrid%d0
				else
					dOverd0 = zero
				end if

				i=i+1

				iM1=i-1; if(iPhi.eq.1)iM1=i-1+sphericalGrid%nPhi
				iP1=i+1; if(iPhi.eq.sphericalGrid%nPhi)iP1=i+1-sphericalGrid%nPhi

				if(iLat.gt.1)then
					k=k+1;AC%j(k)=i-sphericalGrid%nPhi;AC%value(k)=dOverd0*cosLat*cosLatMHalf
				end if

				k=k+1; AC%j(k)=iM1; AC%value(k)=dOverd0 !DBLE(oneOverCosLat2) !1/sin^2(theta_{ij})

				k=k+1; ki=k; AC%j(k)=i; AC%value(k)=-DBLE(dOverd0*(two + cosLat*(cosLatMHalf+cosLatPHalf)) + &
				&(ROverLd02*sinLat2)*cosLat2*dtheta2) !DBLE(-two - two*oneOverCosLat2 - ROverLd02*dtheta2) !-2 - 2/sin^2(theta_{ij}) - (R/Ld0)^2 dtheta^2

				k=k+1; AC%j(k)=iP1; AC%value(k)=dOverd0 !DBLE(oneOverCosLat2) !1/sin^2(theta_{ij})

				if(iLat.lt.sphericalGrid%nLat)then
					k=k+1; AC%j(k)=i+sphericalGrid%nPhi;AC%value(k)=dOverd0*cosLat*cosLatPHalf
				end if

				AC%nentries(i) = k
				!psiOvera=-sinLat*(cosLat**m)*cos(real(m,QR_K)*phi) !psi has units of m^2/s, psiOvera is unitless
				!Qminusf=factor*psiOvera !Qminusf has the units of /sec

				!d2 = (lat-lat0)**2+(phi-phi0)**2
				!QOnGrid(iLat,iPhi) = (10.0D0*two*OMEGA)*DEXP(-d2/sigma02)/(pi*sigma02)

				Qminusf=QOnGrid(iLat,iPhi)*dOverd0 - coriolisTerm;
				!				Qminusf= - coriolisTerm;
				!dtheta^2 R^2 cos^2(lat_{ij})	RHS0
				b(i) = DBLE(dOverd0*cosLat2*dtheta2*R2*Qminusf)
				!if(ISNAN(b(i)))then
				!write(*,*)cosLat2,dtheta2,R2,QOnGrid(iLat,iPhi),coriolisTerm
				!stop
				!end if
				!				if(i.eq.1)then
				!					write(*,*)i,ki,AC%j(1:AC%nentries(i)),AC%value(1:AC%nentries(i)),b(i)
				!					stop
				!				end if
			end do
		end do

		AC%m = nElements
		AC%n = nElements
		AC%rowIndexed = 1

		allocate(sol(nElements))
		allocate(res(nElements))
		sol = zero
		res = zero
		signConstrained=.FALSE.
		!		call GMRES_RESTARTED_OMP_RETURNR(AC,b(1:nElements),sol(1:nElements),signConstrained,&
		!		&res(1:nElements),nElements,GMRES_ninner,GMRES_nouter,GMRES_tol,successFlag,normr)

		!BCG(A,b,x,m,n,nIter,res)

		call BICGSTAB_OMP(AC,b(1:nElements),sol(1:nElements),nElements,GMRES_nouter,GMRES_tol)

		if(.NOT.allocated(psiOnGrid))allocate(psiOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		psiOnGrid=zero

		i=0
		do iLat=1,sphericalGrid%nLat
			do iPhi=1,sphericalGrid%nPhi
				i=i+1
				psiOnGrid(iLat,iPhi)=sol(i)
			end do
		end do

		deallocate(b)
		deallocate(sol)
		deallocate(res)
		deallocate(AC%value)
		deallocate(AC%j)
		deallocate(AC%nentries)

		!		do iLat=1,sphericalGrid%nLat
		!			if(abs(QOnGrid(iLat,1)).gt.1.0E-10)then
		!			write(*,*)psiOnGrid(iLat,sphericalGrid%nPhi)-psiOnGrid(iLat,sphericalGrid%nPhi-1),&
		!			&psiOnGrid(iLat,1)-psiOnGrid(iLat,sphericalGrid%nPhi),psiOnGrid(iLat,2)-psiOnGrid(iLat,1)
		!			write(*,*)QOnGrid(iLat,sphericalGrid%nPhi-1),QOnGrid(iLat,sphericalGrid%nPhi),QOnGrid(iLat,1)&
		!			&,QOnGrid(iLat,2)
		!			write(*,*)
		!		end if
		!		end do
		!		stop

		!computing velocity using finite difference of psi
		if(.NOT.allocated(velocityVorOnGrid))allocate(velocityVorOnGrid(2,sphericalGrid%nLat,sphericalGrid%nPhi))
		velocityVorOnGrid=zero

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = sphericalGrid%nLat/nThreads
		write(*,*)CHUNK

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,sphericalGrid,&
		!$OMP &velocityVorOnGrid,psiOnGrid)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		!		write(*,*)threadID,i1,i2
		if(threadID+1.eq.nThreads)i2=sphericalGrid%nLat
		!		do iLat=1,sphericalGrid%nLat
		do iLat=i1,i2
			cosLat=sphericalGrid%cosLatC(iLat);
			sinLat=sphericalGrid%sinLatC(iLat);
			iLatP1=iLat+1; iLatM1=iLat-1
			if(iLat.eq.1)iLatM1=1
			if(iLat.eq.sphericalGrid%nLat)iLatP1=sphericalGrid%nLat
			do iPhi=1,sphericalGrid%nPhi
				phi=half*sphericalGrid%dphi+sphericalGrid%dphi*real(iPhi-1,QR_K)
				iPhiP1=iPhi+1;iPhiM1=iPhi-1
				if(iPhi.eq.1)iPhiM1=sphericalGrid%nPhi
				if(iPhi.eq.sphericalGrid%nPhi)iPhiP1=1

				velocityVorOnGrid(1,iLat,iPhi)=(psiOnGrid(iLat,iPhiP1)-psiOnGrid(iLat,iPhiM1))&
				&/(two*sphericalGrid%dLat*R*cosLat)
				velocityVorOnGrid(2,iLat,iPhi)=-(psiOnGrid(iLatP1,iPhi)-psiOnGrid(iLatM1,iPhi))&
				&/(two*sphericalGrid%dLat*R)
			end do
		end do
		!$OMP END PARALLEL

		write(*,*)'POISSON SOLVER sim time=',omp_get_wtime()-time1

		return

	end subroutine buildCoefficientsMatrixAndSolve

	subroutine buildCoefficientsMatrix()
		use sphericalGridDataModule
		use sparseMatricesDataModule
		use basicDataStructures
		use simulationParameters
		use lsqrModule
		use GMRESSettings
		use bathymetryDataModule
		!$ use omp_lib
		implicit none
		! here, we solve Laplacian_s(Psi) - f(theta,Lon) Psi/Ld^2 = RHS
		! the coefficient matrix of the operator Laplacian_s(Psi) - f(theta,Lon) Psi/Ld^2
		! in spherical coordinates
		! R^{-2} [1/sin(theta) diff(sin(theta) diff(psi,theta),theta) + 1/sin^2(theta) diff(diff(psi,Lon),Lon)] - Ld^{-2} f(theta,Lon) psi = RHS0
		!
		! [1/sin(theta) diff(sin(theta) diff(psi,theta),theta) + 1/sin^2(theta) diff(diff(psi,Lon),Lon)] - (R/Ld0)^2 f(theta,Lon) psi = R^2 RHS0
		!
		! 1/sin(theta_{ij}) [sin(theta_{i,j+1/2}) [psi(i,j+1) - psi(i,j)]/dLat - sin(theta_{i,j-1/2}) [psi(i,j) - psi(i,j-1)]/dLat]/dLat
		! + 1/sin^2(theta_{ij}) [psi(i+1,j) - 2 psi(i,j) + psi(i-1,j)]/DLon^2 - (R/Ld0)^2 f(i,j) psi(i,j) = R^2 RHS0
		!
		! coefficients
		! psi(i,j) : -1/sin(theta_{ij}) [sin(theta_{i,j+1/2}) + sin(theta_{i,j-1/2}]/dLat^2 - 2/sin^2(theta_{ij})/DLon^2 - (R/Ld0)^2 f(i,j)
		!            -2/dLat^2 - 2/sin^2(theta_{ij})/DLon^2 - (R/Ld0)^2 f(i,j)
		! psi(i-1,j): 1/sin^2(theta_{ij}) 1/DLon^2
		! psi(i+1,j): 1/sin^2(theta_{ij}) 1/DLon^2
		! psi(i,j-1): sin(theta_{i,j-1/2})/sin(theta_{ij}) 1/dLat^2
		! psi(i,j+1): sin(theta_{i,j+1/2})/sin(theta_{ij}) 1/dLat^2
		!
		! if dLat=dLon, the new right hand side is dLat^2 R^2 RHS0
		!
		! psi(i,j) : -2 - 2/sin^2(theta_{ij}) - (R/Ld0)^2 dLat^2 f(i,j)
		! psi(i-1,j): 1/sin^2(theta_{ij})
		! psi(i+1,j): 1/sin^2(theta_{ij})
		! psi(i,j-1): sin(theta_{i,j-1/2})/sin(theta_{ij})
		! psi(i,j+1): sin(theta_{i,j+1/2})/sin(theta_{ij})
		! RHS = dLat^2 R^2 RHS0
		!
		! If using the latitude
		!
		! R^{-2} [1/sin(theta) diff(sin(theta) diff(psi,theta),theta) + 1/sin^2(theta) diff(diff(psi,Lon),Lon)] - Ld^{-2} f(Lon,theta) psi = RHS0
		! R^{-2} [1/cos(lat) diff(cos(lat) diff(psi,lat),lat) + 1/cos^2(theta) diff(diff(psi,Lon),Lon)] - Ld^{-2} f(Lon,theta) psi = RHS0
		!
		! psi(i,j) : -2 - 2/cos^2(lat_{ij}) - (R/Ld0)^2 dLat^2 f(ij)
		! psi(i-1,j): 1/cos^2(lat_{ij})
		! psi(i+1,j): 1/cos^2(lat_{ij})
		! psi(i,j-1): cos(lat_{i,j-1/2})/cos(lat_{ij}) = 1 + (dLat/2) tan(lat_{ij})
		! psi(i,j+1): cos(lat_{i,j+1/2})/cos(lat_{ij}) = 1 - (dLat/2) tan(lat_{ij})
		! RHS = dLat^2 R^2 RHS0
		!
		! Multiplying everybody by cos^2(lat_{ij})
		!
		! psi(i,j) : -2 cos^2(lat_{ij}) - 2 - (R/Ld0)^2 dLat^2 f(ij) cos^2(lat_{ij})
		! psi(i-1,j): 1
		! psi(i+1,j): 1
		! psi(i,j-1): cos(lat_{ij}) cos(lat_{i,j-1/2})
		! psi(i,j+1): cos(lat_{ij}) cos(lat_{i,j+1/2})
		! RHS = dLat^2 R^2 cos^2(lat_{ij})	RHS0
		!
		integer :: i,k,iLat,iLon,iM1,iP1,nElements
		real (kind=QR_K) :: ROverLd02, dLat, dLat2,lat,cosLat,sinLat,tanLat,Lon
		real (kind=QR_K) :: oneOverCosLat2,psiOvera,Qminusf
		real (kind=QR_K) :: cosLatMHalf,cosLatPHalf,cosLat2,sinLat2
		double precision, allocatable :: b(:),res(:),sol(:)
		logical :: signConstrained
		integer :: i1,i2,CHUNK,nThreads,threadID
		integer :: iLatP1,iLatM1,iLonP1,iLonM1
		double precision :: time1
		real (kind=QR_K) :: d2, lat0, Lon0, sigma02, coriolisTerm,dOverd0 , oneOverLd02 , R2

		time1=omp_get_wtime()

		write(*,*)'filling matrix for screened poisson',sphericalGrid%nLat,sphericalGrid%nPhi

		AC%length = (sphericalGrid%nLat-2)*sphericalGrid%nPhi*5+2*sphericalGrid%nPhi*4
		if(allocated(AC%j))deallocate(AC%j)
		if(allocated(AC%value))deallocate(AC%value)
		if(allocated(AC%nentries))deallocate(AC%nentries)
		allocate(AC%j(AC%length))
		allocate(AC%value(AC%length))
		nElements=sphericalGrid%nLat*sphericalGrid%nPhi
		allocate(AC%nentries(nElements))
		AC%value=zero

		BC%length = AC%length
		if(allocated(BC%j))deallocate(BC%j)
		if(allocated(BC%value))deallocate(BC%value)
		if(allocated(BC%nentries))deallocate(BC%nentries)
		allocate(BC%j(BC%length))
		allocate(BC%value(BC%length))
		allocate(BC%nentries(nElements))
		BC%value=zero

		i=0
		k=0
		!a = a_sim ! a is in m^2/s
		!m=m_sim
		!R2=R*R
		!aOverR2=a/(R2)
		!aOverR=a/R
		!aOverLd02=a/(Ld*Ld)
		!factor = -(real((1+m)*(2+m),QR_K)*aOverR2+aOverLd02) ! factor has the unit /sec

		!ADDED HODA: Ld0 =  sqrt(g d0)/(2 Ω)
		!COMMENT: initial back up code Ld0 is evaluated as: 
		R2 = R*R
		ROverLd02=(R/Ld0)**2
		oneOverLd02 = (one/Ld0)**2

		dLat=sphericalGrid%dLat
		dLat2=dLat*dLat
		do iLat=1,sphericalGrid%nLat
			lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
			cosLat=sphericalGrid%cosLatC(iLat);
			cosLat2=cosLat*cosLat
			! oneOverCosLat2=one/cosLat2
			sinLat=sphericalGrid%sinLatC(iLat)
			sinLat2=sinLat*sinLat

			cosLatMHalf=cos(Lat-half*dLat)
			cosLatPHalf=cos(Lat+half*dLat)
			if(iLat.eq.1)cosLatMHalf=zero
			if(iLat.eq.sphericalGrid%nLat)cosLatPHalf=zero

			! tanLat=sinLat/cosLat
			do iLon=1,sphericalGrid%nPhi
				Lon=half*sphericalGrid%dPhi+sphericalGrid%dPhi*real(iLon-1,QR_K)
				if (sphericalGrid%d(iLat,iLon).gt.zero)then
					dOverd0 = sphericalGrid%d(iLat,iLon)/sphericalGrid%d0
				else
					dOverd0 = zero
				end if

				i=i+1
				iM1=i-1; if(iLon.eq.1)iM1=i-1+sphericalGrid%nPhi
				iP1=i+1; if(iLon.eq.sphericalGrid%nPhi)iP1=i+1-sphericalGrid%nPhi

				!!!ALTERED: SCALED!!
				!COMMENT: ψ₁(ƛᵢ,θᵢ₋₁)
				if(iLat.gt.1)then  
					k=k+1;
					AC%j(k)=i-sphericalGrid%nPhi;
					AC%value(k)=dOverd0*cosLat*cosLatMHalf 
					BC%j(k)=i-sphericalGrid%nPhi;BC%value(k)=cosLat*cosLatMHalf 
				end if
				!COMMENT: ψ₁(ƛᵢ₋₁, θᵢ)
				k=k+1;
				AC%j(k)=iM1; AC%value(k)=dOverd0   !DBLE(oneOverCosLat2) !1/sin^2(theta_{ij})
				BC%j(k)=iM1; BC%value(k)=one
				
				!COMMENT: ψ₁(ƛᵢ,θᵢ) !the very original input sin is variable Ld(θ): AC%value(k)=-DBLE(dOverd0*(two + cosLat*(cosLatMHalf+cosLatPHalf)) + (ROverLd02*sinLat2)*cosLat2*dLat2)
				k=k+1; 
				AC%j(k)=i; 
				!!SCALED 8. constant Ld- depth and θ
				AC%value(k)=-DBLE(dOverd0*(two + cosLat*(cosLatMHalf+cosLatPHalf)) + (ROverLd02*sinTheta02)*cosLat2*dLat2) 
				
				!!SCALED Forcing the screening term Ld to be constant 7. 6. 
				! AC%value(k)= - DBLE(dOverd0*(two + cosLat*(cosLatMHalf+cosLatPHalf)) + (oneOverLd02*R2*cosLat2*dLat2)) 
				
				BC%j(k)=i; BC%value(k)=-DBLE(two + cosLat*(cosLatMHalf+cosLatPHalf))

				!COMMENT: ψ₁(ƛᵢ₊₁, θᵢ) 
				k=k+1;
				AC%j(k)=iP1; AC%value(k)=dOverd0   !DBLE(oneOverCosLat2) !1/sin^2(theta_{ij})
				BC%j(k)=iP1; BC%value(k)=one
				!COMMENT: ψ₁(ƛᵢ,θᵢ₊₁)
				if(iLat.lt.sphericalGrid%nLat)then 
					k=k+1;
					AC%j(k)=i+sphericalGrid%nPhi;AC%value(k)=dOverd0*cosLat*cosLatPHalf 
					BC%j(k)=i+sphericalGrid%nPhi;BC%value(k)=cosLat*cosLatPHalf
				end if

				!ALTERED:UNSCALED !!
				! !COMMENT: ψ₁(ƛᵢ,θᵢ₋₁)
				! if(iLat.gt.1)then  
				! 	k=k+1;
				! 	AC%j(k)=i-sphericalGrid%nPhi;AC%value(k)=dOverd0*cosLatMHalf*(one/ (cosLat*dLat2*R2)) 
				! 	BC%j(k)=i-sphericalGrid%nPhi;BC%value(k)=cosLat*cosLatMHalf 
				! end if
				! !COMMENT: ψ₁(ƛᵢ₋₁, θᵢ)
				! k=k+1;
				! AC%j(k)=iM1; AC%value(k)= dOverd0 *(one/(cosLat2*dLat2*R2)) 
				! BC%j(k)=iM1; BC%value(k)=one

				! !COMMENT: ψ₁(ƛᵢ,θᵢ)
				! k=k+1; 
				! AC%j(k)=i; 
				! AC%value(k)=-DBLE(dOverd0* (two *(one/(cosLat2*dLat2*R2)) + (one/(cosLat*dLat2*R2))*(cosLatMHalf+cosLatPHalf)) + (oneOverLd02) ) 
				! BC%j(k)=i; BC%value(k)=-DBLE(two + cosLat*(cosLatMHalf+cosLatPHalf))

				! !COMMENT: ψ₁(ƛᵢ₊₁, θᵢ) 
				! k=k+1;
				! AC%j(k)=iP1; AC%value(k)=dOverd0*(one/(cosLat2*dLat2*R2))  
				! BC%j(k)=iP1; BC%value(k)=one
				! !COMMENT: ψ₁(ƛᵢ,θᵢ₊₁)
				! if(iLat.lt.sphericalGrid%nLat)then 
				! 	k=k+1;
				! 	AC%j(k)=i+sphericalGrid%nPhi;AC%value(k)=dOverd0*cosLatPHalf*(one/(cosLat*dLat2*R2)) 
				! 	BC%j(k)=i+sphericalGrid%nPhi;BC%value(k)=cosLat*cosLatPHalf
				! end if

				AC%nentries(i) = k
				BC%nentries(i) = k
			end do
		end do

		AC%m = nElements
		AC%n = nElements
		AC%rowIndexed = 1

		BC%m = nElements
		BC%n = nElements
		BC%rowIndexed = 1
		write(*,*)'coefficient matrix is built=',omp_get_wtime()-time1
		return

	end subroutine buildCoefficientsMatrix

	subroutine solvePoissonForPsi()
		use sphericalGridDataModule
		use sparseMatricesDataModule
		use basicDataStructures
		use simulationParameters
		use lsqrModule
		use GMRESSettings
		use bathymetryDataModule
		!$ use omp_lib
		implicit none
		! here, we solve Laplacian_s(Psi) - f(theta,Lon) Psi/Ld^2 = RHS
		! the coefficient matrix of the operator Laplacian_s(Psi) - f(theta,Lon) Psi/Ld^2
		! in spherical coordinates
		! R^{-2} [1/sin(theta) diff(sin(theta) diff(psi,theta),theta) + 1/sin^2(theta) diff(diff(psi,Lon),Lon)] - Ld^{-2} f(theta,Lon) psi = RHS0
		!
		! [1/sin(theta) diff(sin(theta) diff(psi,theta),theta) + 1/sin^2(theta) diff(diff(psi,Lon),Lon)] - (R/Ld0)^2 f(theta,Lon) psi = R^2 RHS0
		!
		! 1/sin(theta_{ij}) [sin(theta_{i,j+1/2}) [psi(i,j+1) - psi(i,j)]/dLat - sin(theta_{i,j-1/2}) [psi(i,j) - psi(i,j-1)]/dLat]/dLat
		! + 1/sin^2(theta_{ij}) [psi(i+1,j) - 2 psi(i,j) + psi(i-1,j)]/DLon^2 - (R/Ld0)^2 f(i,j) psi(i,j) = R^2 RHS0
		!
		! coefficients
		! psi(i,j) : -1/sin(theta_{ij}) [sin(theta_{i,j+1/2}) + sin(theta_{i,j-1/2}]/dLat^2 - 2/sin^2(theta_{ij})/DLon^2 - (R/Ld0)^2 f(i,j)
		!            -2/dLat^2 - 2/sin^2(theta_{ij})/DLon^2 - (R/Ld0)^2 f(i,j)
		! psi(i-1,j): 1/sin^2(theta_{ij}) 1/DLon^2
		! psi(i+1,j): 1/sin^2(theta_{ij}) 1/DLon^2
		! psi(i,j-1): sin(theta_{i,j-1/2})/sin(theta_{ij}) 1/dLat^2
		! psi(i,j+1): sin(theta_{i,j+1/2})/sin(theta_{ij}) 1/dLat^2
		!
		! if dLat=dLon, the new right hand side is dLat^2 R^2 RHS0
		!
		! psi(i,j) : -2 - 2/sin^2(theta_{ij}) - (R/Ld0)^2 dLat^2 f(i,j)
		! psi(i-1,j): 1/sin^2(theta_{ij})
		! psi(i+1,j): 1/sin^2(theta_{ij})
		! psi(i,j-1): sin(theta_{i,j-1/2})/sin(theta_{ij})
		! psi(i,j+1): sin(theta_{i,j+1/2})/sin(theta_{ij})
		! RHS = dLat^2 R^2 RHS0
		!
		! If using the latitude
		!
		! R^{-2} [1/sin(theta) diff(sin(theta) diff(psi,theta),theta) + 1/sin^2(theta) diff(diff(psi,Lon),Lon)] - Ld^{-2} f(Lon,theta) psi = RHS0
		! R^{-2} [1/cos(lat) diff(cos(lat) diff(psi,lat),lat) + 1/cos^2(theta) diff(diff(psi,Lon),Lon)] - Ld^{-2} f(Lon,theta) psi = RHS0
		!
		! psi(i,j) : -2 - 2/cos^2(lat_{ij}) - (R/Ld0)^2 dLat^2 f(ij)
		! psi(i-1,j): 1/cos^2(lat_{ij})
		! psi(i+1,j): 1/cos^2(lat_{ij})
		! psi(i,j-1): cos(lat_{i,j-1/2})/cos(lat_{ij}) = 1 + (dLat/2) tan(lat_{ij})
		! psi(i,j+1): cos(lat_{i,j+1/2})/cos(lat_{ij}) = 1 - (dLat/2) tan(lat_{ij})
		! RHS = dLat^2 R^2 RHS0
		!
		! Multiplying everybody by cos^2(lat_{ij})
		!
		! psi(i,j) : -2 cos^2(lat_{ij}) - 2 - (R/Ld0)^2 dLat^2 f(ij) cos^2(lat_{ij})
		! psi(i-1,j): 1
		! psi(i+1,j): 1
		! psi(i,j-1): cos(lat_{ij}) cos(lat_{i,j-1/2})
		! psi(i,j+1): cos(lat_{ij}) cos(lat_{i,j+1/2})
		! RHS = dLat^2 R^2 cos^2(lat_{ij})	RHS0
		!

		integer :: i,k,iLat,iLon,iM1,iP1,ki,m,nElements
		real (kind=QR_K) :: ROverLd02, dLat, dLat2,lat,cosLat,sinLat,tanLat,Lon
		real (kind=QR_K) :: oneOverCosLat2,a,R2,aOverR2,aOverLd02,factor,psiOvera,Qminusf
		real (kind=QR_K) :: cosLatMHalf,cosLatPHalf,cosLat2,aOverR,sinLat2,fij
		double precision, allocatable :: b(:),res(:),sol(:)
		logical :: signConstrained
		integer :: i1,i2,CHUNK,nThreads,threadID
		integer :: iLatP1,iLatM1,iLonP1,iLonM1 , endIdx, startIdx,row
		double precision :: time1
		real (kind=QR_K) :: d2, lat0, Lon0, sigma02, coriolisTerm,dOverd0

		time1=omp_get_wtime()

		write(*,*)'solving the Poisson equation',sphericalGrid%nLat,sphericalGrid%nPhi
		nElements=sphericalGrid%nLat*sphericalGrid%nPhi
		allocate(b(nElements))
		b = zero

		i=0
		R2=R*R
		dLat=sphericalGrid%dLat
		dLat2=dLat*dLat
		do iLat=1,sphericalGrid%nLat
			lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
			cosLat=sphericalGrid%cosLatC(iLat);
			cosLat2=cosLat*cosLat
			sinLat=sphericalGrid%sinLatC(iLat)
			coriolisTerm = two*OMEGA*sinLat !coriolisTerm has the units of /sec

			do iLon=1,sphericalGrid%nPhi
				Lon=half*sphericalGrid%dPhi+sphericalGrid%dPhi*real(iLon-1,QR_K)
				i=i+1
				if (sphericalGrid%d(iLat,iLon).gt.zero)then
					!ADDED: To double check the fraction d/d0 :
					if (iLat == 1 .and. iLon == 1) then
						write(*,*) 'sphericalGrid%d(iLat,iLon)/sphericalGrid%d0', &
									sphericalGrid%d(iLat,iLon) / sphericalGrid%d0
					end if

					dOverd0 = sphericalGrid%d(iLat,iLon)/sphericalGrid%d0
				
					Qminusf= - coriolisTerm! QOnGrid(iLat,iLon)*dOverd0 ;

					!ALTERED: SCALED 
					b(i) = DBLE(dOverd0*cosLat2*dLat2*R2*Qminusf) 
					
					!ALTERED: UNSCALED 
					! b(i) = DBLE(dOverd0*Qminusf) 
					
				else
					b(i) = zero
				end if
				if(ISNAN(b(i)))stop "b nan"
			end do
		end do
		
		allocate(sol(nElements))
		allocate(res(nElements))
		sol = zero
		res = zero
		
		signConstrained=.FALSE.
		!		call GMRES_RESTARTED_OMP_RETURNR(AC,b(1:nElements),sol(1:nElements),signConstrained,&
		!		&res(1:nElements),nElements,GMRES_ninner,GMRES_nouter,GMRES_tol,successFlag,normr)

		!BCG(A,b,x,m,n,nIter,res)
		!write(*,*)'calling BICGSTAB_OMP'

		call BICGSTAB_OMP(AC,b(1:nElements),sol(1:nElements),nElements,GMRES_nouter,GMRES_tol)

		if(.NOT.allocated(psiOnGrid))allocate(psiOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		psiOnGrid=zero

		i=0
		do iLat=1,sphericalGrid%nLat
			do iLon=1,sphericalGrid%nPhi
				i=i+1
				psiOnGrid(iLat,iLon)=sol(i)
			end do
		end do

		!!Psi
		! open(unit=14, file='../Python/output/petsc/scaled/PsiSpherical.csv', status='replace', action='write', form='formatted')
		! write(14,'(A)') 'iLat,iLon,psi'
		! do iLat = 1, sphericalGrid%nLat
		! 	do iLon = 1, sphericalGrid%nPhi
		! 		write(14,'(I0,",",I0,",",ES24.16)') iLat, iLon, psiOnGrid(iLat,iLon)
		! 	end do
		! end do
		! close(14)
		! write(*,*) 'Save PsiSpherical'
		!!b
		! open(unit=15, file='../Python/output/petsc/scaled/b.csv', status='replace', action='write', form='formatted')
		! write(15,'(A)') 'iLat,iLon,b'
		! i = 0
		! do iLat = 1, sphericalGrid%nLat
		! 	do iLon = 1, sphericalGrid%nPhi
		! 		i = i + 1
		! 		write(15,'(I0,",",I0,",",ES24.16)') iLat, iLon, b(i)
		! 	end do
		! end do
		! close(15)
		! write(*,*) 'Saved b.csv'

		!!AC
		! open(unit=16, file='../Python/output/petsc/scaled/AC.csv', status='replace', action='write', form='formatted')
		! write(16,'(A)') 'row,col,value'  ! exact header
		! do row = 1, AC%m
		! 	if (row == 1) then
		! 		startIdx = 1
		! 	else
		! 		startIdx = AC%nentries(row-1) + 1
		! 	end if
		! 	endIdx = AC%nentries(row)
		! 	do k = startIdx, endIdx
		! 		! Write exactly with commas
		! 		write(16,'(I0,",",I0,",",ES24.16)') row, AC%j(k), AC%value(k)
		! 	end do
		! end do
		! close(16)
		! write(*,*) 'Saved AC.csv'

		write(*,*) 'OLD: ψ max, ψ min = ', maxval(psiOnGrid), minval(psiOnGrid)

		deallocate(b)
		deallocate(sol)
		deallocate(res)

		!computing velocity using finite difference of psi
		if(.NOT.allocated(velocityVorOnGrid))allocate(velocityVorOnGrid(2,sphericalGrid%nLat,sphericalGrid%nPhi))
		velocityVorOnGrid=zero
		if(.NOT.allocated(velocityOnGrid))allocate(velocityOnGrid(2,sphericalGrid%nLat,sphericalGrid%nPhi))
		velocityOnGrid=zero

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = sphericalGrid%nLat/nThreads
		write(*,*)CHUNK

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,sphericalGrid,&
		!$OMP &velocityVorOnGrid,velocityOnGrid,psiOnGrid)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		!		write(*,*)threadID,i1,i2
		if(threadID+1.eq.nThreads)i2=sphericalGrid%nLat
		!		do iLat=1,sphericalGrid%nLat
		do iLat=i1,i2
			cosLat=sphericalGrid%cosLatC(iLat);
			sinLat=sphericalGrid%sinLatC(iLat);
			iLatP1=iLat+1; iLatM1=iLat-1
			if(iLat.eq.1)iLatM1=1
			if(iLat.eq.sphericalGrid%nLat)iLatP1=sphericalGrid%nLat
			
			do iLon=1,sphericalGrid%nPhi
				Lon=half*sphericalGrid%dPhi+sphericalGrid%dPhi*real(iLon-1,QR_K)
				iLonP1=iLon+1;iLonM1=iLon-1
				if(iLon.eq.1)iLonM1=sphericalGrid%nPhi
				if(iLon.eq.sphericalGrid%nPhi)iLonP1=1

				velocityVorOnGrid(1,iLat,iLon)=(psiOnGrid(iLat,iLonP1)-psiOnGrid(iLat,iLonM1))&
				&/(two*sphericalGrid%dLat*R*cosLat)
				velocityVorOnGrid(2,iLat,iLon)=-(psiOnGrid(iLatP1,iLon)-psiOnGrid(iLatM1,iLon))&
				&/(two*sphericalGrid%dLat*R)
				velocityOnGrid(1,iLat,iLon)=velocityVorOnGrid(1,iLat,iLon)
				velocityOnGrid(2,iLat,iLon)=velocityVorOnGrid(2,iLat,iLon)
			end do
		end do
		!$OMP END PARALLEL

		write(*,*)'POISSON SOLVER sim time=',omp_get_wtime()-time1

		write(*,*)'u theta  max , min = ',maxval(velocityOnGrid(1,:,:)) , minval(velocityOnGrid(1,:,:))
		write(*,*)'u lambda max , min = ',maxval(velocityOnGrid(2,:,:)) , minval(velocityOnGrid(2,:,:)) 

		write(*,*) 'Velocity by the old subroutine Vmax, Vmin = ', maxval(sqrt(velocityOnGrid(1,:,:)**2+velocityOnGrid(2,:,:)**2)), minval(sqrt(velocityOnGrid(1,:,:)**2 + velocityOnGrid(2,:,:)**2))

		return

	end subroutine solvePoissonForPsi

	subroutine solvePoissonForPhi()
		use sphericalGridDataModule
		use sparseMatricesDataModule
		use basicDataStructures
		use simulationParameters
		use lsqrModule
		use GMRESSettings
		use bathymetryDataModule
		!$ use omp_lib
		implicit none
		! here, we solve Laplacian_s(Phi) = RHS
		! the coefficient matrix of the operator Laplacian_s(Phi)
		! in spherical coordinates
		! R^{-2} [1/sin(theta) diff(sin(theta) diff(Phi,theta),theta) + 1/sin^2(theta) diff(diff(Phi,Lon),Lon)] = RHS0
		!
		! [1/sin(theta) diff(sin(theta) diff(Phi,theta),theta) + 1/sin^2(theta) diff(diff(Phi,Lon),Lon)]  = R^2 RHS0
		!
		! 1/sin(theta_{ij}) [sin(theta_{i,j+1/2}) [Phi(i,j+1) - Phi(i,j)]/dLat - sin(theta_{i,j-1/2}) [Phi(i,j) - Phi(i,j-1)]/dLat]/dLat
		! + 1/sin^2(theta_{ij}) [Phi(i+1,j) - 2 Phi(i,j) + Phi(i-1,j)]/DLon^2 = R^2 RHS0
		!
		! coefficients
		! Phi(i,j) : -1/sin(theta_{ij}) [sin(theta_{i,j+1/2}) + sin(theta_{i,j-1/2}]/dLat^2 - 2/sin^2(theta_{ij})/DLon^2
		!            -2/dLat^2 - 2/sin^2(theta_{ij})/DLon^2
		! Phi(i-1,j): 1/sin^2(theta_{ij}) 1/DLon^2
		! Phi(i+1,j): 1/sin^2(theta_{ij}) 1/DLon^2
		! Phi(i,j-1): sin(theta_{i,j-1/2})/sin(theta_{ij}) 1/dLat^2
		! Phi(i,j+1): sin(theta_{i,j+1/2})/sin(theta_{ij}) 1/dLat^2
		!
		! if dLat=dLon, the new right hand side is dLat^2 R^2 RHS0
		!
		! Phi(i,j) : -2 - 2/sin^2(theta_{ij})
		! Phi(i-1,j): 1/sin^2(theta_{ij})
		! Phi(i+1,j): 1/sin^2(theta_{ij})
		! Phi(i,j-1): sin(theta_{i,j-1/2})/sin(theta_{ij})
		! Phi(i,j+1): sin(theta_{i,j+1/2})/sin(theta_{ij})
		! RHS = dLat^2 R^2 RHS0
		!
		! If using the latitude
		!
		! R^{-2} [1/sin(theta) diff(sin(theta) diff(Phi,theta),theta) + 1/sin^2(theta) diff(diff(Phi,Lon),Lon)]  Phi = RHS0
		! R^{-2} [1/cos(lat) diff(cos(lat) diff(Phi,lat),lat) + 1/cos^2(theta) diff(diff(Phi,Lon),Lon)]  Phi = RHS0
		!
		! Phi(i,j) : -2 - 2/cos^2(lat_{ij})
		! Phi(i-1,j): 1/cos^2(lat_{ij})
		! Phi(i+1,j): 1/cos^2(lat_{ij})
		! Phi(i,j-1): cos(lat_{i,j-1/2})/cos(lat_{ij}) = 1 + (dLat/2) tan(lat_{ij})
		! Phi(i,j+1): cos(lat_{i,j+1/2})/cos(lat_{ij}) = 1 - (dLat/2) tan(lat_{ij})
		! RHS = dLat^2 R^2 RHS0
		!
		! Multiplying everybody by cos^2(lat_{ij})
		!
		! Phi(i,j) : -2 cos^2(lat_{ij}) - 2
		! Phi(i-1,j): 1
		! Phi(i+1,j): 1
		! Phi(i,j-1): cos(lat_{ij}) cos(lat_{i,j-1/2})
		! Phi(i,j+1): cos(lat_{ij}) cos(lat_{i,j+1/2})
		! RHS = dLat^2 R^2 cos^2(lat_{ij})	RHS0
		!

		integer :: i,k,iLat,iLon,iM1,iP1,ki,m,nElements
		real (kind=QR_K) :: ROverLd02, dLat, dLat2,lat,cosLat,sinLat,tanLat,Lon
		real (kind=QR_K) :: oneOverCosLat2,a,R2,aOverR2,aOverLd02,factor,psiOvera,Qminusf
		real (kind=QR_K) :: cosLatMHalf,cosLatPHalf,cosLat2,aOverR,sinLat2,fij
		double precision, allocatable :: b(:),res(:),sol(:)
		logical :: signConstrained
		integer :: i1,i2,CHUNK,nThreads,threadID
		integer :: iLatP1,iLatM1,iLonP1,iLonM1
		double precision :: time1
		real (kind=QR_K) :: d2, lat0, Lon0, sigma02, coriolisTerm,dOverd0

		time1=omp_get_wtime()

		write(*,*)'solving the Poisson equation',sphericalGrid%nLat,sphericalGrid%nPhi
		nElements=sphericalGrid%nLat*sphericalGrid%nPhi
		allocate(b(nElements))
		b = zero

		i=0
		R2=R*R
		dLat=sphericalGrid%dLat
		dLat2=dLat*dLat
		do iLat=1,sphericalGrid%nLat
			lat=-half*pi+half*sphericalGrid%dLat+sphericalGrid%dLat*real(iLat-1,QR_K)
			cosLat=sphericalGrid%cosLatC(iLat);
			cosLat2=cosLat*cosLat
			sinLat=sphericalGrid%sinLatC(iLat)
			coriolisTerm = two*OMEGA*sinLat !coriolisTerm has the units of /sec

			do iLon=1,sphericalGrid%nPhi
				Lon=half*sphericalGrid%dPhi+sphericalGrid%dPhi*real(iLon-1,QR_K)
				i=i+1
				if (sphericalGrid%d(iLat,iLon).gt.zero)then
					b(i) = DBLE(divOnGrid(iLat,iLon))
				else
					b(i) = zero
				end if			
				if(ISNAN(b(i)))stop "b nan"
			end do
		end do

		allocate(sol(nElements))
		allocate(res(nElements))
		sol = zero
		res = zero

		signConstrained=.FALSE.

		call BICGSTAB_OMP(BC,b(1:nElements),sol(1:nElements),nElements,GMRES_nouter,GMRES_tol)

		if(.NOT.allocated(phiOnGrid))allocate(phiOnGrid(sphericalGrid%nLat,sphericalGrid%nPhi))
		phiOnGrid=zero

		i=0
		do iLat=1,sphericalGrid%nLat
			do iLon=1,sphericalGrid%nPhi
				i=i+1
				phiOnGrid(iLat,iLon)=sol(i)
			end do
		end do

		deallocate(b)
		deallocate(sol)
		deallocate(res)

		!computing velocity using finite difference of phi
		if(.NOT.allocated(velocityDivOnGrid))allocate(velocityDivOnGrid(2,sphericalGrid%nLat,sphericalGrid%nPhi))
		velocityDivOnGrid=zero
		if(.NOT.allocated(velocityOnGrid))allocate(velocityOnGrid(2,sphericalGrid%nLat,sphericalGrid%nPhi))
		velocityOnGrid=zero

		nThreads = OMP_get_max_threads()
		call OMP_SET_NUM_THREADS(nThreads)
		CHUNK = sphericalGrid%nLat/nThreads
		write(*,*)CHUNK

		!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CHUNK,sphericalGrid,&
		!$OMP &velocityDivOnGrid,velocityOnGrid,velocityVorOnGrid,phiOnGrid)
		threadID = omp_get_thread_num()
		i1 = threadID*CHUNK+1
		i2=i1-1+CHUNK
		!		write(*,*)threadID,i1,i2
		if(threadID+1.eq.nThreads)i2=sphericalGrid%nLat
		!		do iLat=1,sphericalGrid%nLat
		do iLat=i1,i2
			cosLat=sphericalGrid%cosLatC(iLat);
			sinLat=sphericalGrid%sinLatC(iLat);
			iLatP1=iLat+1; iLatM1=iLat-1
			if(iLat.eq.1)iLatM1=1
			if(iLat.eq.sphericalGrid%nLat)iLatP1=sphericalGrid%nLat
			do iLon=1,sphericalGrid%nPhi
				Lon=half*sphericalGrid%dPhi+sphericalGrid%dPhi*real(iLon-1,QR_K)
				iLonP1=iLon+1;iLonM1=iLon-1
				if(iLon.eq.1)iLonM1=sphericalGrid%nPhi
				if(iLon.eq.sphericalGrid%nPhi)iLonP1=1

				velocityDivOnGrid(1,iLat,iLon)=(phiOnGrid(iLatP1,iLon)-phiOnGrid(iLatM1,iLon))&
				&/(two*sphericalGrid%dLat*R)
				velocityDivOnGrid(2,iLat,iLon)=(phiOnGrid(iLat,iLonP1)-phiOnGrid(iLat,iLonM1))&
				&/(two*sphericalGrid%dLat*R*cosLat)

				velocityOnGrid(1,iLat,iLon)=velocityDivOnGrid(1,iLat,iLon)+velocityVorOnGrid(1,iLat,iLon)
				velocityOnGrid(2,iLat,iLon)=velocityDivOnGrid(2,iLat,iLon)+velocityVorOnGrid(2,iLat,iLon)
			end do
		end do
		!$OMP END PARALLEL

		write(*,*)'POISSON SOLVER for phi sim time=',omp_get_wtime()-time1
		write(*,*) 'MAX PHI , MIN PHI = ', maxval(phiOnGrid) , minval(phiOnGrid)
		write(*,*) 'MAX VELOCITY  TOTAL , MIN VELOCITY TOTAL = ', maxval(sqrt(velocityOnGrid(1,:,:)**2+ velocityOnGrid(2,:,:)**2) ), minval(sqrt(velocityOnGrid(1,:,:)**2+ velocityOnGrid(2,:,:)**2) )

		return

	end subroutine solvePoissonForPhi

	subroutine removeCoefficientsMatrix
		use sparseMatricesDataModule
		implicit none
		deallocate(AC%value)
		deallocate(AC%j)
		deallocate(AC%nentries)

		deallocate(BC%value)
		deallocate(BC%j)
		deallocate(BC%nentries)
		return
	end subroutine removeCoefficientsMatrix

!ADDED SUBROUTINE: Building Coefficient Matrix for ψ₁ and appending the values to OceanGrid. Discretization is shown below.                                         
	subroutine OceanbuildCoefficientsMatrixPsiOne()
		use OceanGridDataModule
		use sparseMatricesDataModule
		use basicDataStructures
		use simulationParameters
		use lsqrModule
		use GMRESSettings
		use bathymetryDataModule
		!$ use omp_lib
		implicit none

		integer :: i,k,iLat,iLon,iM1,iP1,nElements 
		real (kind=QR_K) :: dLat, dLat2,lat,cosLat,sinLat,Lon, R2, ScreeningTerm,Ld, Ld2 ,beta,gamma,alpha
		real (kind=QR_K) :: latMHalf,latPHalf ,cosLatMHalf,cosLatPHalf,cosLat2,sinLat2, sinLatM1, sinLatP1
		double precision, allocatable :: b(:),res(:),sol(:)
		real (kind=QR_K) :: coriolisTerm,coriolisTerm2
		logical :: signConstrained
		integer :: i1,i2,CHUNK,nThreads,threadID
		integer :: iLatP1,iLatM1,iLonP1,iLonM1 ,epsilon
		double precision :: time1
		integer :: row, startIdx, endIdx , k00,k11 , halfLat
		real (kind=QR_K) :: dijOverd0 , theta0, OMEGA2,CORIOLIS,CORIOLIS2, d
		logical :: unscaled_printed , scaled_printed

		time1=omp_get_wtime()
		unscaled_printed = .FALSE.
		scaled_printed = .FALSE.
	
		write(*,*)' Dimensions of the Coefficient Matrix = ',OceanGrid%nLat,OceanGrid%nPhi 
		AC%length = (OceanGrid%nLat-2) * OceanGrid%nPhi*5 + 2*OceanGrid%nPhi*4
	
		write(*,*) ' Length Assigned for the Coefficient Matrix = ', (OceanGrid%nLat-2) * OceanGrid%nPhi*5 + 2*OceanGrid%nPhi*4
	
		if(allocated(AC%j))deallocate(AC%j)
		if(allocated(AC%value))deallocate(AC%value)
		if(allocated(AC%nentries))deallocate(AC%nentries)
		allocate(AC%j(AC%length))
		allocate(AC%value(AC%length))
		nElements=OceanGrid%nLat*OceanGrid%nPhi
		allocate(AC%nentries(nElements))
		AC%value=zero
	
		!!initializing the indices
		i=0; k=0
	
		dLat=OceanGrid%dLat
		dLat2=dLat*dLat
		R2= R*R
	
		!ALTERED: ∇²ψ = 1/(R² cos²θ ) ∂²ψ/∂ƛ² + 1/(R² cosθ) ∂/∂θ (cosθ ∂ψ/∂θ)
		theta0 = pi/four
		OMEGA2 = OMEGA*OMEGA
		CORIOLIS=two*OMEGA*sin(theta0)

		! ALTERED: DOUBLE CHECK IF YOU WANT CONSTANT OR VARYING CORIOLIS
		CORIOLIS2= CORIOLIS*CORIOLIS

		write(*,*)'BUILDING THE COEFFICIENT MATRIX FOR ψ₁'

		do iLat=1,OceanGrid%nLat
			lat = -half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)
			
			cosLat=OceanGrid%cosLatC(iLat);
			cosLat2=cosLat*cosLat
			
			!COMMENT: DEFINING HALF COSINES AT THE N-S POLES
			latPHalf = Lat+half*dLat; cosLatPHalf=cos(latPHalf)  
			latMHalf = Lat -half*dLat;cosLatMHalf=cos(latMHalf)
			if(iLat.eq.1) cosLatMHalf=zero
			if(iLat.eq.OceanGrid%nLat) cosLatPHalf=zero
			
			sinLat=OceanGrid%sinLatC(iLat)
	
			coriolisTerm = two*OMEGA*sinLat
			coriolisTerm2= coriolisTerm*coriolisTerm 
			!ALTERED: DOUBLE CHECK IF YOU WANT CONSTANT OR VARYING CORIOLIS
			! CORIOLIS2 = coriolisTerm2
			do iLon=1,OceanGrid%nPhi
				i=i+1
				Lon=half*OceanGrid%dPhi+OceanGrid%dPhi*real(iLon-1,QR_K)
			
				!COMMENT: flag off depth
				if (OceanGrid%d(iLat,iLon).gt.zero)then
					d = OceanGrid%d(iLat,iLon)
				else
					d = zero
				end if

				!COMMENT: Periodic Boundary Conditions for ƛ
				iM1=i-1; if(iLon.eq.1) iM1=i-1+OceanGrid%nPhi
				iP1=i+1; if(iLon.eq.OceanGrid%nPhi) iP1= i+1-OceanGrid%nPhi

				!COMMENT:ψ₁(ƛᵢ,θᵢ)(center)
				k=k+1;
				! k00 = k
				AC%j(k)=i
				AC%value(k)= d*(-2-(cosLatMHalf + cosLatPHalf)*cosLat) - ((coriolisTerm2)/GRAVITY)*(R2*dLat2*cosLat2)
				
				! -(one/Ld0**2) * (R2*dLat2*cosLat)

				!COMMENT:(left) ψ₁(ƛᵢ₋₁, θᵢ)
				k=k+1;
				AC%j(k)=iM1
				AC%value(k)= d*one 
				!COMMENT:(right) ψ₁(ƛᵢ₊₁, θᵢ)   
				k=k+1; 
				AC%j(k)=iP1
				AC%value(k)= d*one 

				!COMMENT:(bottom) ψ₁(ƛᵢ,θᵢ₋₁) 
				if(iLat.gt.1)then  
					k=k+1;
					AC%j(k)=i - OceanGrid%nPhi
					AC%value(k)= d*cosLatMHalf *cosLat
				end if

				!COMMENT:(top) ψ₁(ƛᵢ,θᵢ₊₁) 
				if(iLat.lt.OceanGrid%nLat)then 
					k=k+1;
					AC%j(k)=i + OceanGrid%nPhi
					AC%value(k)= d*cosLatPHalf *cosLat
				end if

				AC%nentries(i) = k
			end do
		end do
	
		AC%m = nElements
		AC%n = nElements
		AC%rowIndexed = 1
!		write(*,*)AC%m, AC%length, AC%nentries(nElements), k
!		stop

		return 
	end subroutine OceanbuildCoefficientsMatrixPsiOne

!ADDED SUBROUTINE: Solving for ψ₁.
subroutine solvePoissonForPsiOne()
	use OceanGridDataModule
	use OceanAdaptiveGridDataModule 
	use sparseMatricesDataModule
	use basicDataStructures
	use simulationParameters
	use lsqrModule
	use GMRESSettings
	use bathymetryDataModule
	use AnalysisModule !ADDED!
	! use petscModule !not working yet idk if i will use it later
	!$ use omp_lib
	implicit none

	integer :: i,k,iLat,iLon,iM1,iP1,ki,m,nElements
	real (kind=QR_K) :: dLat, dLat2,lat,cosLat,sinLat,Lon, R2 , dphi ,dphi2
	real (kind=QR_K) :: cosLatMHalf,cosLatPHalf,cosLat2,sinLat2 , latPHalf, latMHalf
	double precision, allocatable :: b(:),res(:),sol(:)
	logical :: signConstrained
	integer :: i1,i2,CHUNK,nThreads,threadID , rasberries , row,endIdx,startIdx
	integer :: iLatP1,iLatM1,iLonP1,iLonM1,CorrespondingiLonM1,CorrespondingiLonP1, epsilon, index_equatorial , equatorial_particles, northern_hemisphere_particles, southern_hemisphere_particles
	double precision :: time1 , tol_BlueBerries
	real (kind=QR_K) :: d2, coriolisTerm,dijOverd0 ,d
	
	integer :: successFlag !ADDED 
	double precision :: normr !ADDED
	double precision :: HighestTerm, tmp

	integer :: nRows, nnz, ierrOut
	integer, allocatable :: rowPtr_csr(:), colInd_csr(:) 
	double precision, allocatable :: val_csr(:) 
	double precision, allocatable :: val_coo_dp(:)
	
	real (kind=QR_K) :: counter,alpha,beta, psiave , psimax

	time1=omp_get_wtime()

	write(*,*)'SOLVING THE POISSON EQUATION FOR ψ₁' ,OceanGrid%nLat,OceanGrid%nPhi
	
	nElements= OceanGrid%nLat*OceanGrid%nPhi
	allocate(b(nElements))
	b = zero
	
	i=0 
	R2=R*R
	dLat=OceanGrid%dLat
	dLat2=dLat*dLat
	dphi=OceanGrid%dphi
	dphi2=dphi*dphi

	do iLat=1,OceanGrid%nLat
		lat=-half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)

		cosLat=OceanGrid%cosLatC(iLat);
		cosLat2=cosLat*cosLat
		
		sinLat=OceanGrid%sinLatC(iLat)
		coriolisTerm = two*OMEGA*sinLat

		do iLon=1,OceanGrid%nPhi
			Lon=half*OceanGrid%dPhi+OceanGrid%dPhi*real(iLon-1,QR_K)
			i=i+1
			if (OceanGrid%d(iLat,iLon).gt.zero)then
				d = OceanGrid%d(iLat,iLon)
			else
				d=zero 
			end if

			b(i) = - d *coriolisTerm*(dLat2*R2*cosLat2)

			if(ISNAN(b(i)))stop "b nan"
		end do
	end do

	allocate(sol(nElements))
	allocate(res(nElements))
	sol = zero
	res = zero

	!BUGTRACKING Stabilization                                       
	HighestTerm = zero

	do i = 1, AC%m
		if (i == 1) then
			startIdx = 1
		else
			startIdx = AC%nentries(i-1) + 1
		end if
		endIdx = AC%nentries(i)

		do k = startIdx, endIdx
			if (AC%j(k) == i) then
				HighestTerm = max(HighestTerm, abs(AC%value(k)))
				exit   ! only one diagonal per row
			end if
		end do
	end do

	if (HighestTerm == zero) then
		write(*,*) ' HighestTerm == zero, take it to be one '
		HighestTerm = one
	end if

	! HighestTerm=one

	!COMMENT: divide matrix values by HighestTerm
	do i = 1, AC%nentries(nElements)   !AC%length, 
		AC%value(i) = AC%value(i) / HighestTerm
	end do

	!COMMENT: divide RHS vector
	do i = 1, nElements
		b(i) = b(i) / HighestTerm
	end do

	! signConstrained=.FALSE.
	!		call GMRES_RESTARTED_OMP_RETURNR(AC,b(1:nElements),sol(1:nElements),signConstrained,&
	!		&res(1:nElements),nElements,GMRES_ninner,GMRES_nouter,GMRES_tol,successFlag,normr)

	!BCG(A,b,x,m,n,nIter,res)
	!write(*,*)'calling BICGSTAB_OMP'

	!COMMENT: input:     BICGSTAB_OMP(A,b,x,n,nIter,res) 
	!COMMENT: BICGSTAB target residual=',k,nIter,resnorm,res*normb
	call BICGSTAB_OMP(AC,b(1:nElements),sol(1:nElements),nElements,GMRES_nouter,GMRES_tol)
	write(*,*)'BICGSTAB_OMP returned successfully'

	!COMMENT: Solution of the streamfunction ψ₁
	if(.NOT.allocated(psi1onGrid))allocate(psi1onGrid(OceanGrid%nLat,OceanGrid%nPhi))
	psi1onGrid=zero
	i=0
	do iLat=1,OceanGrid%nLat
		do iLon=1,OceanGrid%nPhi
			i=i+1
			psi1onGrid(iLat,iLon)= sol(i)
		end do
	end do

	! write(*,*)'psi=',psi1onGrid(1,1:OceanGrid%nPhi)
	write(*,*)'POISSON SOLVER for ψ₁ sim time =', omp_get_wtime()-time1
	write(*,*)' '

	write(*,*) 'ψ₁: Max , Min = ', maxval(psi1onGrid), minval(psi1onGrid)

	!ADDED: ______________________________________________________________________________________________________________
	call system("mkdir -p ../Python/output/Fortran/EquationOne/Normalized")
	
	!COMMENT: save the entries for comparison: Save AC, b and sol
	! ! !!RHS
	! open(unit=10, file='../Python/output/Fortran/EquationOne/Normalized/b.csv', status='replace', action='write', form='formatted')
	! write(10,'(A)') 'iLat,iLon,b'
	! i = 0
	! do iLat = 1, OceanGrid%nLat
	! 	do iLon = 1, OceanGrid%nPhi
	! 		i = i + 1
	! 		write(10,'(I0,",",I0,",",ES24.16)') iLat, iLon, b(i)
	! 	end do
	! end do
	! close(10)
	! write(*,*) 'Saved b.csv'

	!!SOLUTION PSI
	open(unit=11, file='../Python/output/Fortran/EquationOne/Normalized/PsiOcean.csv', status='replace', action='write', form='formatted')
	write(11,'(A)') 'iLat,iLon,psi'
	do iLat = 1, OceanGrid%nLat
		do iLon = 1, OceanGrid%nPhi
			write(11,'(I0,",",I0,",",ES24.16)') iLat, iLon, psi1onGrid(iLat,iLon)
		end do
	end do
	close(11)
	write(*,*) 'Save Psi'

	! ! !! Save AC 
	! open(unit=12, file='../Python/output/Fortran/EquationOne/Normalized/AC.csv', status='replace', action='write', form='formatted')
	! write(12,'(A)') 'row,col,value'  ! exact header
	! do row = 1, AC%m
	! 	if (row == 1) then
	! 		startIdx = 1
	! 	else
	! 		startIdx = AC%nentries(row-1) + 1
	! 	end if
	! 	endIdx = AC%nentries(row)
	! 	do k = startIdx, endIdx
	! 		! Write exactly with commas
	! 		write(12,'(I0,",",I0,",",ES24.16)') row, AC%j(k), AC%value(k)
	! 	end do
	! end do
	! close(12)

	! write(*,*) 'Saved AC.csv'

	deallocate(b)
	deallocate(sol)
	deallocate(res)

	return
end subroutine solvePoissonForPsiOne

! subroutine SolvingLibrary()
! 	use OceanGridDataModule
! 	use OceanAdaptiveGridDataModule
! 	use sparseMatricesDataModule
! 	use basicDataStructures
! 	use simulationParameters
! 	use lsqrModule
! 	use GMRESSettings
! 	use bathymetryDataModule
! 	use AnalysisModule !ADDED
! 	use petscksp
! 	!$ use omp_lib
! 	implicit none
! 	integer :: i,k,iLat,iLon,iM1,iP1,ki,m,nElements
! 	real (kind=QR_K) :: dLat, dLat2,lat,cosLat,sinLat,Lon, R2 , dphi ,dphi2
! 	real (kind=QR_K) :: cosLatMHalf,cosLatPHalf,cosLat2,sinLat2 , latPHalf, latMHalf
! 	double precision, allocatable :: b(:),res(:),sol(:)
! 	logical :: signConstrained
! 	integer :: i1,i2,CHUNK,nThreads,threadID , rasberries
! 	integer :: iLatP1,iLatM1,iLonP1,iLonM1,CorrespondingiLonM1,CorrespondingiLonP1, epsilon, index_equatorial , equatorial_particles, northern_hemisphere_particles, southern_hemisphere_particles
! 	double precision :: time1
! 	real (kind=QR_K) :: d2, coriolisTerm,dijOverd0
! 	integer :: successFlag !ADDED
! 	double precision :: normr !ADDED
! 	real (kind=QR_K) :: counter,alpha,beta, psiave , psimax
! 	! PETSc vars
! 	Mat A
! 	Vec pb, px
! 	KSP ksp
! 	PC pc
! 	PetscErrorCode ierr
! 	PetscInt N, rowstart, rowend, col, nnz_row_max
! 	PetscScalar val
! 	PetscReal rtol = 1.0d-6, abstol = 1.0d-8, norm_b
! 	integer :: row, startIdx, endIdx, jidx
! 	time1=omp_get_wtime()
! 	write(*,*)'Solving the Poisson equation for ψ₁: ∇ψ₁ - (1/Ld²) ψ₁ = - f(θ); θ₀=π/4 ; d₀= 1000 m',OceanGrid%nLat,OceanGrid%nPhi
! 	nElements=OceanGrid%nLat*OceanGrid%nPhi
! 	N = nElements  ! Size
! 	allocate(b(nElements))
! 	b = zero
! 	i=0
! 	R2=R*R
! 	dLat=OceanGrid%dLat
! 	dLat2=dLat*dLat
! 	dphi=OceanGrid%dphi
! 	dphi2=dphi*dphi
! 	do iLat=1,OceanGrid%nLat
! 	lat=-half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)
! 	cosLat=OceanGrid%cosLatC(iLat);
! 	cosLat2=cosLat*cosLat
! 	sinLat=OceanGrid%sinLatC(iLat)
! 	coriolisTerm = two*OMEGA*sinLat
! 	do iLon=1,OceanGrid%nPhi
! 	Lon=half*OceanGrid%dPhi+OceanGrid%dPhi*real(iLon-1,QR_K)
! 	i=i+1
! 	if (ScaleEquation) then
! 	!ROCKET: ∂λ²ψ + cosθ ∂θ·(cosθ ∂θ ψ) - (R²cos²θ Δθ²)·(f₀²/g*d₀) ψ = -f(θ)·(R²cos²θ Δθ²)
! 	b(i) = - coriolisTerm * ( dLat2 * R2 * cosLat2 )
! 	else
! 	!ROCKET: 1/(R²cos²θ)·∂²ψ/∂λ² + 1/(R²cosθ)·∂/∂θ·(cosθ ∂ψ/∂θ) - (f₀²/g*d₀) ψ = -f(θ)
! 	b(i) = - coriolisTerm
! 	end if
! 	!ADDED CONDTION: Apply Neumann boundary conditions on the boundary indices
! 	! if (iLat == 1 .or. iLat == OceanGrid%nLat) then
! 	! b(i) = zero
! 	! end if
! 	if(ISNAN(b(i)))stop "b nan"
! 	end do
! 	end do
	
! 	allocate(sol(nElements))
! 	allocate(res(nElements))
! 	sol = zero
! 	res = zero
! 	!TOFIX:
! 	signConstrained=.FALSE.
! 	! Initialize PETSc
! 	call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
! 	! Create matrix A
! 	call MatCreate(PETSC_COMM_WORLD, A, ierr)
! 	call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N, ierr)
! 	call MatSetFromOptions(A, ierr)
! 	nnz_row_max = 5  ! Max nonzeros per row (from stencil)
! 	call MatSeqAIJSetPreallocation(A, nnz_row_max, PETSC_NULL_INTEGER, ierr)
! 	! Fill A from AC (assuming 1-based indices; adjust to 0-based for PETSc)
! 	do row = 1, N
! 	if (row == 1) then
! 		startIdx = 1
! 	else
! 		startIdx = AC%nentries(row-1) + 1
! 	end if
! 	endIdx = AC%nentries(row)
! 	do jidx = startIdx, endIdx
! 		col = AC%j(jidx) - 1  ! 0-based
! 		val = AC%value(jidx)
! 		call MatSetValues(A, 1, row-1, 1, col, val, INSERT_VALUES, ierr)
! 	end do
! 	end do
! 	call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
! 	call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
! 	! Create vectors
! 	call VecCreate(PETSC_COMM_WORLD, pb, ierr)
! 	call VecSetSizes(pb, PETSC_DECIDE, N, ierr)
! 	call VecSetFromOptions(pb, ierr)
! 	call VecDuplicate(pb, px, ierr)
! 	! Set b
! 	do i = 0, N-1
! 	call VecSetValues(pb, 1, i, b(i+1), INSERT_VALUES, ierr)  ! b is 1-based
! 	end do
! 	call VecAssemblyBegin(pb, ierr)
! 	call VecAssemblyEnd(pb, ierr)
! 	! Adaptive tol: abstol based on ||b||
! 	call VecNorm(pb, NORM_2, norm_b, ierr)
! 	abstol = rtol * norm_b * 1.0d-2
! 	! Create KSP
! 	call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
! 	call KSPSetOperators(ksp, A, A, ierr)
! 	call KSPSetType(ksp, KSPBCGS, ierr)  ! BiCGSTAB
! 	call KSPSetTolerances(ksp, rtol, abstol, PETSC_DEFAULT_REAL, 10000, ierr)
! 	call KSPGetPC(ksp, pc, ierr)
! 	call PCSetType(pc, PCJACOBI, ierr)  ! Simple preconditioner; try ILU if needed
! 	! Solve
! 	call KSPSolve(ksp, pb, px, ierr)
! 	! Get residual norm
! 	call KSPGetResidualNorm(ksp, normr, ierr)
! 	write(*,*) 'Residual norm:', normr
! 	! Extract solution to sol
! 	call VecGetOwnershipRange(px, rowstart, rowend, ierr)
! 	do i = rowstart, rowend-1
! 	call VecGetValues(px, 1, i, val, ierr)
! 	sol(i+1) = val
! 	end do
! 	! Cleanup PETSc
! 	call MatDestroy(A, ierr)
! 	call VecDestroy(pb, ierr)
! 	call VecDestroy(px, ierr)
! 	call KSPDestroy(ksp, ierr)
! 	call PetscFinalize(ierr)
! 	!COMMENT: Solution of the streamfunction ψ₁
! 	if(.NOT.allocated(psi1onGrid))allocate(psi1onGrid(OceanGrid%nLat,OceanGrid%nPhi))
! 	psi1onGrid=zero
! 	i=0
! 	do iLat=1,OceanGrid%nLat
! 		do iLon=1,OceanGrid%nPhi
! 		i=i+1
! 			psi1onGrid(iLat,iLon)= sol(i)
! 		end do
! 	end do
! 	write(*,*)'POISSON SOLVER for Psi_1 sim time =', omp_get_wtime()-time1
! 	write(*,*)' '
! 	write(*,*) 'ψ: Max , Min = ', maxval(psi1onGrid), minval(psi1onGrid)
! 	end subroutine SolvingLibrary

!ADDED SUBROUTINE For Calculating the velocity u₁ = (∇ ⨉ ψ₁) r  = ∇ψ₁ ⨉ r ; where r is r hat.
subroutine VelocityFieldsOneOnGrid()
	use OceanGridDataModule
	use sparseMatricesDataModule
	use basicDataStructures
	use simulationParameters
	use lsqrModule
	use GMRESSettings
	use bathymetryDataModule
	!$ use omp_lib
	implicit none
	
	integer :: k, iLat, iLon
	real (kind=QR_K) :: dlat,dphi, dLat2, cosLat, sinLat, sinLatP1, sinLatM1,latPHalf,latMHalf
	real (kind=QR_K) :: cosLatMHalf, cosLatPHalf, cosLat2, sinLat2 ,VelocityMagnitude_PhiOne_average , velocityPsimax, lat,lon ,VelocityMagnitude_Psi_average
	integer :: i1, i2, CHUNK, nThreads, threadID
	integer :: iLatP1, iLatM1, iLonP1, iLonM1 , coconut, cranberries , blueberries, epsilon, index_equatorial,equatorial_particles
 
	write(*,*) ' Computing Velocity Fields u₁(θ,ƛ): u₁ = (∇ ⨉ ψ₁) r  = ∇ψ₁ ⨉ r '

	!ROCKET: uψ₁ (ƛ̂) = ( −(1/R) · (ψ₁(i,j+1) − ψ₁(i,j−1)) / (2Δθ) ƛ̂
	!ROCKET: uψ₁ (θ̂)  = ((1 / (R·cosθ(i,j))) · (ψ₁(i+1,j) − ψ₁(i−1,j)) / (2Δλ) θ̂ )

	if(.NOT.allocated(velocityPsiOneOnGrid))allocate(velocityPsiOneOnGrid(2,OceanGrid%nLat,OceanGrid%nPhi))
	velocityPsiOneOnGrid=zero
	if (.not.allocated(VelocityMagnitude_PsiOne)) allocate(VelocityMagnitude_PsiOne(OceanGrid%nLat, OceanGrid%nPhi))

	VelocityMagnitude_PsiOne = zero

	nThreads = OMP_get_max_threads()
	call OMP_SET_NUM_THREADS(nThreads)
	CHUNK = OceanGrid%nLat/nThreads

	write(*,*)CHUNK

	!$OMP PARALLEL DEFAULT(PRIVATE) &
	!$OMP& SHARED(CHUNK,OceanGrid,&
	!$OMP& velocityPsiOneOnGrid,&
	!$OMP& psi1onGrid,VelocityMagnitude_PsiOne)

	threadID = omp_get_thread_num()
	i1 = threadID*CHUNK+1
	i2=i1-1+CHUNK
	if(threadID+1.eq.nThreads)i2=OceanGrid%nLat

	do iLat=i1,i2
		cosLat=OceanGrid%cosLatC(iLat);
		sinLat=OceanGrid%sinLatC(iLat);

		!COMMENT: Neumann Boundary Conditons
		iLatM1=iLat-1
		if(iLat.eq.1)iLatM1=1
		iLatP1=iLat+1; 
		if(iLat.eq.OceanGrid%nLat)iLatP1=OceanGrid%nLat

		do iLon=1,OceanGrid%nPhi
			iLonM1=iLon-1
			if(iLon.eq.1)iLonM1=OceanGrid%nPhi
			iLonP1=iLon+1;
			if(iLon.eq.OceanGrid%nPhi)iLonP1=1

			!COMMENT: u₁(θ,ƛ) : u₁(1,:,:) ⇢ u₁(θ) , u₁(2,:,:) ⇢ u₁(ƛ)
			velocityPsiOneOnGrid(1,iLat,iLon)=(psi1onGrid(iLat,iLonP1)-psi1onGrid(iLat,iLonM1))/(two* OceanGrid%dphi*R*cosLat) 

			velocityPsiOneOnGrid(2,iLat,iLon)= -(psi1onGrid(iLatP1,iLon)-psi1onGrid(iLatM1,iLon))/(two* OceanGrid%dLat*R)

			!COMMENT: magnitude u₁(θ,ƛ) 
			VelocityMagnitude_PsiOne(iLat,iLon) = sqrt( velocityPsiOneOnGrid(1,iLat,iLon)**2 + velocityPsiOneOnGrid(2,iLat,iLon)**2 )
		end do
	end do
	!$OMP END PARALLEL

	write(*,*) 'Velocity Fields computations are done'
	write(*,*) 

	!ADDED STATEMENT: printing statements for the results min,max,magnitude
	write(*,*) ' '
	write(*,*)' ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯ '
	write(*,*) 'u₁(θ,ƛ): Max , Min = ', maxval(VelocityMagnitude_PsiOne), minval(VelocityMagnitude_PsiOne)

	! velocityPsimax= maxval(VelocityMagnitude_PsiOne)
	write(*,*) ' '
	write(*,*)' ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯ '
	write(*,*) ' '
	write(*,*) 'u₁(θ): Max, Min = ', maxval(velocityPsiOneOnGrid(1,:,:)), ',',minval(velocityPsiOneOnGrid(1,:,:))
	write(*,*) 'u₁(ƛ): Max, Min = ', maxval(velocityPsiOneOnGrid(2,:,:)), ',',minval(velocityPsiOneOnGrid(2,:,:))
	write(*,*)' ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯ '
	write(*,*) ' '

	!ADDED SECTION: The section below was used as bug tracking mechanism. IGNORE or INVESTIGATE.
	!ADDED: ______________________________________________________________________________________________________________

	!BUGTRACK: Allocating three arrays to save and view later 3 regions with degree epsilon: 
	!BUGTRACK: 1. Equatorial band with epsilon range 

	!BUGTRACK: 1. Equatorial band with epsilon range 
	! if(.NOT.allocated(velocitypsi1Band))allocate(velocitypsi1Band(OceanGrid%nLat,OceanGrid%nPhi))
	! velocitypsi1Band=zero

	! epsilon = bandie!twenty !zero!twelve!five
	! index_equatorial= OceanGrid%nLat/two
	! equatorial_particles=0
	! do iLat = index_equatorial - epsilon, index_equatorial + epsilon
	! 	do iLon = 1, OceanGrid%nPhi
	! 		velocitypsi1Band(iLat,iLon) = VelocityMagnitude_PsiOne(iLat, iLon)
	! 		equatorial_particles=equatorial_particles+1 
	! 	end do
	! end do
	! write(*,*) ' '
	! write(*,*)' --------------- '

	! write(*,*) 'Number of particles on the equatorial band of epsilon ', epsilon, 'is ', equatorial_particles, 'with Velocity psi average =', sum(velocitypsi1Band)/size(velocitypsi1Band), 'max Velocity psi = ', maxval(velocitypsi1Band), 'and min value = ', minval(velocitypsi1Band) , 'absolute value Velocity min = ', minval(abs(velocitypsi1Band))

	!BUGTRACK: Any zero values on the band.
	! write(*,*)'any Velocity psi zero values on band ?? -> ', any(velocitypsi1Band==zero)
	! write(*,*)' --------------- '
	! write(*,*) ' '

	!BUGTRACK: Track and count values equal to a certain value ⇢ here velocityPsimax
	!Blueberries coconuts
	! blueberries =0 
	! do iLat = 1,OceanGrid%nLat
	! 	do iLon = 1, OceanGrid%nPhi
	! 		if (VelocityMagnitude_PsiOne(iLat,iLon).eq. velocityPsimax) then 
	! 		lat=-half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)
	! 		lon=half*OceanGrid%dPhi+OceanGrid%dPhi*real(iLon-1,QR_K)
	! 		blueberries = blueberries+1
	! 		write(*,*) 'Number of particles with highest velocity is ', blueberries, 'of velocity =  ', velocityPsimax, ' iLat =',iLat, 'iLon = ', iLon , 'lat val in deg ',lat*180,'long val in deg ',lon*180, 'depth = ', OceanGrid%d(iLat,iLon), 'psi = ', psi1onGrid(iLat,iLon)
	! 		end if
	! 	end do
	! end do

	!BUGTRACK: Track and count values greater and less than a certain value ⇢ here taverage velocity
	! coconut=0 
	! blueberries =0
	! VelocityMagnitude_Psi_average= sum(VelocityMagnitude_PsiOne)/size(VelocityMagnitude_PsiOne) 
	! write(*,*) 'Average Velocity Psi = ', VelocityMagnitude_Psi_average
	! do iLat=1,OceanGrid%nLat 
	! 	do iLon=1,OceanGrid%nPhi 
	! 		if (VelocityMagnitude_PsiOne(iLat,iLon).gt.VelocityMagnitude_Psi_average) then 
	! 			coconut=coconut+1 
	! 		end if
	! 		if (VelocityMagnitude_PsiOne(iLat,iLon).lt.VelocityMagnitude_Psi_average) then 
	! 			blueberries=blueberries+1 
	! 		end if
	! 	end do
	! end do

	! write(*,*) 'Number of particles Greater than the average of Velocity psi = ', coconut
	! write(*,*) 'Number of particles Less than the average of Velocity psi = ', blueberries

	! write(*,*) 'at highest psi:      U_psi(383,718) = ', VelocityMagnitude_PsiOne(383,718) 

end subroutine VelocityFieldsOneOnGrid

!TOFIX: Subroutine (in progress — not yet ready for use or aligned with current workflow). Build the coefficient matrix for ψ₂ and append values to AdaptiveOceanGrid.The discretization scheme is shown below.   
subroutine OceanbuildCoefficientsMatrixPsiTwo()
	use OceanAdaptiveGridDataModule
	use sparseMatricesDataModule
	use basicDataStructures
	use simulationParameters
	use lsqrModule
	use GMRESSettings
	use bathymetryDataModule
	!$ use omp_lib
	implicit none
	
	integer :: i,k,iLat,iLon,iM1,iP1,nElements 
	real (kind=QR_K) :: dLat, dLat2,lat,cosLat,sinLat,Lon, R2, ScreeningTerm,Ld, Ld2 ,beta,gamma,alpha
	real (kind=QR_K) :: latMHalf,latPHalf ,cosLatMHalf,cosLatPHalf,cosLat2,sinLat2, sinLatM1, sinLatP1, R2OverLd02
	double precision, allocatable :: b(:),res(:),sol(:)
	real (kind=QR_K) :: coriolisTerm,coriolisTerm2, Ld02
	logical :: signConstrained
	integer :: i1,i2,CHUNK,nThreads,threadID
	integer :: iLatP1,iLatM1,iLonP1,iLonM1 ,startIdx, endIdx
	double precision :: time1
	real (kind=QR_K) :: dijOverd0
	
	time1=omp_get_wtime()
	
	write(*,*)' Dimensions of the Coefficient Matrix = ',OceanAdaptiveGrid%nLat,OceanAdaptiveGrid%nPhi

	AC%length = (OceanAdaptiveGrid%nLat-2) * OceanAdaptiveGrid%nPhi*5 + 2*OceanAdaptiveGrid%nPhi*4

	write(*,*) ' Length Assigned for the Coefficient Matrix = ', (OceanAdaptiveGrid%nLat-2) * OceanAdaptiveGrid%nPhi*5 + 2*OceanAdaptiveGrid%nPhi*4

	if(allocated(AC%j))deallocate(AC%j)
	if(allocated(AC%value))deallocate(AC%value)
	if(allocated(AC%nentries))deallocate(AC%nentries)
	allocate(AC%j(AC%length))
	allocate(AC%value(AC%length))
	nElements=OceanAdaptiveGrid%nLat*OceanAdaptiveGrid%nPhi
	allocate(AC%nentries(nElements))
	AC%value=zero

	i=0 ; k=0

	dLat=OceanAdaptiveGrid%dLat
	dLat2=dLat*dLat
	R2= R*R

	!COMMENT: This subroutine solves the following equation: 
	!COMMENT:               ∇ψ₂ − (1/Ld²) ψ₂ = − Q(θ,ƛ)*d
	!COMMENT: In Progress this subroutine...
	!ROCKET: The Stencil(ƛ,θ):                                              
	!ROCKET:                     ψ₂(ƛᵢ,θᵢ₊₁)                              
	!ROCKET:                          |                                   
	!ROCKET:    ψ₂(ƛᵢ₋₁, θᵢ) ___ ψ₂(ƛᵢ,θᵢ) ___ ψ₂(ƛᵢ₊₁, θᵢ)              
	!ROCKET:                          |                                   
	!ROCKET:                     ψ₂(ƛᵢ,θᵢ₋₁)                               
	!COMMENT:  * bottom = 
	!COMMENT:  * top = 
	!COMMENT:  * left = right = 

	! Ld02 = Ld0*Ld0
	! R2OverLd02 = (R/Ld0)**2 !R2/Ld0
	! write(*,*) 'Ld0 , Ld02 ,R2 , alpha2 = ', Ld0 ,Ld0*Ld0, R*R , alpha

	write(*,*)' Building the Coefficients Matrix for ψ₂: ∇ψ₂ - (1/Ld²) ψ₂ =  Q(θ,ƛ)*d '
	
	do iLat=1,OceanAdaptiveGrid%nLat
		lat= -half*pi+half*OceanAdaptiveGrid%dLat+OceanAdaptiveGrid%dLat*real(iLat-1,QR_K)
	
		cosLat=OceanAdaptiveGrid%cosLatC(iLat)
		cosLat2=cosLat*cosLat

		latPHalf = Lat+half*dLat
		latMHalf = Lat -half*dLat
		
		cosLatMHalf=cos(latMHalf) 
		cosLatPHalf=cos(latPHalf) 
		if(iLat.eq.1) cosLatMHalf=zero
		if(iLat.eq.OceanAdaptiveGrid%nLat) cosLatPHalf=zero

		sinLat=OceanAdaptiveGrid%sinLatC(iLat)
		coriolisTerm = two*OMEGA*sinLat		!two*OMEGA!*sinLat !f_equatorial! two*OMEGA*sinLat
		coriolisTerm2= coriolisTerm*coriolisTerm  

		do iLon=1,OceanAdaptiveGrid%nPhi
			i=i+1
			Lon=half*OceanAdaptiveGrid%dPhi+OceanAdaptiveGrid%dPhi*real(iLon-1,QR_K)

			if (OceanAdaptiveGrid%d(iLat,iLon).gt.zero) then 
				dijOverd0= OceanAdaptiveGrid%d(iLat,iLon) / OceanAdaptiveGrid%d0 
			else 
				dijOverd0= zero
			end if
			
				iM1=i-1; if(iLon.eq.1) iM1=i-1+OceanAdaptiveGrid%nPhi
				iP1=i+1; if(iLon.eq.OceanAdaptiveGrid%nPhi) iP1=i+1-OceanAdaptiveGrid%nPhi

				if(iLat.gt.1)then !!Latex: psi_{i,j-1} 
					k=k+1;
					AC%j(k)=i-OceanAdaptiveGrid%nPhi
					AC%value(k)=dijOverd0*cosLat*cosLatMHalf
				end if

				k=k+1; !! Latex: psi_{i-1,j}
				AC%j(k)=iM1
				AC%value(k)=dijOverd0

				k=k+1;  !!psi_{i,j} - diagonal term 
				AC%j(k)=i
				AC%value(k)=dijOverd0*(-two- cosLat*(cosLatMHalf+cosLatPHalf)) -((coriolisTerm2*dLat2*R2*cosLat2)/ (GRAVITY*OceanAdaptiveGrid%d0))
				
				!- (dLat2*cosLat2*sinLat2*R2OverLd02)
				
				!((coriolisTerm2*dLat2*R2*cosLat2)/ (GRAVITY*OceanAdaptiveGrid%d0))
				
				k=k+1; !!Latex: psi_{i+1,j}
				AC%j(k)=iP1
				AC%value(k)=dijOverd0
				
				if(iLat.lt.OceanAdaptiveGrid%nLat)then !!Latex: psi_{i,j+1}
					k=k+1;
					AC%j(k)=i+OceanAdaptiveGrid%nPhi
					AC%value(k)=dijOverd0*cosLat*cosLatPHalf
				end if 
				
				AC%nentries(i) = k
				
			! end if
			end do
		end do

		AC%m = nElements
		AC%n = nElements
		AC%rowIndexed = 1

	write(*,*) 'AC Matrix building for Psi_1 completed in ', omp_get_wtime()-time1, ' seconds'

end subroutine OceanbuildCoefficientsMatrixPsiTwo

!TOFIX: Subroutine (in progress — not yet ready for use or aligned with current workflow). Build the coefficient matrix for ϕ₂ and append values to AdaptiveOceanGrid.The discretization scheme is shown below.   
subroutine OceanbuildCoefficientsMatrixPhiOne()
	use OceanGridDataModule
	use sparseMatricesDataModule
	use basicDataStructures
	use simulationParameters
	use lsqrModule
	use GMRESSettings
	use bathymetryDataModule
	!$ use omp_lib
	implicit none
	
	integer :: i,k,iLat,iLon,iM1,iP1,nElements 
	real (kind=QR_K) :: dLat, dLat2,lat,cosLat,sinLat,Lon, R2
	real (kind=QR_K) :: latMHalf,latPHalf,cosLatMHalf,cosLatPHalf,cosLat2,sinLat2, sinLatM1, sinLatP1
	double precision, allocatable :: b(:),res(:),sol(:)
	real (kind=QR_K) :: coriolisTerm,coriolisTerm2, coriolisTermP1, coriolisTermM1
	logical :: signConstrained
	integer :: i1,i2,CHUNK,nThreads,threadID 
	integer :: iLatP1,iLatM1,iLonP1,iLonM1,CorrespondingiLonM1, CorrespondingiLonP1

	double precision :: time1
	!!parameters
	real (kind=QR_K) :: dij, diP1j, diM1j, dijP1, dijM1 
	real (kind=QR_K) :: hs1ij,hs1iM1j, hs1iP1j,hs1ijP1,hs1ijM1 , h1ij, h1iM1j, h1iP1j, h1ijP1, h1ijM1
	real (kind=QR_K) :: psi1onGridij, psi1onGridiP1j, psi1onGridiM1j, psi1onGridijP1,psi1onGridijM1
	
	time1=omp_get_wtime()

	!! the coefficient matrix of: Laplacian_s(Psi) - f^2(theta)/gd(theta,phi) Psi(theta,phi) = -f(theta)
	write(*,*)'Filling the Matrix for Phi1_1 ',OceanGrid%nLat,OceanGrid%nPhi

	BC%length = (OceanGrid%nLat-2) * OceanGrid%nPhi*5 + 2*OceanGrid%nPhi*4
	write(*,*) 'Length Assigned for the BC matrix is ', (OceanGrid%nLat-2) * OceanGrid%nPhi*5 + 2*OceanGrid%nPhi*4

	!! BC for psi1
	if(allocated(BC%j))deallocate(BC%j)
	if(allocated(BC%value))deallocate(BC%value)
	if(allocated(BC%nentries))deallocate(BC%nentries)
	allocate(BC%j(BC%length))
	allocate(BC%value(BC%length))
	nElements=OceanGrid%nLat*OceanGrid%nPhi
	allocate(BC%nentries(nElements))
	BC%value=zero 

	! !!Allocating h1 and hs1
	if(allocated(OceanGrid%hs1))deallocate(OceanGrid%hs1)
	if(.NOT.allocated(OceanGrid%hs1))allocate(OceanGrid%hs1(OceanGrid%nLat,OceanGrid%nPhi))
	OceanGrid%hs1=zero

	if(allocated(OceanGrid%h1))deallocate(OceanGrid%h1)
	if(.NOT.allocated(OceanGrid%h1))allocate(OceanGrid%h1(OceanGrid%nLat,OceanGrid%nPhi))

	OceanGrid%h1=zero
	
	i=0 ; k=0

	dLat=OceanGrid%dLat
	dLat2=dLat*dLat
	R2= R*R

	write(*,*)' Building the Coefficients Matrix for Phi_1 '
	
	do iLat=1,OceanGrid%nLat
		lat= -half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)

		cosLat=OceanGrid%cosLatC(iLat); 
		cosLat2=cosLat*cosLat

		latPHalf = Lat+half*dLat
		latMHalf = Lat -half*dLat

		cosLatMHalf=cos(latMHalf) 
		cosLatPHalf=cos(latPHalf)

		if(iLat.eq.1) cosLatMHalf=zero
		if(iLat.eq.OceanGrid%nLat) cosLatPHalf=zero

		! latP1 = lat + dLat 
		! latM1 = lat - dLat

		! if(latP1>pi/2)then 
		! 	latP1 = pi-latP1 
		! 	ReflonP1 = pi+lon
		! end if
		! if(latM1<-pi/2)then
		! 	latM1 = -pi-latM1 
		! 	ReflonM1 = pi+lon
		! end if

		! if(lonP1>2*pi) then
		! 	ReflonP1 = lonP1 - 2*pi
		! end if
		! if(lonP1<0) then
		! 	ReflonP1 = lonP1 + 2*pi
		! end if

		iLatP1 = iLat + 1
		if(iLatP1.gt.OceanGrid%nLat) iLatP1 = iLat 
		iLatM1 = iLat - 1
		if(iLatM1.lt.1) iLatM1 = iLat

		sinLatP1 = OceanGrid%sinLatC(iLatP1)
		sinLatM1 = OceanGrid%sinLatC(iLatM1)

		sinLat = OceanGrid%sinLatC(iLat)
		
		coriolisTerm = two * OMEGA * sinLat !f_equatorial !two * OMEGA * sinLat

		coriolisTermP1 = two * OMEGA * sinLatP1 !f_equatorial!two * OMEGA * sinLatP1
		coriolisTermM1 = two * OMEGA * sinLatM1 !f_equatorial! two * OMEGA * sinLatM1
		
		do iLon=1,OceanGrid%nPhi
			i=i+1
			Lon=half*OceanGrid%dPhi+OceanGrid%dPhi*real(iLon-1,QR_K)

		!Getting h1 and hs1 values
			if (OceanGrid%d(iLat,iLon).gt.zero) then !if ocean particle
				OceanGrid%hs1(iLat,iLon) = (coriolisTerm*psi1onGrid(iLat,iLon))/GRAVITY
				OceanGrid%h1(iLat,iLon) = OceanGrid%d(iLat,iLon)+ OceanGrid%hs1(iLat,iLon)
	
			end if
			
			iLonP1 = iLon + 1; if(iLon.eq.OceanGrid%nPhi) iLonP1 = 1
			iLonM1 = iLon - 1; if(iLon.eq.1) iLonM1 = OceanGrid%nPhi
			iLatP1 = iLat + 1
	        CorrespondingiLonP1 = iLon
			
			if(iLatP1.gt.OceanGrid%nLat) then
				iLatP1 = iLat 
				CorrespondingiLonP1 = OceanGrid%nPhi/2 +iLon !!meaning pi + lambda
			end if
			!* (theta-1, lambda) or (Reflected theta-1, reflected lambda)
			iLatM1 = iLat - 1
			CorrespondingiLonM1 = iLon
			if(iLatM1.lt.1)then 
				iLatM1 = iLat
				CorrespondingiLonM1 = OceanGrid%nPhi/2 +iLon
			end if

			if(CorrespondingiLonP1 > OceanGrid%nPhi) CorrespondingiLonP1 = 1
			if(CorrespondingiLonM1 > OceanGrid%nPhi) CorrespondingiLonM1 = 1
			if(CorrespondingiLonM1 < 1) CorrespondingiLonM1 = OceanGrid%nPhi
			if(CorrespondingiLonP1 < 1) CorrespondingiLonP1 = OceanGrid%nPhi

			if (OceanGrid%d(iLat,iLon).gt.zero) then
				dij= OceanGrid%d(iLat,iLon) !/ OceanGrid%d0
				psi1onGridij = psi1onGrid(iLat,iLon)
				hs1ij= (coriolisTerm*psi1onGridij)/GRAVITY
				h1ij = (dij + hs1ij) !/ OceanGrid%d0

				! write(*,*) 'Depth on the Grid at (',iLat,',',iLon,') = ', dij
				! write(*,*) 'Free surface on the Grid at (',iLat,',',iLon,') = ', hs1ij
				! write(*,*) 'Height on the Grid at (',iLat,',',iLon,') = ', h1ij

				if (OceanGrid%d(iLat,iLonP1).gt.zero) then 
					diP1j = OceanGrid%d(iLat,iLonP1) !/ OceanGrid%d0
					psi1onGridiP1j = psi1onGrid(iLat,iLonP1)
					hs1iP1j= (coriolisTerm*psi1onGridiP1j)/GRAVITY
					h1iP1j = (diP1j + hs1iP1j) !/ OceanGrid%d0
				else
					diP1j = zero
					psi1onGridiP1j = zero
					hs1iP1j=0
					h1iP1j=0
				end if

				if (OceanGrid%d(iLat,iLonM1).gt.zero) then
					!!Latex d(i-1,j) = d(theta,long-1)
					diM1j = OceanGrid%d(iLat,iLonM1) !/ OceanGrid%d0
					psi1onGridiM1j = psi1onGrid(iLat,iLonM1)
					hs1iM1j= (coriolisTerm*psi1onGridiM1j)/GRAVITY
					h1iM1j = (diM1j + hs1iM1j) !/ OceanGrid%d0
				else
					diM1j = zero
					psi1onGridiM1j = zero
					hs1iM1j=0
					h1iM1j=0
				end if

				if (OceanGrid%d(iLatP1,CorrespondingiLonP1).gt.zero) then
					!!Latex d(i,j+1) = d(theta+1,long)
					dijP1 = OceanGrid%d(iLatP1,CorrespondingiLonP1) !/ OceanGrid%d0
					psi1onGridijP1 = psi1onGrid(iLatP1,CorrespondingiLonP1)
					hs1ijP1= (coriolisTermP1*psi1onGridijP1)/GRAVITY
					h1ijP1 = (dijP1 + hs1ijP1) !/ OceanGrid%d0
				else
					dijP1 = zero
					psi1onGridijP1 = zero
					hs1ijP1=0
					h1ijP1=0
				end if

				if (OceanGrid%d(iLatM1,CorrespondingiLonM1).gt.zero) then
					!!Latex d(i,j-1) = d(theta-1,long)
					dijM1 = OceanGrid%d(iLatM1,CorrespondingiLonM1) !/ OceanGrid%d0
					psi1onGridijM1 = psi1onGrid(iLatM1,CorrespondingiLonM1)
					hs1ijM1= (coriolisTermM1*psi1onGridijM1)/GRAVITY
					h1ijM1 = (dijM1 + hs1ijM1) !/ OceanGrid%d0
				else
					dijM1 = zero
					psi1onGridijM1 = zero
					hs1ijM1=0
					h1ijM1=0
				end if	
			else
				dij      = zero 
				diP1j    = zero
				diM1j    = zero
				dijP1    = zero
				dijM1    = zero
				psi1onGridij   = zero
				psi1onGridiP1j = zero
				psi1onGridiM1j = zero
				psi1onGridijP1 = zero
				psi1onGridijM1 = zero
				hs1ij   = zero
				hs1iP1j = zero
				hs1iM1j = zero
				hs1ijP1 = zero
				hs1ijM1 = zero
				h1ij    = zero
				h1iP1j  = zero
				h1iM1j  = zero
				h1ijP1  = zero
				h1ijM1  = zero
			end if

			iM1=i-1; if(iLon.eq.1) iM1=i-1+OceanGrid%nPhi
			iP1=i+1; if(iLon.eq.OceanGrid%nPhi) iP1=i+1-OceanGrid%nPhi

			if(iLat.gt.1)then !!Latex: phi_{i,j+1} 
				k=k+1;
				BC%j(k)=i-OceanGrid%nPhi
				BC%value(k)=  h1ij*cosLat*cosLatPHalf + (cosLat2*quarter)*(h1ijP1-h1ijM1)
			end if

			k=k+1; !! Latex: phi_{i-1,j}
			BC%j(k)=iM1
			BC%value(k)= h1ij -quarter*(h1iP1j - h1iM1j)

			k=k+1;  !!phi_{i,j} - diagonal term
			BC%j(k)=i
			BC%value(k)= h1ij*(-two - cosLat*cosLatPHalf - cosLat*cosLatMHalf)
			
			k=k+1; !!Latex: phi_{i+1,j}
			BC%j(k)=iP1
			BC%value(k)= h1ij +quarter*(h1iP1j - h1iM1j)
			
			if(iLat.lt.OceanGrid%nLat)then !!Latex: phi_{i,j-1}
				k=k+1;
				BC%j(k)=i+OceanGrid%nPhi
				BC%value(k)= h1ij*coslat*cosLatMHalf -(cosLat2*quarter )* (h1ijP1 - h1ijM1)
			end if 
			
			BC%nentries(i) = k
		end do
	end do

	BC%m = nElements
	BC%n = nElements
	BC%rowIndexed = 1

	write(*,*) 'BC Matrix building for Phi_1 completed in ', omp_get_wtime()-time1, ' seconds'

	write(*,*) 'hmax , hmin = ',  &
	maxval(OceanGrid%h1,  mask = OceanGrid%d > zero), &
	minval(OceanGrid%h1,  mask = OceanGrid%d > zero)

	write(*,*) 'hsmax , hsmin = ', &
	maxval(OceanGrid%hs1, mask = OceanGrid%d > zero), &
	minval(OceanGrid%hs1, mask = OceanGrid%d > zero)

	write(*,*) 'dmax , dmin = ', maxval(OceanGrid%d, mask = OceanGrid%d > zero), minval(OceanGrid%d, mask = OceanGrid%d > zero)

end subroutine OceanbuildCoefficientsMatrixPhiOne

!TOFIX: Subroutine (in progress — not yet ready for use or aligned with current workflow). 
!TOFIX: Purpose: Solving for ψ₂.
subroutine solvePoissonForPsiTwo()
	use OceanAdaptiveGridDataModule 
	use sparseMatricesDataModule
	use basicDataStructures
	use simulationParameters
	use lsqrModule
	use GMRESSettings
	use bathymetryDataModule
	!$ use omp_lib
	implicit none

	integer :: i,k,iLat,iLon,iM1,iP1,ki,m,nElements
	real (kind=QR_K) :: dLat, dLat2,lat,cosLat,sinLat,Lon, R2 , dphi ,dphi2
	real (kind=QR_K) :: cosLatMHalf,cosLatPHalf,cosLat2,sinLat2 , latPHalf, latMHalf
	double precision, allocatable :: b(:),res(:),sol(:)
	logical :: signConstrained
	integer :: i1,i2,CHUNK,nThreads,threadID , rasberries
	integer :: iLatP1,iLatM1,iLonP1,iLonM1,CorrespondingiLonM1,CorrespondingiLonP1, epsilon, index_equatorial , equatorial_particles, northern_hemisphere_particles, southern_hemisphere_particles
	double precision :: time1
	real (kind=QR_K) :: d2, coriolisTerm,dijOverd0, term1,term2

	real (kind=QR_K) :: counter,alpha,beta, psiave , psimax

	time1=omp_get_wtime()

	write(*,*)'Solving the Poisson equation for ψ₁: ∇ψ₁ - (1/Ld²) ψ₁ = - Q*d  ',OceanAdaptiveGrid%nLat,OceanAdaptiveGrid%nPhi
	nElements=OceanAdaptiveGrid%nLat*OceanAdaptiveGrid%nPhi
	allocate(b(nElements))
	b = zero

	i=0
	R2=R*R
	dLat=OceanAdaptiveGrid%dLat
	dLat2=dLat*dLat

	dphi=OceanAdaptiveGrid%dphi
	dphi2=dphi*dphi

	!! For the RHS
	do iLat=1,OceanAdaptiveGrid%nLat
		lat=-half*pi+half*OceanAdaptiveGrid%dLat+OceanAdaptiveGrid%dLat*real(iLat-1,QR_K)

		cosLat=OceanAdaptiveGrid%cosLatC(iLat);
		cosLat2=cosLat*cosLat
		sinLat=OceanAdaptiveGrid%sinLatC(iLat)

		coriolisTerm = two*OMEGA*sinLat !two*OMEGA !*sinLat !f_equatorial!two*OMEGA*sinLat

		do iLon=1,OceanAdaptiveGrid%nPhi
			Lon=half*OceanAdaptiveGrid%dPhi+OceanAdaptiveGrid%dPhi*real(iLon-1,QR_K)
			i=i+1
			if (OceanAdaptiveGrid%d(iLat,iLon).gt.zero)then
				!!Boundary Conditions
				dijOverd0 = OceanAdaptiveGrid%d(iLat,iLon)/OceanAdaptiveGrid%d0 

				!!psi_2: = Qd 
				b(i) = dijOverd0*R2*cosLat2*dLat2*coriolisTerm
				
				!OceanAdaptiveGrid%Q(iLat,iLon) * OceanAdaptiveGrid%d(iLat,iLon)

				! if(b(i).gt.0)write(*,*)'b(i) is not 0 = ', b(i)
			else
				b(i) = zero
			end if
			if(ISNAN(b(i)))stop "b nan"
		end do
	end do

	allocate(sol(nElements))
	allocate(res(nElements))
	sol = zero
	res = zero

	signConstrained=.FALSE.
	!		call GMRES_RESTARTED_OMP_RETURNR(AC,b(1:nElements),sol(1:nElements),signConstrained,&
	!		&res(1:nElements),nElements,GMRES_ninner,GMRES_nouter,GMRES_tol,successFlag,normr)

	!BCG(A,b,x,m,n,nIter,res)
	!write(*,*)'calling BICGSTAB_OMP'

	call BICGSTAB_OMP(AC,b(1:nElements),sol(1:nElements),nElements,GMRES_nouter,GMRES_tol)

	if(.NOT.allocated(psi2onGrid))allocate(psi2onGrid(OceanAdaptiveGrid%nLat,OceanAdaptiveGrid%nPhi))
	psi2onGrid=zero

	write(*,*)'BICGSTAB_OMP returned successfully'

	i=0
	do iLat=1,OceanAdaptiveGrid%nLat
		do iLon=1,OceanAdaptiveGrid%nPhi
			i=i+1
			psi2onGrid(iLat,iLon)= sol(i)
			! write(*,*)'Psi_1 on the Grid at (',iLat,',',iLon,') = ', psi2onGrid(iLat,iLon)
		end do
	end do

	psimax=maxval(psi2onGrid) 

	write(*,*)'POISSON SOLVER for Psi 2 sim time =', omp_get_wtime()-time1

	write(*,*) 'Psi max, Psi min = ', maxval(psi2onGrid), minval(psi2onGrid)

	psiave = sum(abs(psi2onGrid))/size(psi2onGrid)
	write(*,*) 'Psi average = ', psiave
	
	deallocate(b)
	deallocate(sol)
	deallocate(res)

	return

end subroutine solvePoissonForPsiTwo

subroutine solvePoissonForPhiOne()
	use OceanGridDataModule
	use sparseMatricesDataModule
	use basicDataStructures
	use simulationParameters
	use lsqrModule
	use GMRESSettings
	use bathymetryDataModule
	!$ use omp_lib
	implicit none
	
	integer :: i, k, iLat, iLon, iM1, iP1, ki, m, nElements
	real (kind=QR_K) :: dLat, dLat2, lat, cosLat, sinLat, Lon, R2, sinLatP1, sinLatM1, dphi
	real (kind=QR_K) :: cosLatMHalf, cosLatPHalf, cosLat2, sinLat2
	double precision, allocatable :: b(:), res(:), sol(:)
	logical :: signConstrained
	integer :: i1, i2, CHUNK, nThreads, threadID
	integer :: iLatP1, iLatM1, iLonP1, iLonM1, CorrespondingiLonP1, CorrespondingiLonM1
	double precision :: time1
	real (kind=QR_K) :: d2, coriolisTerm,coriolisTermP1,coriolisTermM1
	!! new
	real (kind=QR_K) :: diP1j ,dijP1, diM1j,dijM1
	real (kind=QR_K) ::  psi1onGridiP1j, psi1onGridijP1,psi1onGridijM1,psi1onGridiM1j
	real(kind = QR_K) :: h1ijM1, h1ijP1 , h1iP1j, h1iM1j ,hs1ijM1, hs1ijP1 , hs1iP1j, hs1iM1j
	
	real(kind = QR_K) :: bmax, bave , dh1_dlambda, dh1_dtheta
	integer :: nActiveCells, counterIssam 

	time1=omp_get_wtime()

	write(*,*)'solving the Poisson equation for Phi_1 ',OceanGrid%nLat,OceanGrid%nPhi
	nElements=OceanGrid%nLat*OceanGrid%nPhi

	allocate(b(nElements))
	b = zero

	bmax = zero
	bave = zero 
	nActiveCells = 0
	counterIssam = 0
	
	if(.NOT.allocated(GradH1))allocate(GradH1(2,OceanGrid%nLat,OceanGrid%nPhi))
	GradH1=zero
	if(.NOT.allocated(uDotGradH1))allocate(uDotGradH1(OceanGrid%nLat,OceanGrid%nPhi))
	uDotGradH1=zero

	i=0
	R2=R*R
	dLat=OceanGrid%dLat
	dLat2=dLat*dLat
	dphi = OceanGrid%dphi
	do iLat=1,OceanGrid%nLat
		lat=-half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)

		cosLat=OceanGrid%cosLatC(iLat);
		sinLat=OceanGrid%sinLatC(iLat)

		coriolisTerm = two*OMEGA*sinLat !f_equatorial !two*OMEGA*sinLat
		
		!!index lat and lon: 
		iLatP1 = iLat + 1
		if(iLatP1.gt.OceanGrid%nLat) iLatP1 = iLat 
		iLatM1 = iLat - 1
		if(iLatM1.lt.1) iLatM1 = iLat

		sinLatP1 = OceanGrid%sinLatC(iLatP1)
		sinLatM1 = OceanGrid%sinLatC(iLatM1)

		coriolisTermP1 = two * OMEGA * sinLatP1 !f_equatorial!two * OMEGA * sinLatP1
		coriolisTermM1 = two * OMEGA * sinLatM1 !f_equatorial!two * OMEGA * sinLatM1

		do iLon=1,OceanGrid%nPhi
			Lon=half*OceanGrid%dPhi+OceanGrid%dPhi*real(iLon-1,QR_K)

			!!BC on iLon pm 1
			iLonP1=iLon+1; iLonM1=iLon-1
			if(iLon.eq.1)iLonM1=OceanGrid%nPhi
			if(iLon.eq.OceanGrid%nPhi)iLonP1=1

			!!Reflection:
			!! (theta+1, lambda) or (Reflected theta+1, reflected lambda)
			iLatP1 = iLat + 1 ; CorrespondingiLonP1 = iLon
			if(iLatP1.gt.OceanGrid%nLat) then
				iLatP1 = iLat 
				CorrespondingiLonP1 = OceanGrid%nPhi/2 +iLon
			end if
			!! (theta-1, lambda) or (Reflected theta-1, reflected lambda)
			iLatM1 = iLat - 1; CorrespondingiLonM1 = iLon
			if(iLatM1.lt.1)then 
				iLatM1 = iLat
				CorrespondingiLonM1 = OceanGrid%nPhi/2 +iLon
			end if
			if(CorrespondingiLonP1 > OceanGrid%nPhi)CorrespondingiLonP1 = 1
			if(CorrespondingiLonM1 < 1) CorrespondingiLonM1 = OceanGrid%nPhi

			i=i+1
			if (OceanGrid%d(iLat,iLon).gt.zero)then

				if (OceanGrid%d(iLat,iLonP1).gt.zero) then 
					!!Latex (i+1,j) = (theta,long+1)
					diP1j = OceanGrid%d(iLat,iLonP1) !/ OceanGrid%d0
					psi1onGridiP1j = psi1onGrid(iLat,iLonP1)
					hs1iP1j= coriolisTerm*psi1onGridiP1j/GRAVITY
					h1iP1j = (diP1j+ hs1iP1j) !/ OceanGrid%d0
				else
					diP1j = zero
					psi1onGridiP1j = zero
					hs1iP1j=0
					h1iP1j=0
				end if
				if (OceanGrid%d(iLat,iLonM1).gt.zero) then
					!!Latex (i-1,j) = (theta,long-1)
					diM1j = OceanGrid%d(iLat,iLonM1) !/ OceanGrid%d0
					psi1onGridiM1j = psi1onGrid(iLat,iLonM1)
					hs1iM1j= coriolisTerm*psi1onGridiM1j/GRAVITY
					h1iM1j = (diM1j + hs1iM1j) !/ OceanGrid%d0
				else
					diM1j = zero
					psi1onGridiM1j = zero
					hs1iM1j=0
					h1iM1j=0
				end if

				if (OceanGrid%d(iLatP1,CorrespondingiLonP1).gt.zero) then
					!!Latex (i,j+1) = (theta+1,long)
					dijP1 = OceanGrid%d(iLatP1,CorrespondingiLonP1) !/ OceanGrid%d0
					psi1onGridijP1 = psi1onGrid(iLatP1,CorrespondingiLonP1)
					hs1ijP1= coriolisTermP1*psi1onGridijP1/GRAVITY
					h1ijP1 = (dijP1 + hs1ijP1) !/ OceanGrid%d0
				else
					dijP1 = zero
					psi1onGridijP1 = zero
					hs1ijP1=0
					h1ijP1=0
				end if

				if (OceanGrid%d(iLatM1,CorrespondingiLonM1).gt.zero) then
					!!Latex (i,j-1) = (theta-1,long)
					dijM1 = OceanGrid%d(iLatM1,CorrespondingiLonM1) !/ OceanGrid%d0
					psi1onGridijM1 = psi1onGrid(iLatM1,CorrespondingiLonM1)
					hs1ijM1= coriolisTermM1*psi1onGridijM1/GRAVITY
					h1ijM1 = (dijM1 + hs1ijM1) !/ OceanGrid%d0
				else
					dijM1 = zero
					psi1onGridijM1 = zero
					hs1ijM1=0
					h1ijM1=0
				end if
				!!Solve	the Poisson equation for Phi_1
				b(i) =  -cosLat *(psi1OnGridiP1j - psi1OnGridiM1j) *(h1ijP1 - h1ijM1)+cosLat * (psi1OnGridijP1 - psi1OnGridijM1)*(h1iP1j - h1iM1j)

				GradH1(1,iLat,iLon)= ((h1ijP1 - h1ijM1) / (two* dLat))!*(R2*cosLat2*dphi**2)
				GradH1(2,iLat,iLon) = ((h1iP1j - h1iM1j) / (two* dphi))!*(R2*cosLat2*dphi**2)

				uDotGradH1(iLat,iLon) = b(i)
				
				!!b01.txt: saving all the values and iterations
				! write(*,*) 'ilat = ',iLat,'ilon = ',iLon, 'b(i) = ', b(i), 'h1ijP1 = ', h1ijP1, 'h1ijM1 = ', h1ijM1, 'h1iP1j = ', h1iP1j, 'h1iM1j = ', h1iM1j, 'psi1onGridiP1j = ', psi1onGridiP1j, 'psi1onGridiM1j = ', psi1onGridiM1j, 'psi1onGridijP1 = ', psi1onGridijP1, 'psi1onGridijM1 = ', psi1onGridijM1, 'cosLat = ', cosLat

				! DEB UG GING
				! if(abs(b(i)).gt.bmax) bmax= abs(b(i))
				! if(abs(b(i)).gt.1.0*2826324.634321)then 
				! 	! write(*,*) 'ilat = ',iLat,'ilon = ',iLon
				! 	counterIssam = counterIssam + 1
				! end if 	
				! bave = bave + abs(b(i))
				! nActiveCells = nActiveCells + 1
			
			else
				b(i) = zero
			end if
			if(ISNAN(b(i)))stop "b nan"
		end do
	end do

	! bave = bave/nActiveCells
	! write(*,*)'-----------------------------------'
	! write(*,*) 'bmax , bave = ',bmax, bave
	! write(*,*) 'counterIssam , nActiveCells= ',counterIssam, nActiveCells
	! stop	 

	allocate(sol(nElements))
	allocate(res(nElements))
	sol = zero
	res = zero

	signConstrained=.FALSE.

	call BICGSTAB_OMP(BC,b(1:nElements),sol(1:nElements),nElements,GMRES_nouter,GMRES_tol)

	if(.NOT.allocated(phi1onGrid))allocate(phi1onGrid(OceanGrid%nLat,OceanGrid%nPhi))
	phi1onGrid=zero

	i=0
	do iLat=1,OceanGrid%nLat
		do iLon=1,OceanGrid%nPhi
			i=i+1
			phi1onGrid(iLat,iLon)=sol(i)
			! write(*,*) 'Phi_1 on the Grid at (',iLat,',',iLon,') = ', phi1onGrid(iLat,iLon)
		end do
	end do
	write(*,*)'POISSON SOLVER for phi1 sim time=',omp_get_wtime()-time1

	write(*,*) 'Phi max, Phi min = ', maxval(phi1onGrid) , minval(phi1onGrid)

	deallocate(b)
	deallocate(sol)
	deallocate(res)
	
	return

end subroutine solvePoissonForPhiOne

subroutine VelocityFieldsTwoOnGrid()
	use OceanAdaptiveGridDataModule
	use sparseMatricesDataModule
	use basicDataStructures
	use simulationParameters
	use lsqrModule
	use GMRESSettings
	use bathymetryDataModule
	!$ use omp_lib
	implicit none
	
	integer :: k, iLat, iLon
	real (kind=QR_K) :: dlat,dphi, dLat2, cosLat, sinLat, sinLatP1, sinLatM1,latPHalf,latMHalf
	real (kind=QR_K) :: cosLatMHalf, cosLatPHalf, cosLat2, sinLat2 ,VelocityMagnitude_PhiOne_average , velocityPsimax, lat,lon ,VelocityMagnitude_Psi_average
	integer :: i1, i2, CHUNK, nThreads, threadID
	integer :: iLatP1, iLatM1, iLonP1, iLonM1 , coconut, cranberries , blueberries
 
	!!Computing velocity using finite difference of phi
	write(*,*) ' Computing Velocity Fields '

	!!Deltas

	if(.NOT.allocated(velocityPsiTwoOnGrid))allocate(velocityPsiTwoOnGrid(2,OceanAdaptiveGrid%nLat,OceanAdaptiveGrid%nPhi))
	velocityPsiTwoOnGrid=zero

	if(.NOT.allocated(velocityTwoOnGrid))allocate(velocityTwoOnGrid(2,OceanAdaptiveGrid%nLat,OceanAdaptiveGrid%nPhi))
	velocityTwoOnGrid=zero

	if (.not.allocated(VelocityMagnitude_PsiTwo)) allocate(VelocityMagnitude_PsiTwo(OceanAdaptiveGrid%nLat, OceanAdaptiveGrid%nPhi))

	if (.not.allocated(VelocityMagnitude_TotTwo)) allocate(VelocityMagnitude_TotTwo(OceanAdaptiveGrid%nLat, OceanAdaptiveGrid%nPhi))

	!ROCKET: uφ₁ (ƛ̂) =  (1 / (R·cosθ(i,j))) · (φ₁(i+1,j) − φ₁(i−1,j)) / (2Δλ) ƛ̂
	!ROCKET: uφ₁ (θ̂) = (1/R) · (φ₁(i,j+1) − φ₁(i,j−1)) / (2Δθ) θ̂ 

	VelocityMagnitude_PsiTwo = zero
	VelocityMagnitude_TotTwo = zero

	nThreads = OMP_get_max_threads()
	call OMP_SET_NUM_THREADS(nThreads)
	CHUNK = OceanAdaptiveGrid%nLat/nThreads

	write(*,*)CHUNK

	!$OMP PARALLEL DEFAULT(PRIVATE) &
	!$OMP& SHARED(CHUNK,OceanAdaptiveGrid,&
	!$OMP& velocityTwoOnGrid,velocityPsiTwoOnGrid,&
	!$OMP& psi2onGrid,VelocityMagnitude_PsiTwo,VelocityMagnitude_TotTwo)

	threadID = omp_get_thread_num()
	i1 = threadID*CHUNK+1
	i2=i1-1+CHUNK
	!		write(*,*)threadID,i1,i2
	if(threadID+1.eq.nThreads)i2=OceanAdaptiveGrid%nLat
	!		do iLat=1,OceanAdaptiveGrid%nLat

	do iLat=i1,i2
		cosLat=OceanAdaptiveGrid%cosLatC(iLat);
		! write(*,*) 'COSINE = ', cosLat

		sinLat=OceanAdaptiveGrid%sinLatC(iLat);

		iLatP1=iLat+1; iLatM1=iLat-1
		if(iLat.eq.1)iLatM1=1
		if(iLat.eq.OceanAdaptiveGrid%nLat)iLatP1=OceanAdaptiveGrid%nLat

		do iLon=1,OceanAdaptiveGrid%nPhi
			iLonP1=iLon+1;iLonM1=iLon-1
			if(iLon.eq.1)iLonM1=OceanAdaptiveGrid%nPhi
			if(iLon.eq.OceanAdaptiveGrid%nPhi)iLonP1=1

			! !!u_psi -checked
			velocityPsiTwoOnGrid(1,iLat,iLon)=(psi2onGrid(iLat,iLonP1)-psi2onGrid(iLat,iLonM1))/(two* OceanAdaptiveGrid%dphi*R*cosLat) 

			velocityPsiTwoOnGrid(2,iLat,iLon)= -(psi2onGrid(iLatP1,iLon)-psi2onGrid(iLatM1,iLon))/(two* OceanAdaptiveGrid%dLat*R)

			VelocityMagnitude_PsiTwo(iLat,iLon) = sqrt( velocityPsiTwoOnGrid(1,iLat,iLon)**2 + velocityPsiTwoOnGrid(2,iLat,iLon)**2 )

			!!Total velocity: u = u_phi + u_psi
			velocityTwoOnGrid(1,iLat,iLon)=velocityPsiTwoOnGrid(1,iLat,iLon) 

			velocityTwoOnGrid(2,iLat,iLon)=velocityPsiTwoOnGrid(2,iLat,iLon)

			VelocityMagnitude_TotTwo(iLat,iLon) = sqrt(velocityTwoOnGrid(1,iLat,iLon)**2 + velocityTwoOnGrid(2,iLat,iLon)**2 ) 
		end do
	end do
	!$OMP END PARALLEL
	
	write(*,*) '* Velocity Field computations are done *'
	write(*,*)
	! write(*,*) 'Velocity Magnitude Total Max, Velocity Magnitude Total Min = ',maxval(VelocityMagnitude_TotTwo), minval(VelocityMagnitude_TotTwo) 

	!!Psi
	write(*,*) 'Velocity Magnitude Psi Max, Velocity Magnitude Psi Min = ',maxval(VelocityMagnitude_PsiTwo), minval(VelocityMagnitude_PsiTwo) 
	
end subroutine VelocityFieldsTwoOnGrid

!TOFIX: (I do not advise using this subroutine yet still in progress and i was experimenting so...)
!TOFIX: Subroutine designed for getting the values of ζ : relative vorticity.
subroutine Compute_RelativeVorticity_Zeta()
	use OceanGridDataModule
	use sparseMatricesDataModule
	use basicDataStructures
	use simulationParameters
	use lsqrModule
	use GMRESSettings
	use bathymetryDataModule
	!$ use omp_lib
	implicit none
	
	integer :: k, iLat, iLon
	real (kind=QR_K) :: dlat,dphi,dphi2, dLat2, cosLat, sinLat, sinLatP1, sinLatM1,latPHalf,latMHalf, lat, coriolisTerm,R2
	real (kind=QR_K) :: cosLatMHalf, cosLatPHalf, cosLat2, sinLat2, meow
	integer :: i1, i2, CHUNK, nThreads, threadID
	integer :: iLatP1, iLatM1, iLonP1, iLonM1

	!!Computing velocity using finite difference of phi

	!! For Relative velocity Comparisons
	if (.NOT.allocated(zetaCurlVelocity)) allocate(zetaCurlVelocity(OceanGrid%nLat,OceanGrid%nPhi))
	zetaCurlVelocity = zero
	if (.NOT.allocated(zetaCurlVelocityPsi)) allocate(zetaCurlVelocityPsi(OceanGrid%nLat,OceanGrid%nPhi))
	zetaCurlVelocityPsi = zero
	if (.NOT.allocated(zetaLaplacianPsi)) allocate(zetaLaplacianPsi(OceanGrid%nLat,OceanGrid%nPhi))
	zetaLaplacianPsi = zero

	if (.NOT.allocated(zetaMinusF)) allocate(zetaMinusF(OceanGrid%nLat,OceanGrid%nPhi))
	zetaMinusF = zero

	R2 = R*R
	dlat = OceanGrid%dLat
	dphi = OceanGrid%dphi
	dLat2 = dlat*dlat
	dphi2 = dphi*dphi

	do iLat = 1, OceanGrid%nLat
		lat= -half*pi+half*OceanGrid%dLat+OceanGrid%dLat*real(iLat-1,QR_K)

		!conditions
		iLatP1 = iLat + 1; iLatM1 = iLat - 1
		if(iLatP1.gt.OceanGrid%nLat) iLatP1 = iLat 
		if(iLatM1.lt.1) iLatM1 = iLat
	
		cosLat=OceanGrid%cosLatC(iLat); !!cos^2 theta_j
		cosLat2=cosLat*cosLat

		!*Conditions on Latitude with 1/2: cos (theta+1/2) = cos(+ -pi/2) =0 automatically
		!========================================================
		latPHalf = Lat+half*dLat
		latMHalf = Lat -half*dLat
		
		cosLatMHalf=cos(latMHalf) !!cos(theta_{i,j-1/2})
		cosLatPHalf=cos(latPHalf) !!cos(theta_{i,j+1/2})

		if(iLat.eq.1) cosLatMHalf=zero
		if(iLat.eq.OceanGrid%nLat) cosLatPHalf=zero
		!========================================================

		sinLat=OceanGrid%sinLatC(iLat) !!sin theta_j
		coriolisTerm = two*OMEGA*sinLat 
		
		do iLon=1,OceanGrid%nPhi
			!conditions
			iLonP1=iLon+1;iLonM1=iLon-1
			if(iLon.eq.1)iLonM1=OceanGrid%nPhi
			if(iLon.eq.OceanGrid%nPhi)iLonP1=1

			!*zeta from laplacian of psi
			zetaLaplacianPsi(iLat,iLon) = (1/R2*cosLat2)*((psi1onGrid(iLat,iLonP1) -2*psi1onGrid(iLat,iLon) +psi1onGrid(iLat,iLonM1))/dphi2) + (1/R2*cosLat)*( (cosLatPHalf*(psi1onGrid(iLatP1,iLon) -psi1onGrid(iLat,iLon)) - cosLatMHalf*(psi1onGrid(iLat,iLon) - psi1onGrid(iLatM1,iLon)) ) /dLat2)

			!*zeta from curl of velocity
			zetaCurlVelocityPsi(iLat,iLon) = (1/R*cosLat)*( (velocityPsiOneOnGrid(1,iLat,iLonP1) - velocityPsiOneOnGrid(1,iLat,iLonM1))/(two*dphi)) - (1/R*cosLat)*((cosLatPHalf*(velocityPsiOneOnGrid(2,iLatP1,iLon) - velocityPsiOneOnGrid(2,iLat,iLon))- (cosLatMHalf*(velocityPsiOneOnGrid(2,iLat,iLon)-velocityPsiOneOnGrid(2,iLatM1,iLon)))) /(two*dlat))

			!* zeta from curl of velocity
			zetaCurlVelocity(iLat,iLon) = (1/R*cosLat)*( (velocityOnGrid(1,iLat,iLonP1) - velocityOnGrid(1,iLat,iLonM1))/(two*dphi)) - (1/R*cosLat)*((cosLatPHalf*(velocityOnGrid(2,iLatP1,iLon) - velocityOnGrid(2,iLat,iLon))- (cosLatMHalf*(velocityOnGrid(2,iLat,iLon)-velocityOnGrid(2,iLatM1,iLon)))) /(two*dlat))

			zetaMinusF(iLat,iLon) = - coriolisTerm

		end do
	end do
	write(*,*) '------------------------------------------------------------------------------------------------------------------------'
	write(*,*) 'Relative vorticity from Laplacian Psi1 : Maximum value , Minimum value = ', maxval(zetaLaplacianPsi), minval(zetaLaplacianPsi)
	write(*,*) 'Relative vorticity from Curl of velocity u_tot_1 : Maximum value , Minimum value = ', maxval(zetaCurlVelocity), minval(zetaCurlVelocity)

	write(*,*) 'Relative vorticity from Curl of velocity u_psi_1 : Maximum value , Minimum value = ', maxval(zetaCurlVelocityPsi), minval(zetaCurlVelocityPsi)

	write(*,*) 'Relative vorticity =? -f : Maximum value , Minimum value = ', maxval(zetaMinusF), minval(zetaMinusF)
	write(*,*) 

end subroutine Compute_RelativeVorticity_Zeta

end module poissonSolverModule

