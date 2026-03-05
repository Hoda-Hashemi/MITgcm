!> \todo
!! - CGNE_NONNEG must include GMRES_RESTARTED, RRGMRES, LSQR, CGLS, CGNR, CGNE with the Armijo condition
!! satisfied
!! - Is it possible to implement NL-BCG ?
!!
!! \defgroup OptimizationTools Optimization Tools
!!
!! \defgroup NES Solvers for Normal Equations
!! \ingroup OptimizationTools
!! \brief This module contains a collection of subroutines to solve the problem A'A x = A' b where A is an m by n matrix
!! \author Issam Lakkis
!!  \version 1.0
!!  \date    August 2009
!!
!! \defgroup SPD Solvers for Symmetric Positive Definite Linear Systems
!! \ingroup OptimizationTools
!! \brief This module contains a collection of subroutines to solve the problem A x = b where A is a symmetric positive definite square matrix
!! \author Issam Lakkis
!!  \version 1.0
!!  \date    August 2009
!!
!! \defgroup PD Solvers for Positive Definite Linear Systems
!! \ingroup OptimizationTools
!! \brief This module contains a collection of subroutines to solve the problem A x = b where A is a positive definite square matrix
!! \author Issam Lakkis
!!  \version 1.0
!!  \date    August 2009
!!
!! \defgroup SLS Solvers for Square Linear Systems
!! \ingroup OptimizationTools
!! \brief This module contains a collection of subroutines to solve the problem A x = b where A is a square matrix
!! \author Issam Lakkis
!!  \version 1.0
!!  \date    August 2009
!!
!! \defgroup NL Non-Linear Solvers Satisfying sign condition
!! \ingroup OptimizationTools
!! \brief This module contains a collection of subroutines to solve the problem A'A x = A' b where A is an m by n matrix with the condition
!!  x(i) b(i) >= 0 by assuming x(i) = sign(b(i)) g(x), where g(x) is a non-negative function of x. The solution is obtained by minimizing a
!! nonlinear function
!! \author Issam Lakkis
!!  \version 1.0
!!  \date    August 2009
!!
!! \defgroup LSIGN Linear Solvers Satisfying sign condition
!! \ingroup OptimizationTools
!! \brief This module contains a collection of subroutines to solve the problem A'A x = A' b where A is an m by n matrix with the condition
!!  x(i) b(i) >= 0. This is done by invoking an outer loop that projects the solution obtained from the inner loop to satisfy the sign condition.
!! \author Issam Lakkis
!! \version 1.0
!! \date    August 2009
!<

!Y!art 	Algebraic reconstruction technique (Kaczmarz’s method)
!Y!cgls 	Computes the least squares solution based on k steps of the conjugate gradient algorithm
!Y!lsqrb 	Computes the least squares solution based on k steps of the LSQR algorithm
!Y!rrgmres Range-restricted GMRES algorithm for square systems
!N!maxent Computes the maximum entropy regularized solution
!N!mr2 	Solution of symmetric indeﬁnite problems by MR-II
!N!nu 	Computes the solution based on k steps of Brakhage’s iterative ν -method
!N!splsqr Computes an approximate standard-form Tikhonov solution via the subspace preconditioned LSQR algorithm (SP-LSQR)

!Y! Conjugate Gradient
!Y! BiConjugate Gradient
!Y! GMRES
!N! SIMPLEX

!pcgls 		Same as cgls, but for general-form regularization
!plsqrb 		Same as lsqrb, but for general-form regularization
!pmr2 		Same as mr2, but for general-form regularization
!pnu 		Same as nu, but for general-form regularization
!prrgmres	Same as rrgmres, but for general-form regularization

module lsqrModule
  implicit none
! solve A' A x = A' b in the least square sense - for inverse problems

  public :: BICGSTAB_OMP
  public :: GMRES_RESTARTED_OMP_RETURNR

  private :: dot
  private :: dotOMPNEW
  private :: AxOMP

contains

  !#########################################################################################################################################
  !--------------------------------------------------------------------------------
  !> Bi-Conjugate Gradient Stabilized Method for solving  \f$ {\mathbf A} \vec{x} = \vec{b} \f$. Matrix  \f$ {\mathbf A} \f$ is \b square \b non-symmetric.
  !!
  !! The CGS algorithm is based on squaring the residual polynomial,and ,in cases of irregular
  !! convergence, this may lead to substantial build-up of rounding errors, or possibly even
  !! overﬂow. The Biconjugate Gradient Stabilized (BICGSTAB) algorithm is a variation of
  !! CGS which was developed to remedy this difficulty. \n
  !!
  !! Iterations stop if one of the following two conditions is satisfied: \n
  !! (1) maximim number of iterations (niter) is exceeded or \n
  !! (2) convergence is reached according to : \f$||{\mathbf A} \vec{x} - \vec{b}||_2 <= \mbox{res} \, ||\vec{b}||_2 \f$
  !!
  !! \b Reference: \n
  !! - Iterative Methods for Sparse Linear Equations by Y. Saad, 2000, Algorithm 7.6 Page 219. \n
  !! - H. A. Van der Vorst. Bi-CGSTAB. A Fast and Smoothly Converging varient of Bi-CG for the Solution of non-Symmetric Linear Systems.
  !! SIAM Journal on Scientific and Statistical Compution, 12:631-644, 1992.
  !!
  !! \ingroup SLS
  !!
  !<
  subroutine BICGSTAB_OMP(A,b,x,n,nIter,res)
  use sparseMatricesDataModule
  implicit none
  type(sprsmatrixA), intent(in) :: A !< coefficients matrix of size ( n x n), stored in sparse format
  integer,  intent(in)  :: n !< number of unknowns (number of columns of matrix A)
  integer,  intent(in)  :: nIter !< maximum number of iterations
  double precision, intent(in) :: b(n) !< vector of size n
  double precision, intent(inout) :: x(n) !< solution vector of size n
  double precision, intent(in) :: res !< parameter for convergence. See above.
  integer :: k
  double precision :: r(n),r0(n),p(n),s(n),v(n),As(n),rho,rhoPrev,dp0,dp1,dp2
  double precision :: normb,alpha,beta,omega,resnorm,regularityMin

  call AxOMP(A,x,r0,0); r0 = b - r0; r = r0 !r0 = b - A x0
  !$OMP WORKSHARE
  rhoPrev = dot_product(r,r0)
  !$OMP END WORKSHARE
  resnorm = DSQRT(DABS(rhoPrev))
  p = r0
  !$OMP WORKSHARE
  normb = dot_product(b,b)
  !$OMP END WORKSHARE
  normb = DSQRT(normb)

!write(*,*)'xxxxx',resnorm,normb

  k = 1
  regularityMin = 1.0D0
  do while(k.lt.nIter.and.resnorm.gt.res*normb)
  	k=k+1
  	call AxOMP(A,p,v,0);
    !$OMP WORKSHARE
    dp0 = dot_product(v,r0)
    !$OMP END WORKSHARE

    alpha = rhoPrev/dp0
  	s = r - alpha*v
  	call AxOMP(A,s,As,0);
    !$OMP WORKSHARE
    dp1 = dot_product(As,s)
    !$OMP END WORKSHARE
    !$OMP WORKSHARE
    dp2 = dot_product(As,As)
    !$OMP END WORKSHARE

    omega = dp1/dp2
  	x = x + alpha*p + omega*s
  	r = s - omega*As
    !$OMP WORKSHARE
    rho = dot_product(r,r0)
    !$OMP END WORKSHARE

  	beta = (rho/rhoPrev)*(alpha/omega); rhoPrev = rho
  	p = r + beta*(p - omega*v)
  	resnorm = DSQRT(DABS(rho))
    write(*,*)k,nIter,resnorm,res*normb
  end do
  call AxOMP(A,x,r,0); r = b - r;
  !$OMP WORKSHARE
  dp1 = dot_product(r,r)
  !$OMP END WORKSHARE
  resnorm = DSQRT(dp1)
  write(*,*)'BICGSTAB target residual=',k,nIter,resnorm,res*normb
 
  return
  end subroutine BICGSTAB_OMP

  ! !ALTERED !ADDED residuals
  ! subroutine BICGSTAB_OMP_HODA(A,b,x,n,nIter,res)
  !   use sparseMatricesDataModule
  !   implicit none
  !   type(sprsmatrixA), intent(in) :: A !< coefficients matrix of size ( n x n), stored in sparse format
  !   integer,  intent(in)  :: n !< number of unknowns (number of columns of matrix A)
  !   integer,  intent(in)  :: nIter !< maximum number of iterations
  !   double precision, intent(in) :: b(n) !< vector of size n
  !   double precision, intent(inout) :: x(n) !< solution vector of size n
  !   double precision, intent(in) :: res !< parameter for convergence. See above.
  !   integer :: k
  !   double precision :: r(n),r0(n),p(n),s(n),v(n),As(n),rho,rhoPrev,dp0,dp1,dp2
  !   double precision :: normb,alpha,beta,omega,resnorm,regularityMin
  
  !   !ADDED: variables for infinity-norm residual
  !   double precision :: r_inf, b_inf
  
  !   call AxOMP(A,x,r0,0); r0 = b - r0; r = r0 !r0 = b - A x0
  !   !$OMP WORKSHARE
  !   rhoPrev = dot_product(r,r0)
  !   !$OMP END WORKSHARE
  !   resnorm = DSQRT(DABS(rhoPrev))
  !   p = r0
  !   !$OMP WORKSHARE
  !   normb = dot_product(b,b)
  !   !$OMP END WORKSHARE
  !   normb = DSQRT(normb)
  
  !   !ADDED: compute b_inf once (∞-norm of RHS)
  !   b_inf = maxval(abs(b))
  !   if (b_inf == 0.0d0) then 
  !     write(*,*) '**** b_inf == 0.0d0 ' 
  !     b_inf = 1.0d0
  !   end if
  
  ! !write(*,*)'xxxxx',resnorm,normb
  
  !   k = 1
  !   regularityMin = 1.0D0
  !   do while(k.lt.nIter.and.resnorm.gt.res*normb)
  !     k=k+1
  !     call AxOMP(A,p,v,0);
  !     !$OMP WORKSHARE
  !     dp0 = dot_product(v,r0)
  !     !$OMP END WORKSHARE
  
  !     alpha = rhoPrev/dp0
  !     s = r - alpha*v
  !     call AxOMP(A,s,As,0);
  !     !$OMP WORKSHARE
  !     dp1 = dot_product(As,s)
  !     !$OMP END WORKSHARE
  !     !$OMP WORKSHARE
  !     dp2 = dot_product(As,As)
  !     !$OMP END WORKSHARE
  
  !     omega = dp1/dp2
  !     x = x + alpha*p + omega*s
  !     r = s - omega*As
  !     !$OMP WORKSHARE
  !     rho = dot_product(r,r0)
  !     !$OMP END WORKSHARE
  
  !     beta = (rho/rhoPrev)*(alpha/omega); rhoPrev = rho
  !     p = r + beta*(p - omega*v)
  !     resnorm = DSQRT(DABS(rho))
  
  !     !ADDED: compute infinity-norm of current residual
  !     call AxOMP(A, x, r, 0)     ! r = A * x
  !     r = b - r                  ! r = b - A * x
  !     r_inf = maxval(abs(r)) / b_inf
  
  !     write(*,*)k, nIter, resnorm, res*normb, 'inf_norm=', r_inf !MODIFIED: added r_inf
  !   end do
  
  !   call AxOMP(A,x,r,0); r = b - r;
  !   !$OMP WORKSHARE
  !   dp1 = dot_product(r,r)
  !   !$OMP END WORKSHARE
  !   resnorm = DSQRT(dp1)
  
  !   !ADDED: final infinity-norm residual
  !   r_inf = maxval(abs(r)) / b_inf
  
  !   write(*,*)'BICGSTAB target residual=',k,nIter,resnorm,res*normb, 'inf_norm=', r_inf !MODIFIED
  
  !   return residuals
  ! end subroutine BICGSTAB_OMP_HODA
  
  !--------------------------------------------------------------------------------
  !> Restarted GMRES Method for solving  \f$ {\mathbf A} \vec{x} = \vec{b} \f$. Matrix  \f$ {\mathbf A} \f$ is \b Square.
  !!
  !! The algorithm utilizes the Modiﬁed Gram-Schmidt orthogonalization in the Arnoldi process. \n
  !! The algorithm provides the option to enforce non-negativity by setting signFlag = 1. It simply projects
  !! solution onto non-negative space in the outer iterations as proposed in reference 2.
  !! Iterations stop if one of the following two conditions is satisfied: \n
  !! (1) maximim number of iterations (niter) is exceeded or \n
  !! (2) convergence is reached according to : \f$||{\mathbf A} \vec{x} - \vec{b}||_2 <= \mbox{res} \, ||\vec{b}||_2 \f$
  !!
  !! \b Reference: \n
  !! - Iterative Methods for Sparse Linear Equations by Y. Saad, 2000, Algorithm 6.11 Page 167. \n
  !! - Calvetti et al. Non-negativity and iterative methods for ill-posed problems. Inverse Problems 20 (2004) 1747-1758
  !!
  !! \ingroup SLS
  !!
  !<
  subroutine GMRES_RESTARTED_OMP_RETURNR(A,b,x,signConstrained,r,n,mMax,nIter,res,success,normr)
  use sparseMatricesDataModule
  !$ use omp_lib
  implicit none

  type(sprsmatrixA), intent(in) :: A !< coefficients matrix of size ( n x n), stored in sparse format
  integer,  intent(in)  :: n !< number of unknowns (number of columns of matrix A)
  integer,  intent(in)  :: nIter !< maximum number of iterations
  integer,  intent(in)  :: mMax !< maximum number of (inner) iterations to take.  0 < m <= m. Typically m is a small integer
  logical, intent(in) :: signConstrained
  double precision, intent(inout) :: b(n) !< vector of size n
  double precision, intent(out) :: r(n) !< vector of size n
  double precision, intent(inout) :: x(n) !< solution vector of size n
  double precision, intent(in) :: res !< parameter for convergence. See above.
  double precision, intent(inout) :: normr !< is used as cutoff circulation when in
  double precision :: w(n),x1(n),xOpt(n),rOpt(n)
  double precision :: hi,hiP1,dotproduct
  double precision :: d,gj,gjP1,resIter,normb,normrPrev,normmin,normb0
  integer :: i,j,iter,totiter,success,totitermax,m,leave,mm
  double precision :: wtimeTot,cutoffCirculation,xMax,cutoffFraction
!  double precision gg(n)

  !double precision, allocatable :: H(mm+1,mm),v(n,mm+1),y(mm),g(mm+1),hj(mm),hjP1(mm),s(mm),c(mm)
  double precision, allocatable :: H(:,:),v(:,:),y(:),g(:),hj(:),hjP1(:),s(:),c(:)
  double precision :: normminGLOBAL,xOptGlobal(n),rOptGlobal(n)
  integer :: resetFlag=1
  !cutoffCirculation=normr

  cutoffFraction=1.0D-3
  cutoffCirculation=-1.0D0
  wtimeTot = 0.0D0

  !normb0 = DSQRT(dot(b,b,n))
  !$OMP WORKSHARE
  normb0 = dot_product(b,b)
  !$OMP END WORKSHARE
  normb0 = DSQRT(normb0)

  b(1:n)=b(1:n)/normb0
  normb = 1.0D0
  mm=5
  !mm=mMax
  normminGLOBAL = 1.0D30

  do while(mm.le.mMax)
    allocate(H(mm+1,mm))
    allocate(v(n,mm+1))
    allocate(y(mm))
    allocate(g(mm+1))
    allocate(hj(mm))
    allocate(hjP1(mm))
    allocate(s(mm))
    allocate(c(mm))
    normmin=1.0D30
    m=mm

    resIter = 1.0D6
    iter = 0; totiter = 0;

    if(resetFlag.eq.1.or.mm.eq.5)then
      x(1:n) = 0.0D0
      r=b
      !normr = DSQRT(dot(r,r,n))
      !$OMP WORKSHARE
      normr = dot_product(r,r)
      !$OMP END WORKSHARE
      normr = DSQRT(normr)
    else
      x = xOptGlobal
      r = rOptGlobal
      normr=normminGLOBAL
    end if

    normrPrev = normr
    totiterMax = nIter*mm
    leave=0

    !write(*,*)'starting with',normr,normb

    !write(*,*)'##################################################################'
    do while(normmin.gt.res*normb.and.leave.eq.0)
      iter = iter+1
      v(1:n,1)=r/normr; H(1:m+1,1:m) = 0.0D0
      g(1:m+1) = 0; g(1) = normr
      !j = 0 ;
      do j=1,m
        totiter = totiter+1

        !		wtime = omp_get_wtime()
        call AxOMP(A,v(1:n,j),w,0)
        !call AxA(A,v(1:n,j),w,0)
        !		wtime = omp_get_wtime() - wtime
        !		wtimeTot = wtimeTot+wtime

        !do i=max(1,j-k+1),j   !for incomplete orthogonalization
        do i=1,j
          !H(i,j) = dot(w,v(1:n,i),n)
          !$OMP WORKSHARE
          H(i,j) = dot_product(w,v(1:n,i))
          !$OMP END WORKSHARE
          w = w - H(i,j)*v(1:n,i)
        end do

        !H(j+1,j)=DSQRT(dot(w,w,n))
        !$OMP WORKSHARE
        dotproduct = dot_product(w,w)
        !$OMP END WORKSHARE
        H(j+1,j) = DSQRT(dotproduct)

        if(H(j+1,j).eq.0.0D0) then
          m=j; exit
        end if
        v(1:n,j+1) = w/H(j+1,j)

        ! apply j-1 Givens rotations to jth column of H
        if(j.gt.1) then
          do i=1,j-1
            hi = c(i)*H(i,j)+s(i)*H(i+1,j); hiP1 = -s(i)*H(i,j)+c(i)*H(i+1,j)
            H(i,j) = hi; H(i+1,j) = hiP1
          end do
        end if

        ! apply jth Givens rotation to H and g
        d = DSQRT(H(j,j)**2+H(j+1,j)**2);s(j) = H(j+1,j)/d;c(j) = H(j,j)/d
        !new jth and j+1 rows of H after multiplying it form the left with a Givens rotation
        hj(1:m) = c(j)*H(j,1:m)+s(j)*H(j+1,1:m); hjP1(1:m) = -s(j)*H(j,1:m)+c(j)*H(j+1,1:m)
        H(j,1:m) = hj(1:m); H(j+1,1:m) = hjP1(1:m)
        ! applying m rotations to g
        gj =  c(j)*g(j)+s(j)*g(j+1); gjP1 = -s(j)*g(j)+c(j)*g(j+1); g(j)=gj; g(j+1)=gjP1
        resIter = DABS(gjP1)
      end do

      ! getting y
      y(m) =  g(m)/H(m,m)

      do j=m-1,1,-1
        !dotproduct = dot(H(j,j+1:m),y(j+1:m),m-j)
        !$OMP WORKSHARE
        dotproduct = dot_product(H(j,j+1:m),y(j+1:m))
        !$OMP END WORKSHARE
        y(j)=(g(j)-dotproduct)/H(j,j)
      end do

      ! this is matrix multiplication
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
      !$OMP DO SCHEDULE(STATIC)
      do i=1,n
        w(i) =dot(v(i,1:m),y(1:m),m)
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      x1 = x+w

  if(signConstrained)then
  ! computing cutoffCirculation if not initialized (i.e. given a negative value)
      if(cutoffCirculation.le.0.0D0)then
        xMax=-1.0D30
        do i=1,n;
          if(DABS(x1(i)).gt.xMax)xMax=DABS(x1(i))
        end do;
        cutoffCirculation=cutoffFraction*xMax;
      end if

      do i=1,n
        !    if(x1(i)*b(i).lt.0.0D0)x1(i)=0.0D0;
        if(x1(i)*b(i).lt.0.0D0.and.DABS(x1(i)).gt.cutoffCirculation)then
          if(b(i).lt.0)x1(i)=cutoffCirculation;
          if(b(i).gt.0)x1(i)=-cutoffCirculation;
        end if
      end do;
  end if

      !	wtime = omp_get_wtime()
      call AxOMP(A,x1,r,0);
      !call AxA(A,x1,r,0);
      !	wtime = omp_get_wtime() - wtime
      !	wtimeTot = wtimeTot+wtime

      r=b-r;

      !normr=DSQRT(dot(r,r,n))
      !$OMP WORKSHARE
      normr = dot_product(r,r)
      !$OMP END WORKSHARE
      normr=DSQRT(normr)

      !write(*,*)'GMRES',m,normr,normminGLOBAL,res,(normr-normrPrev)/normrPrev,cutoffCirculation

      !	if(normr.gt.normrPrev.and.m.eq.1)leave=1
      if((normr-normrPrev)/normrPrev.gt.-0.01D0.and.m.eq.1)leave=1
      if(m.eq.1.and.iter.gt.100)leave=1

      if(normr.lt.normmin) then
        normmin=normr
        xOpt=x1
        rOpt=r
      end if
      if(normr.lt.normrPrev) then
        x=x1
      else
        x = xOpt
        r = rOpt
      end if

      !if(m.gt.5)then
      !  m=m-5
      !else
      if((DABS((normr-normrPrev)/normr).lt.0.1D0.or.normr.gt.normrPrev).and.m.gt.1)then
        if(m.gt.10)then
          m=m-5
        else
          m=m-1
        end if
        iter=0
      end if
      !end if
      !m=mm
      !write(*,*)m
      normrPrev = normr
      write(*,*)iter,normmin,res*normb,leave 

    end do

    write(*,*)
    write(*,*)'GMRES',m,mm,normr,normminGLOBAL,res, cutoffCirculation

    deallocate(H)
    deallocate(v)
    deallocate(y)
    deallocate(g)
    deallocate(hj)
    deallocate(hjP1)
    deallocate(s)
    deallocate(c)

    if(normmin.lt.normminGLOBAL)then
      normminGLOBAL=normmin
      xOptGlobal=xOpt
      rOptGlobal=rOpt
      mm=mm+5
    else
      mm=mMax+5
    end if

    if(normminGLOBAL.le.res)then
      mm=mMax+5
    end if
  end do

  !write(*,*)'##################################################################'
  !write(*,*)'Global normin=',normminGLOBAL
  !write(*,*)'##################################################################'

  x=xOptGlobal
  !r=rOptGlobal

  !wtime = omp_get_wtime()
  call AxOMP(A,x,r,0)
  !call AxA(A,x,r,0)
  !wtime = omp_get_wtime() - wtime
  !wtimeTot = wtimeTot+wtime

  r=b-r;
  x(1:n)=x(1:n)*normb0
  b(1:n)=b(1:n)*normb0
  r(1:n)=r(1:n)*normb0
  !normr=DSQRT(dot(r,r,n))
  !$OMP WORKSHARE
  normr = dot_product(r,r)
  !$OMP END WORKSHARE
  normr=DSQRT(normr)
  !write(*,*)normb0, normr, normr/normb0, res
  !write(*,*) 'sum |b|, normb0, normr , normr/normb0=',sum(DABS(b)), normb0, normr, normr/normb0
  !write(*,*)normr/normb0,res
  !write(*,*)totiter,totiterMax
  !write(*,*)leave
  !write(*,*)'time of Ax',wtimeTot

  return
  end subroutine GMRES_RESTARTED_OMP_RETURNR

  !----------------------
  double precision function dot(a,b,n)
  implicit none
  double precision :: b(*), a(*)
  integer :: n, i

  dot = 0.0D0
  do i=1,n
  	dot = dot+a(i)*b(i)
  end do
  return
  end function dot

  !----------------------------------------------------
  ! SUBROUTINE AxOMP: created on Nov. 7 2011: open mp version of subroutine Ax
  ! y = y + Ax
  !> Ax computes \f$ \vec{y} = \vec{y}_0 + {\mathbf A} \vec{x} \f$
  !!
  !! If flag = 0 then \f$ \vec{y}_0 = 0 \f$ else \f$ \vec{y}_0 = \vec{y} \f$
  !<
  SUBROUTINE AxOMP(A,x,y,flag)
  	use sparseMatricesDataModule
  	!$ use omp_lib
  	implicit none
  	type(sprsmatrixA), intent(in) :: A !< Sparse Matrix A in Sparse Storage Format
    double precision, intent(in)    :: x(*) !< vector x
  	double precision, intent(inout) :: y(A%m) !< vector y
  	integer, intent(in) :: flag !< see above
  	integer(KIND=8) :: inx,i
  	double precision :: t
  	integer(KIND=8)  :: nentries(A%m+1)
  !	double precision :: wtime

  !wtime = omp_get_wtime()

  if(flag.eq.0)y(1:A%m) = 0.0D0

  !!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(A,y,x)
  !!$OMP DO SCHEDULE(STATIC)
  !	do inx = 1, A%length
  !		!OMP ATOMIC
  !		y(A%i(inx)) = y(A%i(inx)) + A%value(inx)*x(A%j(inx))
  !	end do
  !!$OMP END DO
  !!$OMP END PARALLEL

  nentries(1) = 0
  do i=2,A%m+1
  	nentries(i) = A%nentries(i-1)
  end do

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,inx,t)
  !$OMP DO SCHEDULE(STATIC)
  do i=1,A%m
  	t = 0.0D0
  	do inx=nentries(i)+1,nentries(i+1)
  		t = t + A%value(inx)*x(A%j(inx))
  	end do
  	y(i) = t
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !wtime = omp_get_wtime() - wtime
  !write(*,*)'Ax time in AxOMP :', wtime
  	return
  END SUBROUTINE AxOMP

  double precision function dotOMPNEW(a,b,n,nThreads)
  !$ use omp_lib
  implicit none
  double precision :: b(n), a(n)
  integer :: n, i
  double precision :: S(nThreads)
	integer :: i1,i2,CHUNK,nThreads,threadID

  CHUNK = n/nThreads
!  write(*,*)CHUNK
  S=0.0D0

  !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(n,a,b,S)
  threadID = omp_get_thread_num()
  i1 = threadID*CHUNK+1
  i2=i1-1+CHUNK
  if(threadID+1.eq.nThreads)i2=n
  do i=i1,i2
  	S(threadID+1) = S(threadID+1)+a(i)*b(i)
  end do
  !$OMP END PARALLEL

  dotOMPNEW = sum(S)

  return
end function dotOMPNEW

end module LSQRmodule

