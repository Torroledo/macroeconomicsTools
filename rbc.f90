!********************************************************************************
!
! RBC.f90 Solves standard real business cycle model with inelastic labor supply.
!
! Written by Franz Hamann
! Adapated to Fortran by Ivan Torroledo 
!
!********************************************************************************

program rbc

  use basic
  use tools
  use solve
  use def_types
  
  implicit none
  
  !***********************************************************************
  ! Declarations
  !***********************************************************************
  real :: start, finish
  
  !******* Parameters
  
  real :: beta = 0.98, mu = 2., alpha = 0.3, delta = 0.1, sigma = 0.35, rho = 0.5
  integer :: iz, ik, i, j
  
  !******* Approximate z with nz discrete states AR(1)
  
  integer, parameter :: nz = 9
  real :: z(nz), prob(nz,nz)
  
  !******* Construct state space

  integer, parameter :: nk = 100
  real :: k(nk)
  
  real, dimension(nz*nk*nk,3) :: grid
  real, dimension(nz*nk*nk) :: z_m, k_m, kp_m

  integer ::  m = nk, n = nk*nz

  !******* Construct the return function, u

  real :: c(nz*nk*nk), u_p(nz*nk*nk), u(nz*nk,nk)
  
  !******* Construct the transition matrix

  type(sparse_matrix) :: P
  
  !******* Solve the DP problem using policy iteration

  real :: v_p(nk*nz), x(nk*nz), pstar(nz*nk,nz*nk), v(nz,nk)
  real :: d(nz*nk)
  
  !***********************************************************************
  ! Executable
  !***********************************************************************

  call cpu_time(start)	

  write (*,*) ' '
  write (*,*) 'Real business cycles model'
  write (*,*) ' '
  
  !******* Parameters

  !******* Approximate z with nz discrete states AR(1)

  call rouwenhorst(z,prob, nz,0.,rho,sigma)
  z = exp(z)
  
  !******* Construct state space

  k = linspace(0.1,20.,nk)
  grid = gridmake3(z,nz,k,nk,k,nk)
  
  z_m  = grid(:,1)
  k_m  = grid(:,2)
  kp_m = grid(:,3)

  m = nk
  n = nk*nz

  !******* Construct the return function, u
  
  c = z_m*k_m**alpha + (1-delta)*k_m - kp_m

  do i=1,n*m
    if (c(i)>=0)then 
      u_p(i) = (c(i)**(1-mu))/(1-mu)
    else
      u_p(i) = -999999999999999.
    end if
  end do 
  u = reshape(u_p, (/n,m/))
  !******* Construct the transition matrix

  call kroneckerS(P, eye(m),m, repmat(prob,nz,nz,m,1), nz*m,nz) 

  call solvedp(v_p,x,pstar, u,n,m,P,beta,1,0.)  
  call print_vector_dat('policy.dat',x,n)
  !d = ergdist(pstar,size(pstar,1));
  call cpu_time(finish)	
  write(*,*) 'time ' ,finish-start
end program rbc
