module solve

  use basic
  use def_types
  use tools
  implicit none
  
contains
  
  !***********************************************************************
  ! SOLVEDP
  !***********************************************************************
  subroutine solvedp(v,x,pstar, f,n,m,P,beta,alg,v_ini)
  !***********************************************************************
  ! SOLVEDP solves a dynamic programming problem with optimal stopping
  ! Usage: 
  !	call solvedp(v,x,pstar, f,n,m,P,beta,alg,v_ini)
    
  ! INPUTS
  !     f          : reward function (n X m)
  !     n          : # of states 
  !     m          : # of actions
  !     beta       : discount factor
  !     alg        : algorithm to use (1:policy & 2:value)
  !     v_ini      : initial value for value function 
    
  ! OUTPUTS
  !     v          : value function   
  !     x          : policy function
  !     pstar      : optimal transition matrix
  !************************************************************************

    real :: start, finish 
    
    !Local 
    integer :: info, ipiv(n)
    integer, parameter :: maxit = 5000
    real, parameter :: tol = 10e-6
    logical, parameter :: prtiters = .true.

    integer :: i, j, it
    real :: vold(n), fstar(n)
    real :: change
    real :: a(n,n)
    
    !Dummy 
    real, intent(out) :: v(n), x(n)
    real, intent(out) :: pstar(n,n)

    integer, intent(in) :: n, m, alg
    real, intent(in) :: f(n,m), beta, v_ini

    type(sparse_matrix), intent(in) :: P

    !************ Perform policy or value function iterations
    
    v = reshape1(v_ini, n)   
    
    !************
    ! alg = 1 :: Policy
    ! alg = 2 :: Value
    !************
    
    select case (alg)
    case(1) ! Policy
       if (prtiters)then
          write(*,*) 'Solving Bellman by policy function iteration'
       end if
       do it = 1,maxit
	  vold  = v
	  call valmax(v,x, vold,f,P,beta, n,m)
!  	  call print_vector(x,100)
	  call valpol(pstar,fstar, x,f,P,beta, n,m)
	  !**************************************
	  !* Solving with LAPACK
	  !* (eye(n)-beta*pstar)*v=fstar	 
 	  a = (eye(n)-beta*pstar) 
 	  call sgetrf ( n, n, a, n, ipiv, info )
 	  call sgetrs ( 'n', n, 1, a, n, ipiv, fstar, n, info )
!  	  call print_matrix(a(1:100,1:8),100,8)
!  	  write (*,*) shape(a) 
 	  v = fstar
!    	  call print_vector(v(1:100),100)
	  !**************************************
 	  change = norm2(v-vold)
	  if (prtiters)then
              write(*,*) it, ' ' , change 
	  end if
	  if (change<tol)then
             exit
	  end if
       end do
    case(2) ! Value
       if (prtiters)then
          write(*,*) 'Solving Bellman by value function iteration'
       end if
       do it = 1,maxit
          vold = v
          call valmax(v,x, vold,f,P,beta, n,m)
          change = norm2(v-vold)
          if (prtiters)then
             write(*,*) it, ' ' , change 
          end if
          if (change<tol)then
             exit
          end if
       end do
       if (change>tol)then
          write(*,*) 'Failure to converge in solvedp'
       end if
    end select
  end subroutine solvedp
  
  !***********************************************************************
  ! VALMAX
  !***********************************************************************
  subroutine valmax(v,x, vold,f,P,beta, n,m)
        
    !Local
    real ::  h(n,m)
    !Dummy
    real, intent(out) :: v(n), x(n)
    
    integer, intent(in) ::  n, m
    real, intent(in) :: beta, f(n,m), vold(n)
    
    type(sparse_matrix) :: P
    h = reshape(sparse_matmul(P,size(P%values),vold,n),(/n,m/))
!     call print_matrix(f,100,8)
    v = maxval(f + beta * reshape(sparse_matmul(P,size(P%values),vold,n), (/n, m/)),2)
    x = maxloc(f + beta * reshape(sparse_matmul(P,size(P%values),vold,n), (/n, m/)),2)

  end subroutine valmax
  
  !***********************************************************************
  ! VALPOL
  !***********************************************************************
  subroutine valpol(pstar,fstar, x,f,P,beta, n,m)

     implicit none
    
    !Local
    integer :: j
    real :: a(n*m,n)
    real :: ff(n*m,1)
    !Dummy
    real, intent(out) :: pstar(n,n),fstar(n)

    integer, intent(in) :: n, m
    type(sparse_matrix),intent(in) :: P

    real, intent(in) :: x(n), f(n,m), beta
    
    a = sparsetomatrix(P,size(P%d1),n*m,n)
    ff = reshape(f,(/n*m,1/))
!     call print_matrix(a(1,1),1,1)
    do j = 1,n	  
      fstar(j) = ff(n*(int(x(j))-1)+j,1)
      pstar(j,:) = a(n*(int(x(j))-1)+j,:)
    end do
!      call print_matrix((/P%values(1:100),real(P%d1(1:100)),real(P%d2(1:100))/),100,3)
  end subroutine valpol
end module solve
