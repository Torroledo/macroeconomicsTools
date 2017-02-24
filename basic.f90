module basic
  
  implicit none

contains

  function repmat(mat,n_row,n_col)
  !**********************************************************
  ! REPMAT repets a matrix in a desired shape
  ! Usage: 
  !	matrix = repmat(mat, n_row, n_col)
  
  ! INPUTS
  !	mat	: real value matrix
  !	n_row	: # times to repeat mat in row direction 
  !	n_col	: # times to repeat mat in column direction 
  
  ! OUTPUTS	
  !	matrix	: real value matrix 
  !
  !***********************************************************
	integer :: i, j, k, l, n_row, n_col, dim(2)
    real :: mat(:,:)
    real, allocatable :: repmat(:,:)

	dim = shape(mat)
	
	allocate(repmat(dim(1)*n_row,dim(2)*n_col))

    repmat = 0.

	do i = 1,n_row
		do j = 1,n_col
		repmat((i-1)*dim(1)+1:dim(1)*i+1,(j-1)*dim(2)+1:dim(2)*j+1) = mat
		end do
	enddo
	
  end function repmat
  
  function eye(m)
  !**********************************************************
  ! EYE make a indentity m X m matrix 
  ! Usage: 
  !	identity = eye(m)
  
  ! INPUTS
  !	m	    : size of the square matrix 
  
  ! OUTPUTS	
  !	identity    : identity matrix  
  !
  !***********************************************************      
    integer :: i, m
    real :: eye(m,m)
    
    eye = 0. 
    do i = 1,m
       eye(i,i) = 1.
    end do
    
  end function eye
  
  !**********************************************************
  ! ONES
  !**********************************************************  
  function ones(m,n) result(matrix)
  !**********************************************************
  ! ONES make n x m matrix with ones in all entries   
  ! Usage: 
  !	result = ones(m,n)
  
  ! INPUTS
  !	m     	: # of rows 
  !	n		: # of columns 
    
  ! OUTPUTS	
  !	result     	: n x m matrix
  !
  !***********************************************************
        
    integer :: i, j, m, n
    real    :: matrix(m,n)
    
    matrix = 1.
    
  end function ones

  !**********************************************************
  ! LOAD DATA
  !**********************************************************  
  function loaddata(unit, namefile, rows,cols) result(data)
  !**********************************************************
  ! LOADDATA load numeric data  
  ! Usage: 
  !	data = loaddata(unit, namefile, rows,cols)
  
  ! INPUTS
  !	unit     	: integer indicator for each file 
  !	namefile	: string value with the name file
  !	rows	        : # of rows in the data 
  !	cols	        : # of columns in the data
    
  ! OUTPUTS	
  !	data     	: real value matrix 
  !
  !***********************************************************
    !Local
    
    !Dummy
    integer :: unit, rows, cols
    character (len=*) :: namefile

    real, dimension(rows,cols) :: data

    open(unit, file=namefile, action='read')
    read(unit,*) data
    close(unit)
    
  end function loaddata

  !**********************************************************
  ! LINSPACE
  !**********************************************************  
  function linspace(init, endd, points) result(line)
  !**********************************************************
  ! LINSPACE makes a row-vector equal spaciated   
  ! Usage: 
  !	line = linspace(init, end, points)
  
  ! INPUTS
  !	init    : initial value 
  !	end	: final value 
  !	points  : # of points in vector
    
  ! OUTPUTS	
  !	line    : row-vector with the data
  !
  !***********************************************************
    
    !Local
    integer :: i 
    
    !Dummy
    integer :: points
    real :: init, endd

    real :: line(points)
    
    do i = 1, points
       line(i) = init + (endd-init)*(i-1)/(points-1)
    end do
  end function linspace

  !**********************************************************
  ! RESHAPE1
  !**********************************************************  
  function reshape1(vini, size) result(v)
  !**********************************************************
  ! RESHAPE1 reshapes a scalar in a vector with scalar in the
  !          in the first position   
  ! Usage: 
  !	v = reshape1(vini, size)
    
  ! INPUTS
  !	vini     : scalar value 
  !	size	 : size of the vector
 
  ! OUTPUTS	
  !	data     : vector 
  !
  !***********************************************************
    
    !Local
    integer :: i

    !Dummy
    integer :: size
    real :: vini, v(size) 
    
    v = 0.
    v(1) = vini
    
  end function reshape1
  
  !**********************************************************
  ! PRINT VECTOR BASH
  !**********************************************************  
  subroutine print_vector(v)
  !**********************************************************
  ! PRINT_VECTOR prints a vector in the terminal
  ! Usage: 
  !	call print_vector(v)

  ! INPUTS
  !	v     	: vector to print 
  !  
  !***********************************************************
    !Local
    integer :: i
    
    !Dummy
    real, intent(in) :: v(:)

    do i = 1,size(v)
       write(*,*) v(i)
    end do
  end subroutine print_vector

  !**********************************************************
  ! PRINT MATRIX BASH
  !**********************************************************  
  subroutine print_matrix(v)
  !**********************************************************
  ! PRINT_MATRIX prints matrix in the terminal  
  ! Usage: 
  !	call print_matrix(v)
    
  ! INPUTS
  !	v     : matrix to print 
  !
  !***********************************************************
      
    !Local
    integer :: i
    
    !Dummy
    real, intent(in) :: v(:,:)

    do i = 1,size(v,1)
       write(*,*) v(i,:)
    end do
  end subroutine print_matrix

    !**********************************************************
  ! PRINT MATRIX DATA FILE (.dat)
  !**********************************************************  
  subroutine print_matrix_dat(namefile,v)
  !**********************************************************
  ! PRINT_MATRIX_DAT print matrix in a dat file   
  ! Usage: 
  !	call print_matrix_dat(namefile,v)
    
  ! INPUTS
  !     namefile : name of the file    
  !	v        : matrix to print 
  !
  !***********************************************************
    !Local
    integer :: i
    
    !Dummy
    real, intent(in) :: v(:,:)
    character (len=*) :: namefile

    open(1, file=namefile, action='write', status='replace')

    do i = 1,size(v,1)
       write(1,*) v(i,:)
    end do
    close(1)
  end subroutine print_matrix_dat

  !**********************************************************
  ! PRINT VECTOR DATA FILE (.dat)
  !**********************************************************  
  subroutine print_vector_dat(namefile,v) 
  !**********************************************************
  ! PRINT_VECTOR_DAT prints vector in a dat file  
  ! Usage: 
  !	call print_matrix_dat(namefile,v,n)
    
  ! INPUTS
  !	namefile     : name of the dat file 
  !	v            : vector to print  
  !     n            : size of the vector
  !
  !***********************************************************
    !Local
    integer :: i
    
    !Dummy
    real, intent(in) :: v(:)
    character (len=*) :: namefile

    open(1, file=namefile, action='write', status = 'replace')

    do i = 1,size(v)
       write(1,*) v(i)
    end do
    close(1)
    
  end subroutine print_vector_dat

  !**********************************************************
  ! LOG TO REAL 
  !********************************************************** 
  function logtoreal(v) result(vout)
  !**********************************************************
  ! LOGTOREAL converts a logical vector into a real vector  
  ! Usage: 
  !	vout = logtoreal(v,n)
    
  ! INPUTS
  !	v     : # logical vector 
  !	n     : # size of the vector
    
  ! OUTPUTS	
  !	vout  : real vector 
  !***********************************************************
    !Local
    integer :: i 
    !Dummy
    integer :: n
    logical :: v(:)

    real, allocatable :: vout(:)
    allocate(vout(size(v)))
    
    vout = 0.
    do i = 1,size(v)
       if(v(i))then
          vout(i) = 1.
       end if
    end do
  end function logtoreal
  
  !**********************************************************
  ! FIND 
  !********************************************************** 
  function find(mask) result(vout)
  !**********************************************************
  ! FIND gets the position in the vector according with the mask  
  ! Usage: 
  !	vout = find(mask,n)
    
  ! INPUTS
  !	mask  : # the logical vector with a mask  
  !	n     : # size of the vector 
    
  ! OUTPUTS	
  !	vout  : integer vector with the positions 
  !***********************************************************
      
    !Local
    integer :: i 
    !Dummy
    logical :: mask(:)
    
    integer, allocatable :: vout(:)
    allocate(vout(size(mask)))
    vout = 0
    do i =1,size(mask)
      if(mask(i))then 
		vout(i) = i
      end if
    end	do
  end function find
  
  !**********************************************************
  ! PRINT SCALAR IN DATA FILE  
  !********************************************************** 
  subroutine print_scalar_dat(namefile,v) 
  !**********************************************************
  ! PRINT_SCALAR_DAT prints a scalar into a data file
  ! Usage: 
  !	call print_scalar_dat(namefile,v) 
    
  ! INPUTS
  !	namefile	: name of the file with the extension format included
  !	v			: scalar
    
  ! OUTPUTS	
  !
  !***********************************************************
    !Local
    
    !Dummy
    real, intent(in) :: v
    character (len=*) :: namefile
    
    open(1, file=namefile, action='write',status='replace')
    write(1,*) v
    close(1)
  end subroutine print_scalar_dat
  
  !**********************************************************
  ! STR
  !********************************************************** 
  function str(k)
  !**********************************************************
  ! STR converts a integer to string
  ! Usage: 
  !	number = str(num)
    
  ! INPUTS
  !	num		: number to convert 
    
  ! OUTPUTS	
  ! number 	: number in string format
  !***********************************************************
  	!Local
  	!Dummy
    integer, intent(in) :: k
    character(20) :: str
    
    write (str, *) k
    str = adjustl(str)
  end function
end module basic
