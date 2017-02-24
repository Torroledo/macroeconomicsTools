module def_types

  implicit none
  
  type :: sparse_matrix
     real, allocatable :: values(:)
     integer, allocatable, dimension(:) :: d1,d2
  end type sparse_matrix
  
  type :: sparse_vector
     real, allocatable :: values(:)
     integer, allocatable, dimension(:) :: d1
  end type sparse_vector
  
end module def_types
