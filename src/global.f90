module global

  use materials, only: material_type
  use particle,  only: particle_type
  use tally,     only: tally_type

  implicit none
  save

  ! list all types
  type(material_type)           :: mat
  type(particle_type)           :: neut
  type(tally_type), allocatable :: tal(:)

  ! list history information
  integer :: nhistories
  integer :: seed
  integer :: source_type

contains

end subroutine
