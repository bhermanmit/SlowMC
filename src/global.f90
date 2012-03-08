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

  ! list history input information
  integer :: nhistories
  integer :: seed
  integer :: source_type

  ! list global vars that are set during run
  integer :: eidx   ! energy index for cross sections
  integer :: ntals  ! number of tallies in problem

contains

end module global 
