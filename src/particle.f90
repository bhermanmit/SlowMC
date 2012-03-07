module particle

  implicit none
  private

  type, public :: particle_type

    real(8) :: E      ! particle's energy
    logical :: alive  ! am i alive?

  end type particle_type

contains

end module particle
