module materials

  implicit none
  private

  type, public ::  material_type

    type(iso_type), allocatable :: isotopes(:)   ! 1-D array of isotopes in mat
    integer                     :: nisotopes     ! number of isotopes in mat
    real(8), allocatable        :: totalxs(:,:)  ! array of macroscopic tot xs

  end type material_type

  type :: iso_type

    real(8)               :: N          ! number density
    real(8)               :: A          ! atomic weight
    real(8)               :: alpha      ! (A-1)^2/(A+1)^2
    real(8), allocatable  :: xs_capt(:) ! capture micro xs
    real(8), allocatable  :: xs_scat(:) ! scattering micro xs

    logical               :: thermal = .false. ! thermal scatterer
    type(thermal_type)    :: thermal_lib

  end type iso_type

  type :: thermal_type

    real(8) :: kTvec  ! vector of kT values
    real(8) :: cdf    ! matrix of cdf values
! put other stuff like spacing and things in here

  end type thermal_type

contains

!==============================================================================
! SET_UP_MATERIALS
! Doxygen comment
!==============================================================================

  subroutine set_up_material(this,nisotopes)

    ! allocate isotopes array
    if (.not. allocated(this%isotopes)) allocate(this%isotopes(nisotopes))

    ! set size
    this%nisotopes = nisotopes

  end subroutine set_up_material

!==============================================================================
! LOAD_ISOTOPE
! Doxygen comment
!==============================================================================

  subroutine load_isotope(this,path,N,A)

    ! open up hdf5 file

    ! read size of vector

    ! allocate all xs vectors

    ! zero out xs vectors

    ! read in xs

    ! set other parameters

  end subroutine load_isotope

end module materials
