module materials

  implicit none
  private

  type, public ::  material_type

    type(source_type)           :: source        ! the source of neutrons
    type(iso_type), allocatable :: isotopes(:)   ! 1-D array of isotopes in mat
    integer                     :: nisotopes     ! number of isotopes in mat
    real(8), allocatable        :: totalxs(:,:)  ! array of macroscopic tot xs

  end type material_type

  type :: iso_type

    real(8)               :: N           ! number density
    real(8)               :: A           ! atomic weight
    real(8)               :: alpha       ! (A-1)^2/(A+1)^2
    real(8), allocatable  :: xs_capt(:)  ! capture micro xs
    real(8), allocatable  :: xs_scat(:)  ! scattering micro xs

    logical               :: thermal     ! thermal scatterer
    type(thermal_type)    :: thermal_lib ! thermal library

  end type iso_type

  type :: thermal_type

    real(8), allocatable  :: kTvec(:)  ! vector of kT values
    real(8), allocatable  :: Erat(:,:) ! energy
    real(8)               :: cdf_width ! width of cdf interval from 0 to 1 
    
  end type thermal_type

  type :: source_type

    real(8), allocatable :: E         ! energy range for fission source
    real(8)              :: cdf_width ! width of cdf bins from 0 to 1

  end type source type

contains

!===============================================================================
! SET_UP_MATERIALS
! Doxygen comment
!===============================================================================

  subroutine set_up_material(this,nisotopes)

    ! allocate isotopes array
    if (.not. allocated(this%isotopes)) allocate(this%isotopes(nisotopes))

    ! set size
    this%nisotopes = nisotopes

  end subroutine set_up_material

!===============================================================================
! LOAD_ISOTOPE
! Doxygen comment
!===============================================================================

  subroutine load_isotope(this,path,N,A)

    ! open up hdf5 file

    ! read size of vector

    ! allocate all xs vectors

    ! zero out xs vectors

    ! read in xs

    ! set other parameters

    ! check for thermal scattering kernel and load that

  end subroutine load_isotope

!===============================================================================
! LOAD_SOURCE
! Doxygen comment
!===============================================================================

  subroutine load_source(this,source_type)

    type(material_type) :: this        ! a material 
    integer             :: source_type ! 0 - fixed, 1 - fission

  end subroutine load_source

!===============================================================================
! COMPUTE_TOTXS
! Doxygen comment
!===============================================================================

  subroutine compute_totxs(this)

    ! formal variables
    type(material_type) :: this  ! a material

    ! local variables
    integer                 :: i   ! loop counter
    type(iso_type), pointer :: iso ! pointer to current isotope

    ! zero out total xs
    this%totalxs = 0.0_8

    ! begin loop over isotopes
    do i = 1:this%nisotopes

      ! set pointer to isotope
      iso => this%isotope(i)

      ! multiply microscopic cross section by number density and append
      this%totalxs = this%totalxs + iso%N*(iso%xs_capt + iso%xs_scat)

    end do

  end subroutine compute_totxs

end module materials
