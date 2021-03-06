!==============================================================================!
! MODULE: global 
!
!> @author Bryan Herman
!>
!> @brief Contains all of the global variables
!==============================================================================!

module global

  use materials, only: material_type
  use particle,  only: particle_type
  use tally,     only: tally_type
  use timing,    only: Timer

  implicit none
  save

  ! version information
  integer :: VERSION_MAJOR   = 0
  integer :: VERSION_MINOR   = 1
  integer :: VERSION_RELEASE = 1

  ! list all types
  type(particle_type)              :: neut
  type(material_type), allocatable :: mat(:)
  type(tally_type), allocatable    :: tal(:)

  ! list history input information
  integer :: nhistories
  integer :: seed
  integer :: source_type

  ! list global vars that are set during run
  integer :: eidx        ! energy index for cross sections
  integer :: n_tallies   ! number of tallies
  integer :: n_materials ! n_materials
  integer :: res_iso     ! resonant isotope id in material 1
  real(8) :: Dancoff     ! lattice Dancoff factor (C)
  real(8) :: radius      ! radius of fuel pin

  ! set max and min energy
  real(8) :: emin = 1e-11_8
  real(8) :: emax = 20.0_8

  ! kT value base on 300K
  real(8) :: kT = 8.6173324e-5_8*300*1.0e-6_8

  ! set nu value
  real(8) :: nubar = 2.455_8

  ! timers
  type(Timer) :: time_init
  type(Timer) :: time_run

  ! analog counters for k-inf
  integer :: n_abs=0.0_8
  integer :: n_fiss=0.0_8
  real(8) :: ana_kinf_mean = 0.0_8
  real(8) :: ana_kinf_std  = 0.0_8
  real(8) :: col_kinf_mean = 0.0_8
  real(8) :: col_kinf_std  = 0.0_8

contains

!==============================================================================
! ALLOCATE_PROBLEM
!> @brief allocates global variables for calculation
!==============================================================================

  subroutine allocate_problem()

    ! formal variables

    ! allocate tallies
    if (.not.allocated(tal)) allocate(tal(n_tallies))
    if (.not.allocated(mat)) allocate(mat(n_materials))

  end subroutine allocate_problem

!===============================================================================
! DEALLOCATE_PROBLEM
!> @brief deallocates global variables
!===============================================================================

  subroutine deallocate_problem()

    use materials, only: deallocate_material
    use tally,     only: deallocate_tally

    ! local variables
    integer :: i  ! loop counter

    ! deallocate within materials
    do i = 1,n_materials

      ! deallocate material
      call deallocate_material(mat(i))

    end do

    ! deallocate material variable
    if (allocated(mat)) deallocate(mat)

    ! deallocate within tallies
    do i = 1,n_tallies

      ! deallocate tally
      call deallocate_tally(tal(i))

    end do

    ! deallocate tally variable
    if (allocated(tal)) deallocate(tal)

  end subroutine deallocate_problem

!===============================================================================
! COMPUTE_MACRO_CROSS_SECTIONS
!> @brief routine that handles the call to compute macro cross sections
!===============================================================================

  subroutine compute_macro_cross_sections()

    use materials, only: compute_macroxs

    ! local variables
    integer :: i  ! loop counter

    ! begin loop over materals
    do i = 1,n_materials

      ! call routine to compute xs
      call compute_macroxs(mat(i))

    end do

  end subroutine compute_macro_cross_sections

!===============================================================================
! ADD_TO_TALLIES
!> @brief routine that adds temporary value to tallies
!===============================================================================

  subroutine add_to_tallies()

    use tally, only: add_to_tally

    ! local variables
    integer :: i            ! loop counter
    real(8) :: fact = 1.0_8 ! multiplier factor
    real(8) :: totxs        ! total macroscopic xs of material
    real(8) :: mubar        ! average cosine scattering angle

    ! compute macroscopic cross section
    totxs = sum(mat(neut%region)%totalxs(eidx,:))

    ! begin loop over tallies
    do i = 1,n_tallies

      ! set multiplier
      select case(tal(i)%react_type)

        ! flux only
        case(0)
          fact = 1.0_8

        ! total 
        case(1)
          fact = totxs

        ! absorption
        case(2)
          fact = sum(mat(neut%region)%absorxs(eidx,:))

        ! scattering
        case(3)
          fact = sum(mat(neut%region)%scattxs(eidx,:))

        ! nufission
        case(4)
          fact = nubar*sum(mat(neut%region)%fissixs(eidx,:))

        ! fission
        case(5)
          fact = sum(mat(neut%region)%fissixs(eidx,:))

        ! diffusion coefficient
        case(6)
          fact = 1._8/(3._8*sum(mat(neut%region)%transxs(eidx,:)))

        ! transport 
        case(7)
          fact = sum(mat(neut%region)%transxs(eidx,:))

        ! micro capture
        case(8)
          if (neut%region == tal(i)%region) then
            fact = mat(neut%region)%isotopes(tal(i)%isotope)%xs_capt(eidx)
          else
            cycle
          end if 

        ! default is flux tally
        case DEFAULT
          fact = 1.0_8

      end select

      ! call routine to add tally
      call add_to_tally(tal(i),fact,totxs,neut%E,neut%region)

    end do

  end subroutine add_to_tallies

!===============================================================================
! BANK_TALLIES
!> @brief routine that record temporary history information in tallies 
!===============================================================================

  subroutine bank_tallies()

    use tally, only: bank_tally

    ! local variables
    integer :: i  ! loop counter

    ! begin loop over tallies
    do i = 1,n_tallies

      ! call routine to bank tally
      call bank_tally(tal(i))

    end do

  end subroutine bank_tallies

!===============================================================================
! FINALIZE_TALLIES
!> @brief routine that calls another routine to compute tally statistics
!===============================================================================

  subroutine finalize_tallies()

    use tally, only: calculate_statistics

    ! local variables
    integer :: i ! loop counter
    integer :: j ! loop counter

    ! begin loop over tallies
    do i = 1,n_tallies

      ! call routine to compute statistics
      call calculate_statistics(tal(i),nhistories)

      ! normalize by volumes and histories if flux tally
      if (tal(i)%flux_tally .or. tal(i)%dv) then
        do j = 1,n_materials
          tal(i)%mean(:,j) = tal(i)%mean(:,j) / mat(j)%vol
        end do
      end if

    end do

    ! compute k_inf
    ana_kinf_mean = dble(n_fiss)*nubar/dble(nhistories)
    ana_kinf_std  = nubar*sqrt((dble(n_fiss)/dble(nhistories) -                &
   &                (dble(n_fiss)/dble(nhistories))**2)/dble(nhistories-1))
    col_kinf_mean = sum(tal(n_tallies)%mean)
    col_kinf_std  = sum(tal(n_tallies)%std)

  end subroutine finalize_tallies

end module global
