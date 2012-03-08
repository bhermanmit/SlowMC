module physics

  implicit none
  private
  public :: perform_physics

contains

!===============================================================================
! PERFORM_PHYSICS
! Doxygen comment
!===============================================================================

  subroutine perform_physics()

    use global, only: neut,eidx

    integer :: isoidx   ! isotope index
    integer :: reactid  ! reaction id

    ! get energy index
    eidx = get_eidx()

    ! sample isotope
    isoidx = sample_isotope()

    ! sample reaction in isotope
    reactid = sample_reaction(isoidx)

    ! perform reaction
    if (reactid == 1) then ! absorption
      neut%alive = .FALSE.
    else if (reactid == 2) then ! scattering
      call elastic_scattering(isoidx)
    else
      write(*,*) "Something is wrong after isotope sampling"
      stop
    end if

  end subroutine perform_physics

!===============================================================================
! GET_EIDX
! Doxygen comment
!===============================================================================

  function get_eidx() result(eidx)

    ! formal variables
    integer :: eidx

  end function get_eidx

!===============================================================================
! SAMPLE_ISOTOPE
! Doxygen comment
!===============================================================================

  function sample_isotope() result(isoidx)

    use global, only: mat,eidx

    ! formal variables
    integer :: isoidx  ! the index of the isotope sampled

    ! local variables
    real(8), allocatable :: pmf(:)  ! probability mass function
    real(8), allocatable :: cdf(:)  ! cumulative distribution function
    real(8)              :: rn      ! sampled random number
    integer              :: i       ! iteration counter

    ! allocate pmf and cdf
    if(.not. allocated(pmf)) allocate(pmf(mat%nisotopes+1))
    if(.not. allocated(cdf)) allocate(cdf(mat%nisotopes+1))

    ! set both to zero
    pmf = 0.0_8
    cdf = 0.0_8

    ! create pmf at that energy index
    pmf(2:size(pmf)) = mat%totalxs(eidx,:)/sum(mat%totalxs(eidx,:))

    ! create cdf from pmf
    do i = 1,size(pmf)
      cdf(i) = sum(pmf(1:i))
    end do

    ! sample random number
    rn = rand(0)

    ! do linear table search on cdf to find which isotope
    do i = 1,size(cdf)
      if (rn <= cdf(i)) then
        isoidx = i - 1
        exit
      end if
    end do

    ! deallocate pmf and cdf
    if (allocated(pmf)) deallocate(pmf)
    if (allocated(cdf)) deallocate(cdf)

  end function sample_isotope

!===============================================================================
! SAMPLE_REACTION
! Doxygen comment
!===============================================================================

  function sample_reaction(isoidx) result(reactid)

    use global, only: mat,eidx

    ! formal variables
    integer :: isoidx   ! the sampled isotope index
    integer :: reactid  ! the id of the reaction type

    ! local variables
    real(8) :: pmf(3)  ! probability mass function
    real(8) :: cdf(3)  ! cumulative distribution function
    real(8) :: rn      ! sampled random number
    integer :: i       ! iteration counter

    ! set up pmf
    pmf = (/0.0_8,mat%isotopes(isoidx)%xs_capt(eidx),                          &
   &                                       mat%isotopes(isoidx)%xs_scat(eidx)/)

    ! normalize pmf
    pmf = pmf / sum(pmf)

    ! compute cdf
    do i = 1,3
      cdf(i) = sum(pmf(1:i))
    end do

    ! sample random number
    rn = rand(0)

    ! perform linear table search
    do i = 1,3
      if (rn < cdf(i)) then
        reactid = i - 1
        exit
      end if
    end do

  end function sample_reaction

!===============================================================================
! ELASTIC_SCATTERING
! Doxygen comment
!===============================================================================

  subroutine elastic_scattering(isoidx)

    use global, only: neut,mat

    ! formal variables
    integer :: isoidx ! isotope sampled index

    ! local variables
    real(8) :: rn ! sampled random number

    ! sample random number
    rn = rand(0)

    ! perform asymptotic elastic scattering
    neut%E = neut%E - neut%E*(1-mat%isotopes(isoidx)%alpha)*rn;

  end subroutine elastic_scattering

end module physics
