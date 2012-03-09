!==============================================================================!
! MODULE: physics
!
!> @author Bryan Herman
!>
!> @brief Contains routines to model the physics of the problem
!==============================================================================!

module physics

  implicit none
  private
  public :: sample_source,perform_physics,get_eidx

contains

!===============================================================================
! SAMPLE_SOURCE
!> @brief routine to sample source from cdf
!===============================================================================

  subroutine sample_source()

    use global, only: mat,neut

    ! local variables
    integer :: idx ! index for sampling
    real(8) :: rn  ! sampled random number

    ! sample a random number
    rn = rand(0)

    ! compute index in cdf
    idx = ceiling(rn / mat%source%cdf_width) + 1

    ! bounds checker
    if (idx > size(mat%source%E)) then
      write(*,*) 'Bounds error on source samplings'
      write(*,*) 'Random number:',rn
      write(*,*) 'Index Location:',idx
      stop
    end if

    ! extract that E and set it to neutron
    neut%E = mat%source%E(idx)

  end subroutine sample_source

!===============================================================================
! PERFORM_PHYSICS
!> @brief high level routine to perform transport physics
!===============================================================================

  subroutine perform_physics()

    use global, only: neut

    integer :: isoidx   ! isotope index
    integer :: reactid  ! reaction id

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
!> @brief function to compute the index in unionized energy grid
!===============================================================================

  function get_eidx(E) result(eidx)

    use global, only: mat

    ! formal variables
    real(8)             :: E    ! neutron's energy
    integer             :: eidx ! the energy index

    ! compute index
    eidx = ceiling((log10(E) - log10(mat%E_min))/mat%E_width) + 1

    ! check bounds
    if (eidx == 0 .or. eidx >=mat%npts) then
      write(*,*) 'Energy index out of bounds!'
      write(*,*) 'Energy:',E
      write(*,*) 'Width:',mat%E_width
      stop
   end if

  end function get_eidx

!===============================================================================
! SAMPLE_ISOTOPE
!> @brief function to sample interaction isotope
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

    ! check iso
    if (isoidx == 0) then
      isoidx = 1 
    end if

    ! deallocate pmf and cdf
    if (allocated(pmf)) deallocate(pmf)
    if (allocated(cdf)) deallocate(cdf)

  end function sample_isotope

!===============================================================================
! SAMPLE_REACTION
!> @brief function to sample reaction type
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
!> @brief routine to perform thermal/asymptotic elastic scattering physics 
!===============================================================================

  subroutine elastic_scattering(isoidx)

    use global, only: neut,mat,kT

    ! formal variables
    integer :: isoidx ! isotope sampled index

    ! local variables
    integer :: i     ! iteration counter
    integer :: idx   ! index in cdf vector
    integer :: kTidx ! index in kT vector
    real(8) :: rn    ! sampled random number
    real(8) :: EkT   ! energy / kT
    real(8) :: Eint  ! interpolated E value
    real(8), allocatable :: Evec(:)

    ! sample random number
    rn = rand(0)

    ! check for thermal scattering
    if (neut%E < 4e-6_8 .and. mat%isotopes(isoidx)%thermal) then

      ! get index in cdf
      idx = ceiling(rn/mat%isotopes(isoidx)%thermal_lib%cdf_width)

      ! check index
      if (idx == 0) idx = 1

      ! preallocate energy vector
      if (.not.allocated(Evec))                                                &
     & allocate(Evec(size(mat%isotopes(isoidx)%thermal_lib%kTvec)))

      ! set possible energy ratios vector
      Evec = mat%isotopes(isoidx)%thermal_lib%Erat(idx,:)

      ! get energy in kT units
      EkT = neut%E/kT

      ! find index in kT space
      do i = 1,size(mat%isotopes(isoidx)%thermal_lib%kTvec)
        if (EkT < mat%isotopes(isoidx)%thermal_lib%kTvec(i)) then
          kTidx = i
          exit
        end if
      end do 

      ! interpolate on energy value
      if (kTidx == 1) then
        neut%E = Evec(kTidx)
      else
        ! perform linear interplation on kT value
        Eint = Evec(kTidx-1) + (EkT -                                          &
       &       mat%isotopes(isoidx)%thermal_lib%kTvec(kTidx-1))*((Evec(kTidx)  &
       &     - Evec(kTidx-1))/(mat%isotopes(isoidx)%thermal_lib%kTvec(kTidx)   &
       &     - mat%isotopes(isoidx)%thermal_lib%kTvec(kTidx-1)))

       ! multiply by incoming energy
       neut%E = neut%E*Eint

      end if

      ! deallocate energy vector
      if(allocated(Evec)) deallocate(Evec)

    else

      ! perform asymptotic elastic scattering
      neut%E = neut%E - neut%E*(1-mat%isotopes(isoidx)%alpha)*rn;

    end if

  end subroutine elastic_scattering

end module physics
