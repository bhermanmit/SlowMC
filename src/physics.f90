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
    idx = ceiling(rn / mat(1)%source%cdf_width) + 1

    ! bounds checker
    if (idx > size(mat(1)%source%E)) then
      write(*,*) 'Bounds error on source samplings'
      write(*,*) 'Random number:',rn
      write(*,*) 'Index Location:',idx
      stop
    end if

    ! extract that E and set it to neutron
    neut%E = mat(1)%source%E(idx)

  end subroutine sample_source

!===============================================================================
! PERFORM_PHYSICS
!> @brief high level routine to perform transport physics
!===============================================================================

  subroutine perform_physics()

    use global, only: neut,add_to_tallies,n_abs,n_fiss

    ! sample region
    neut%region = sample_region()

    ! a collision has now occurred in a region at an energy, add to tally
    call add_to_tallies()

    ! sample isotope
    neut%isoidx = sample_isotope(neut%region)

    ! sample reaction in isotope
    neut%reactid = sample_reaction(neut%region,neut%isoidx)

    ! perform reaction
    if (neut%reactid == 1 .or. neut%reactid == 2) then ! absorption
      neut%alive = .FALSE.
      if (neut%reactid == 1) n_fiss = n_fiss+1
      n_abs = n_abs + 1
    else if (neut%reactid == 3) then ! scattering
      call elastic_scattering(neut%region,neut%isoidx)
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
    eidx = ceiling((log10(E) - log10(mat(1)%E_min))/mat(1)%E_width) + 1

    ! check bounds
    if (eidx == 0 .or. eidx >=mat(1)%npts) then
      write(*,*) 'Energy index out of bounds!'
      write(*,*) 'Energy:',E
      write(*,*) 'Width:',mat(1)%E_width
      stop
   end if

  end function get_eidx

!===============================================================================
! SAMPLE_REGION
!> @brief function to sample region where interaction occurs
!===============================================================================

  function sample_region() result(region)

    use global, only: n_materials,Dancoff,res_iso,radius,mat,eidx,neut

    ! formal variables
    integer :: region ! region of interaction

    ! local variables
    real(8) :: Pff   ! fuel-to-fuel collision probability
    real(8) :: Pfm   ! fuel-to-moderator collision probability
    real(8) :: Pmf   ! moderator-to-fuel collision probability
    real(8) :: A     ! A factor
    real(8) :: a1    ! alpha 1 factor
    real(8) :: a2    ! alpha 2 factor
    real(8) :: b     ! beta factor
    real(8) :: sig_e ! macro escape cross section
    real(8) :: sig_t ! macro total cross section of resonant isotope
    real(8) :: rn    ! sampled random number

    ! set region number
    region = 1

    ! check for more than 1 material
    if (n_materials == 2) then

      ! calculate A
      A = (1.0_8-Dancoff)/Dancoff

      ! calculate alpha 1
      a1 = ((5._8*A+6._8)-sqrt(A**2+36._8*A+36._8))/(2._8*(A+1._8))

      ! calculate alpha 2
      a2 = ((5._8*A+6._8)+sqrt(A**2+36._8*A+36._8))/(2._8*(A+1._8))

      ! calculate beta
      b = (((4._8*A+6._8)/(A+1._8)) - a1)/(a2 - a1)

      ! calculate macro escape cross section
      sig_e = 1._8/(2._8*radius)

      ! get macro total cross section of resonant isotope
      sig_t = sum(mat(1)%totalxs(eidx,:))

      ! compute fuel-to-fuel collision probability (Carlviks two-term)
      Pff = (b*sig_t)/(a1*sig_e + sig_t) + ((1-b)*sig_t)/(a2*sig_e + sig_t)

      ! compute fuel-to-moderator collision probability
      Pfm = 1._8 - Pff

      ! using reciprocity compute moderator-to-fuel collision probability
      Pmf = Pfm * (sum(mat(1)%totalxs(eidx,:))*mat(1)%vol) /                   &
     &            (sum(mat(2)%totalxs(eidx,:))*mat(2)%vol)

      ! sample random number
      rn = rand(0)

      ! figure out what region currently in and sample accordingly
      if (neut%region == 1) then
        if (rn < Pfm) then
          region = 2
        else
          region = 1
        end if
      else if (neut%region == 2) then
        if (rn < Pmf) then
          region = 1
        else
          region = 2
        end if
      else
        write(*,*) 'Cant find neutron!'
        stop
      end if

    end if
!print *,'Pff:',Pff,'Pfm:',Pfm,'Pmf:',Pmf,'Energy:',neut%E,sig_t,sum(mat(2)%totalxs(eidx,:)) 
!read*
  end function sample_region

!===============================================================================
! SAMPLE_ISOTOPE
!> @brief function to sample interaction isotope
!===============================================================================

  function sample_isotope(region) result(isoidx)

    use global, only: mat,eidx

    ! formal variables
    integer :: region  ! region of interaction
    integer :: isoidx  ! the index of the isotope sampled

    ! local variables
    real(8), allocatable :: pmf(:)  ! probability mass function
    real(8), allocatable :: cdf(:)  ! cumulative distribution function
    real(8)              :: rn      ! sampled random number
    integer              :: i       ! iteration counter

    ! allocate pmf and cdf
    if(.not. allocated(pmf)) allocate(pmf(mat(region)%nisotopes+1))
    if(.not. allocated(cdf)) allocate(cdf(mat(region)%nisotopes+1))

    ! set both to zero
    pmf = 0.0_8
    cdf = 0.0_8

    ! create pmf at that energy index
    pmf(2:size(pmf)) = mat(region)%totalxs(eidx,:) /                           &
   &                   sum(mat(region)%totalxs(eidx,:))

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

  function sample_reaction(region,isoidx) result(reactid)

    use global, only: mat,eidx

    ! formal variables
    integer :: region   ! region of interaction
    integer :: isoidx   ! the sampled isotope index
    integer :: reactid  ! the id of the reaction type

    ! local variables
    real(8) :: pmf(4)  ! probability mass function
    real(8) :: cdf(4)  ! cumulative distribution function
    real(8) :: rn      ! sampled random number
    integer :: i       ! iteration counter

    ! set up pmf
    pmf = (/0.0_8,mat(region)%isotopes(isoidx)%xs_fiss(eidx),                  &
   &              mat(region)%isotopes(isoidx)%xs_capt(eidx),                  &
                  mat(region)%isotopes(isoidx)%xs_scat(eidx)/)

    ! normalize pmf
    pmf = pmf / sum(pmf)

    ! compute cdf
    do i = 1,4
      cdf(i) = sum(pmf(1:i))
    end do

    ! sample random number
    rn = rand(0)

    ! perform linear table search
    do i = 1,4
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

  subroutine elastic_scattering(region,isoidx)

    use global, only: neut,mat,kT

    ! formal variables
    integer :: region ! region of interaction
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
    if (neut%E < 4e-6_8 .and. mat(region)%isotopes(isoidx)%thermal) then

      ! get index in cdf
      idx = ceiling(rn/mat(region)%isotopes(isoidx)%thermal_lib%cdf_width)

      ! check index
      if (idx == 0) idx = 1

      ! preallocate energy vector
      if (.not.allocated(Evec))                                                &
     & allocate(Evec(size(mat(region)%isotopes(isoidx)%thermal_lib%kTvec)))

      ! set possible energy ratios vector
      Evec = mat(region)%isotopes(isoidx)%thermal_lib%Erat(idx,:)

      ! get energy in kT units
      EkT = neut%E/kT

      ! find index in kT space
      do i = 1,size(mat(region)%isotopes(isoidx)%thermal_lib%kTvec)
        if (EkT < mat(region)%isotopes(isoidx)%thermal_lib%kTvec(i)) then
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
       &mat(region)%isotopes(isoidx)%thermal_lib%kTvec(kTidx-1))*((Evec(kTidx) &
       &- Evec(kTidx-1))/(mat(region)%isotopes(isoidx)%thermal_lib%kTvec(kTidx)&
       &- mat(region)%isotopes(isoidx)%thermal_lib%kTvec(kTidx-1)))

       ! multiply by incoming energy
       neut%E = neut%E*Eint

      end if

      ! deallocate energy vector
      if(allocated(Evec)) deallocate(Evec)

    else

      ! perform asymptotic elastic scattering
      neut%E = neut%E - neut%E*(1-mat(region)%isotopes(isoidx)%alpha)*rn;

    end if

  end subroutine elastic_scattering

end module physics
