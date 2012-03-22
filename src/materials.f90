!==============================================================================!
! MODULE: materials
!
!> @author Bryan Herman
!>
!> @brief Contains information about the isotopics of problem
!==============================================================================!

module materials

  implicit none
  private
  public :: setup_material,load_source,load_isotope,compute_macroxs,           &
 &          deallocate_material

  type :: source_type

    real(8), allocatable :: E(:)      ! energy range for fission source
    real(8)              :: cdf_width ! width of cdf bins from 0 to 1

  end type source_type

  type :: thermal_type

    integer               :: kTsize    ! size of kT vector
    integer               :: cdfsize   ! size of cdf
    real(8), allocatable  :: kTvec(:)  ! vector of kT values
    real(8), allocatable  :: Erat(:,:) ! energy
    real(8)               :: cdf_width ! width of cdf interval from 0 to 1 
    
  end type thermal_type

  type :: iso_type

    real(8)               :: N           ! number density
    real(8)               :: A           ! atomic weight
    real(8)               :: alpha       ! (A-1)^2/(A+1)^2
    real(8)               :: mubar       ! average cosine scattering angle
    real(8), allocatable  :: xs_capt(:)  ! capture micro xs
    real(8), allocatable  :: xs_scat(:)  ! scattering micro xs
    real(8), allocatable  :: xs_fiss(:)  ! fission micro xs
    character(len=255)    :: name        ! name of isotope

    logical               :: thermal     ! thermal scatterer
    type(thermal_type)    :: thermal_lib ! thermal library

  end type iso_type

  type, public ::  material_type

    type(source_type)           :: source        ! the source of neutrons
    type(iso_type), allocatable :: isotopes(:)   ! 1-D array of isotopes in mat
    integer                     :: nisotopes     ! number of isotopes in mat
    integer                     :: curr_iso      ! the current isotope
    integer                     :: npts          ! number of points in energy
    real(8)                     :: E_width       ! width of energy interval
    real(8)                     :: E_min         ! min energy
    real(8)                     :: E_max         ! max energy
    real(8)                     :: vol           ! volume of region
    real(8), allocatable        :: totalxs(:,:)  ! array of macroscopic tot xs
    real(8), allocatable        :: scattxs(:,:)  ! array of macroscopic scat xs
    real(8), allocatable        :: absorxs(:,:)  ! array of macroscopic abs xs
    real(8), allocatable        :: captuxs(:,:)  ! array of macroscopic capt xs
    real(8), allocatable        :: fissixs(:,:)  ! array of macroscopic fiss xs

  end type material_type

contains

!===============================================================================
! SET_UP_MATERIALS
!> @brief routine that initializes the materials
!===============================================================================

  subroutine setup_material(this,emin,emax,nisotopes,vol)

    ! formal variables
    type(material_type) :: this      ! a material
    real(8)             :: emin      ! minimum energy to consider
    real(8)             :: emax      ! maximum energy to consider
    real(8)             :: vol       ! volume of material
    integer             :: nisotopes ! number of isotopes

    ! set number of isotopes
    this%nisotopes = nisotopes

    ! set volume
    this%vol = vol

    ! allocate isotopes array
    if (.not. allocated(this%isotopes)) allocate(this%isotopes(this%nisotopes))

    ! set up current isotope index
    this%curr_iso = 1

    ! set energy bounds
    this%E_min = emin
    this%E_max = emax

  end subroutine setup_material

!===============================================================================
! LOAD_ISOTOPE
!> @brief routine that loads isotope properties, xs, etc. into memory
!===============================================================================

  subroutine load_isotope(this,N,A,path,thermal,name)

    use hdf5

    ! formal variables
    type(material_type),target :: this    ! a material
    real(8)                    :: N       ! number density
    real(8)                    :: A       ! atomic weight
    character(len=255)         :: path    ! path to isotope
    character(len=255)         :: name    ! name of isotope
    logical                    :: thermal ! contains a thermal lib

    ! local variables
    integer                        :: error        ! hdf5 error 
    integer(HID_T)                 :: hdf5_file    ! hdf5 file id
    integer(HID_T)                 :: dataset_id   ! hdf5 dataset id
    integer(HSIZE_T), dimension(1) :: dim1         ! dimension of hdf5 var
    integer(HSIZE_T), dimension(2) :: dim2         ! dimension of hdf5 var
    integer                        :: vecsize      ! vector size
    type(thermal_type), pointer    :: therm

    ! display to user
    write(*,*) 'Loading isotope: ',trim(name)

    ! set parameters
    this%isotopes(this%curr_iso)%N = N
    this%isotopes(this%curr_iso)%A = A
    this%isotopes(this%curr_iso)%mubar = 2._8/(3._8*A) 
    this%isotopes(this%curr_iso)%alpha = ((A-1._8)/(A+1._8))**2 
    this%isotopes(this%curr_iso)%thermal = thermal
    this%isotopes(this%curr_iso)%name = name

    ! open up hdf5 file
    call h5fopen_f(trim(path),H5F_ACC_RDWR_F,hdf5_file,error)

    ! read size of vector
    call h5dopen_f(hdf5_file,"/vecsize",dataset_id,error)
    dim1 = (/1/)
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,vecsize,dim1,error)
    call h5dclose_f(dataset_id,error)

    ! allocate all xs vectors
    if (.not.allocated(this%isotopes(this%curr_iso)%xs_scat))                  &
   &         allocate(this%isotopes(this%curr_iso)%xs_scat(vecsize))
    if (.not.allocated(this%isotopes(this%curr_iso)%xs_capt))                  &
   &         allocate(this%isotopes(this%curr_iso)%xs_capt(vecsize))
    if (.not.allocated(this%isotopes(this%curr_iso)%xs_fiss))                  &
   &         allocate(this%isotopes(this%curr_iso)%xs_fiss(vecsize))


    ! keep the size
    this%npts = vecsize

    ! zero out xs vectors
    this%isotopes(this%curr_iso)%xs_scat = 0.0_8
    this%isotopes(this%curr_iso)%xs_capt = 0.0_8
    this%isotopes(this%curr_iso)%xs_fiss = 0.0_8

    ! read in xs
    call h5dopen_f(hdf5_file,"/xs_scat",dataset_id,error)
    dim1 = (/vecsize/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,                               &
   &     this%isotopes(this%curr_iso)%xs_scat,dim1,error)
    call h5dclose_f(dataset_id,error)
    call h5dopen_f(hdf5_file,"/xs_capt",dataset_id,error)
    dim1 = (/vecsize/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,                               &
   &     this%isotopes(this%curr_iso)%xs_capt,dim1,error)
    call h5dclose_f(dataset_id,error)
    call h5dopen_f(hdf5_file,"/xs_fiss",dataset_id,error)
    dim1 = (/vecsize/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,                               &
   &     this%isotopes(this%curr_iso)%xs_fiss,dim1,error)
    call h5dclose_f(dataset_id,error)

    ! get energy interval width
    call h5dopen_f(hdf5_file,"/E_width",dataset_id,error)
    dim1 = (/1/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%E_width,dim1,error)
    call h5dclose_f(dataset_id,error)

    ! check for thermal scattering kernel and load that
    if (this%isotopes(this%curr_iso)%thermal) then

      ! set pointer
      therm => this%isotopes(this%curr_iso)%thermal_lib

      ! load sizes
      call h5dopen_f(hdf5_file,"/kTsize",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,therm%kTsize,dim1,error)
      call h5dclose_f(dataset_id,error)
      call h5dopen_f(hdf5_file,"/cdfsize",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,therm%cdfsize,dim1,error)
      call h5dclose_f(dataset_id,error)

      ! read in cdf width
      call h5dopen_f(hdf5_file,"/cdf_width",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,therm%cdf_width,dim1,error)
      call h5dclose_f(dataset_id,error)

      ! preallocate vectors
      if(.not.allocated(therm%kTvec)) allocate(therm%kTvec(therm%kTsize))
      if(.not.allocated(therm%Erat)) allocate(therm%Erat(therm%cdfsize,therm%kTsize))

      ! read in vectors
      call h5dopen_f(hdf5_file,"/kT",dataset_id,error)
      dim1 = (/therm%kTsize/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,therm%kTvec,dim1,error)
      call h5dclose_f(dataset_id,error)
      call h5dopen_f(hdf5_file,"/Erat",dataset_id,error)
      dim2 = (/therm%cdfsize,therm%kTsize/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,therm%Erat,dim2,error)
      call h5dclose_f(dataset_id,error)

    end if

    ! close hdf5 file
    call h5fclose_f(hdf5_file,error)

    ! increment isotope counter
    this%curr_iso = this%curr_iso + 1

  end subroutine load_isotope

!===============================================================================
! LOAD_SOURCE
!> @brief routine to load fission source into memory
!===============================================================================

  subroutine load_source(this,source_type,source_path)

    use hdf5 

    ! formal variables
    type(material_type) :: this        ! a material
    integer             :: source_type ! 0 - fixed, 1 - fission
    character(len=255)  :: source_path ! path to source file

    ! local variables
    integer                        :: error        ! hdf5 error 
    integer(HID_T)                 :: hdf5_file    ! hdf5 file id
    integer(HID_T)                 :: dataset_id   ! hdf5 dataset id
    integer(HSIZE_T), dimension(1) :: dim1         ! dimension of hdf5 var
    integer                        :: vecsize      ! vector size for fission


    ! check for fission source
    if (source_type == 1) then

      ! open the fission source file
      call h5fopen_f(trim(source_path),H5F_ACC_RDWR_F,hdf5_file,error)

      ! open dataset and read in vector size
      call h5dopen_f(hdf5_file,"/vecsize",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,vecsize,dim1,error)
      call h5dclose_f(dataset_id,error)

      ! open dataset and read in width of cdf interval
      call h5dopen_f(hdf5_file,"/cdf_width",dataset_id,error)
      dim1 = (/1/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%source%cdf_width,dim1,  &
     &               error)
      call h5dclose_f(dataset_id,error)

      ! preallocate vectors in source object
      if(.not.allocated(this%source%E)) allocate(this%source%E(vecsize))

      ! open dataset and read in energy vector
      call h5dopen_f(hdf5_file,"/E",dataset_id,error)
      dim1 = (/vecsize/)
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%source%E,dim1,error)
      call h5dclose_f(dataset_id,error)

      ! close the file
      call h5fclose_f(hdf5_file,error)

    end if

  end subroutine load_source

!===============================================================================
! COMPUTE_MACROXS
!> @brief routine to pre-compute macroscopic cross sections 
!===============================================================================

  subroutine compute_macroxs(this)

    ! formal variables
    type(material_type),target :: this ! a material

    ! local variables
    integer                 :: i   ! loop counter
    type(iso_type), pointer :: iso ! pointer to current isotope

    ! allocate xs arrays
    if (.not.allocated(this%totalxs))                                          &
   &                            allocate(this%totalxs(this%npts,this%nisotopes))
    if (.not.allocated(this%scattxs))                                          &
   &                            allocate(this%scattxs(this%npts,this%nisotopes))
    if (.not.allocated(this%absorxs))                                          &
   &                            allocate(this%absorxs(this%npts,this%nisotopes))
    if (.not.allocated(this%captuxs))                                          &
   &                            allocate(this%captuxs(this%npts,this%nisotopes))
    if (.not.allocated(this%fissixs))                                          &
   &                            allocate(this%fissixs(this%npts,this%nisotopes))


    ! zero out total xs
    this%totalxs = 0.0_8

    ! begin loop over isotopes
    do i = 1,this%nisotopes

      ! set pointer to isotope
      iso => this%isotopes(i)

      ! multiply microscopic cross section by number density and append
      this%captuxs(:,i) = iso%N*(iso%xs_capt)
      this%fissixs(:,I) = iso%N*(iso%xs_fiss)
      this%scattxs(:,i) = iso%N*(iso%xs_scat)
      this%absorxs(:,i) = iso%N*(iso%xs_capt + iso%xs_fiss)
      this%totalxs(:,i) = iso%N*(iso%xs_capt + iso%xs_fiss + iso%xs_scat)

    end do

  end subroutine compute_macroxs

!===============================================================================
! DEALLOCATE_MATERIAL
!> @brief routine to deallocate a material
!===============================================================================

  subroutine deallocate_material(this)

    ! formal variables
    type(material_type) :: this ! a material

    ! local variables
    integer :: i ! loop counter

    ! deallocate source information
    if (allocated(this%source%E)) deallocate(this%source%E)

    ! begin loop over isotopes for deallocation
    do i = 1,this%nisotopes

      ! deallocate thermal library
      if (allocated(this%isotopes(i)%thermal_lib%kTvec)) deallocate            &
     &             (this%isotopes(i)%thermal_lib%kTvec)
      if (allocated(this%isotopes(i)%thermal_lib%Erat)) deallocate             &
     &             (this%isotopes(i)%thermal_lib%Erat)

      ! deallocate xs
      if (allocated(this%isotopes(i)%xs_scat)) deallocate                      &
     &             (this%isotopes(i)%xs_scat)
      if (allocated(this%isotopes(i)%xs_capt)) deallocate                      &
     &             (this%isotopes(i)%xs_capt)
      if (allocated(this%isotopes(i)%xs_fiss)) deallocate                      &
     &             (this%isotopes(i)%xs_fiss)
    end do

    ! deallocate isotopes
    if (allocated(this%isotopes)) deallocate(this%isotopes) 

    ! deallocate macro xs
    if (allocated(this%totalxs)) deallocate(this%totalxs)
    if (allocated(this%scattxs)) deallocate(this%scattxs)
    if (allocated(this%absorxs)) deallocate(this%absorxs)
    if (allocated(this%captuxs)) deallocate(this%captuxs)
    if (allocated(this%fissixs)) deallocate(this%fissixs)

  end subroutine deallocate_material

end module materials
