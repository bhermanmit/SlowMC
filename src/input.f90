!==============================================================================!
! MODULE: input
!
!> @author Bryan Herman
!>
!> @brief Handles reading in the input xml file and intializing global vars
!==============================================================================!

module input

  implicit none
  private
  public read_input

contains

!===============================================================================
! READ_INPUT
!> @brief Reads the input xml file and sets global variables
!===============================================================================

  subroutine read_input

    use global,           only: nhistories,seed,source_type,mat,emin,emax,     &
   &                            allocate_problem,tal,n_tallies,n_materials,    &
   &                            res_iso,Dancoff,radius
    use materials,        only: setup_material,load_source,load_isotope
    use tally,            only: set_user_tally,set_spectrum_tally
    use xml_data_input_t

    ! local variables
    logical                        :: file_exists  ! see if file exists
    character(len=255)             :: filename     ! filename to open
    real(8)                        :: N            ! temp number dens
    real(8)                        :: A            ! temp atomic weight
    real(8)                        :: vol          ! volume of region
    character(len=255)             :: path         ! path to isotope file
    character(len=255)             :: name         ! name of isotope
    logical                        :: thermal      ! contains thermal lib
    integer                        :: i            ! iteration counter
    integer                        :: j            ! iteration counter
    integer                        :: nisotopes    ! number of isotopes in mat
    integer                        :: react_type   ! reaction type
    integer                        :: isotope=0    ! isotope for micro mult
    real(8), allocatable           :: Ebins(:)     ! tally energy bins

    ! check for input file
    filename = "input.xml"
    inquire(FILE=trim(filename), EXIST=file_exists)
    if (.not. file_exists) then
      write(*,*) 'Cannot read input file!'
      stop
    else

      ! tell user
      write(*,'(A/)') "Reading INPUT XML file..."

    end if

    ! read in input file
    call read_xml_file_input_t(trim(filename))

    ! read in settings
    nhistories = settings_%histories
    seed = settings_%seed
    source_type = settings_%source_type

    ! get size of materials
    n_materials = size(materials_%material)

    ! get size of tallies
    if (.not.associated(tallies_%tally)) then
      n_tallies = 1
    else
      n_tallies = size(tallies_%tally) + 1
    end if

    ! allocate problem
    call allocate_problem()

    ! begin loop around materials
    do i = 1,n_materials

      ! get number of isotopes and volume
      nisotopes = size(materials_%material(i)%nuclides)
      vol = materials_%material(i)%V

      ! set homogeneous volume
      if (trim(materials_%material(i)%type)=='homogeneous') vol = 1.0_8

      ! set up the material object
      call setup_material(mat(i),emin,emax,nisotopes,vol)

      ! begin loop over isotope materials
      do j = 1,mat(i)%nisotopes

        ! check volumes and number densities
        if (trim(materials_%material(i)%type)=='homogeneous') then

          ! set volume to 1 and don't adjust n dens
          N = materials_%material(i)%nuclides(j)%N

        else if (trim(materials_%material(i)%type)=='fuel') then

          ! don't adjust n dens
          N = materials_%material(i)%nuclides(j)%N

          ! check volume
          if (abs(vol - 0.0_8) < 1e-10_8) then
            write(*,*) 'Please enter a physical fuel volume!'
            stop
          end if

        else

          ! check volume
          if (abs(vol - 0.0_8) < 1e-10_8) then
            write(*,*) 'Please enter a physical moderator volume!'
            stop
          end if

          ! adjust number density by volume weighting
          N = materials_%material(i)%nuclides(j)%N*                            &
         &   (materials_%material(i)%nuclides(j)%V/vol)

        end if

        ! extract other info
        A = materials_%material(i)%nuclides(j)%A
        path = materials_%material(i)%nuclides(j)%path
        thermal = materials_%material(i)%nuclides(j)%thermal
        name = materials_%material(i)%nuclides(j)%name 

        ! load the isotope into memory
        call load_isotope(mat(i),N,A,path,thermal,name)

        ! check for resonant isotope in material 1
        if (trim(materials_%material(i)%type)=='fuel' .and.                    &
       &    trim(settings_%res_iso) == trim(name)) then

          ! get Dancoff factor and resonant isotope
          res_iso = j
          Dancoff = settings_%Dancoff
          radius = settings_%radius

        end if

      end do

    end do

    ! begin loop over tallies
    do i = 1,n_tallies-1

      ! set reaction type
      select case(trim(tallies_%tally(i)%type))
        case('flux')
          react_type = 0
        case('absorption')
          react_type = 1
        case('scattering')
          react_type = 2
        case('nufission')
          react_type = 3
        case('micro_capture')
          react_type = 4
          isotope = tallies_%tally(i)%isotope
        case DEFAULT
          react_type = 0
      end select

      ! preallocate Ebins
      if(.not. allocated(Ebins)) allocate(Ebins(size(tallies_%tally(i)%Ebins)))

      ! set Ebins
      Ebins = tallies_%tally(i)%Ebins

      ! set up user-defined tallies
      call set_user_tally(tal(i),Ebins,size(Ebins),react_type,isotope,n_materials)

      ! deallocate Ebins
      if(allocated(Ebins)) deallocate(Ebins)

    end do

    ! set up spectrum tally
    call set_spectrum_tally(tal(n_tallies),emax,emin,n_materials)

    ! load the source
    call load_source(mat(1),source_type,settings_%source_path)

  end subroutine read_input

end module input
