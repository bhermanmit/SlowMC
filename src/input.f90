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
   &                            allocate_problem,tal,n_tallies,n_materials
    use materials,        only: setup_material,load_source,load_isotope
    use tally,            only: set_user_tally,set_spectrum_tally
    use xml_data_input_t

    ! local variables
    logical                        :: file_exists  ! see if file exists
    character(255)                 :: filename     ! filename to open
    real(8)                        :: N            ! temp number dens
    real(8)                        :: A            ! temp atomic weight
    character(255)                 :: path         ! path to isotope file
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

      nisotopes = size(materials_%material(i)%nuclides)

      ! set up the material object
      call setup_material(mat(i),emin,emax,nisotopes)

      ! begin loop over isotope materials
      do j = 1,mat(i)%nisotopes

        ! extract info
        N = materials_%material(i)%nuclides(j)%N
        A = materials_%material(i)%nuclides(j)%A
        path = materials_%material(i)%nuclides(j)%path
        thermal = materials_%material(i)%nuclides(j)%thermal

        ! load the isotope into memory
        call load_isotope(mat(i),N,A,path,thermal)

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
        case('micro_capture')
          react_type = 3
          isotope = tallies_%tally(i)%isotope
        case DEFAULT
          react_type = 0
      end select

      ! preallocate Ebins
      if(.not. allocated(Ebins)) allocate(Ebins(size(tallies_%tally(i)%Ebins)))

      ! set Ebins
      Ebins = tallies_%tally(i)%Ebins

      ! set up user-defined tallies
      call set_user_tally(tal(i),Ebins,size(Ebins),react_type,isotope)

      ! deallocate Ebins
      if(allocated(Ebins)) deallocate(Ebins)

    end do

    ! set up spectrum tally
    call set_spectrum_tally(tal(n_tallies),emax,emin)

    ! load the source
    call load_source(mat(1),source_type,settings_%source_path)

  end subroutine read_input

end module input
