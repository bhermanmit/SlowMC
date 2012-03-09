module input

  implicit none
  private
  public read_input

contains

!===============================================================================
! READ_INPUT
! Doxygen comment
!===============================================================================

  subroutine read_input

    use global,           only: nhistories,seed,source_type,mat,emin,emax
    use materials,        only: setup_material,load_source,load_isotope
    use xml_data_input_t

    ! local variables
    logical                        :: file_exists  ! see if file exists
    character(255)                 :: filename     ! filename to open
    real(8)                        :: N            ! temp number dens
    real(8)                        :: A            ! temp atomic weight
    character(255)                 :: path         ! path to isotope file
    logical                        :: thermal      ! contains thermal lib
    integer                        :: i            ! iteration counter

    ! check for input file
    filename = "input.xml"
    inquire(FILE=trim(filename), EXIST=file_exists)
    if (.not. file_exists) then
      write(*,*) 'Cannot read input file!'
      stop
    else

      ! tell user
      write(*,*) "Reading INPUT XML file..."

    end if

    ! read in input file
    call read_xml_file_input_t(trim(filename))

    ! read in settings
    nhistories = settings_%histories
    seed = settings_%seed
    source_type = settings_%source_type
    mat%nisotopes = settings_%nisotopes

    ! load the source
    call load_source(mat,source_type,settings_%source_path)

    ! set up the material object
    call setup_material(mat,emin,emax)

    ! begin loop over isotope materials
    do i = 1,size(material_)

      ! extract info
      N = material_(i)%N
      A = material_(i)%A
      path = material_(i)%path
      thermal = material_(i)%thermal

      ! load the isotope into memory
      call load_isotope(mat,N,A,path,thermal)

    end do

  end subroutine read_input

end module input
