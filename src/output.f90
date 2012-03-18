!==============================================================================!
! MODULE: output 
!
!> @author Bryan Herman
!>
!> @brief Contains routines for outputtting major info to user
!==============================================================================!

module output

  implicit none
  private
  public :: print_heading,write_output

contains

!===============================================================================
! PRINT_HEADING
!> @brief prints the code heading and run information
!===============================================================================

  subroutine print_heading()

    use global, only: VERSION_MAJOR,VERSION_MINOR,VERSION_RELEASE

    ! local variables
    character(len=10) :: today_date
    character(len=8)  :: today_time

    ! write header
    write(*, FMT='(/9(A/))') &
  & ' .d8888b.  888                        888b     d888  .d8888b.    ',       &
  & 'd88P  Y88b 888                        8888b   d8888 d88P  Y88b   ',       &
  & 'd88P  Y88b 888                        8888b   d8888 d88P  Y88b   ',       &
  & 'd88P  Y88b 888                        8888b   d8888 d88P  Y88b   ',       &
  & ' "Y888b.   888  .d88b.  888  888  888 888Y88888P888 888          ',       &
  & '    "Y88b. 888 d88""88b 888  888  888 888 Y888P 888 888          ',       &
  & '      "888 888 888  888 888  888  888 888  Y8P  888 888    888   ',       &
  & 'Y88b  d88P 888 Y88..88P Y88b 888 d88P 888   "   888 Y88b  d88P   ',       &
  & ' "Y8888P"  888  "Y88P"   "Y8888888P"  888       888  "Y8888P"    '

    ! Write version information
    write(*, FMT=*)                                                            &
 &    '     Developed At:  Massachusetts Institute of Technology'
    write(*, FMT='(6X,"Version:",7X,I1,".",I1,".",I1)')                        &
 &        VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE

    ! Write the date and time
    call get_today(today_date, today_time)
    write(*, FMT='(6X,"Date/Time:",5X,A,1X,A)')                                &
 &      trim(today_date), trim(today_time)


    ! write out divider
    write(*,FMT='(A/)') '------------------------------------------------------'

  end subroutine print_heading

!===============================================================================
! WRITE_OUTPUT
!> @brief routine that writes timing info and hdf5 file
!===============================================================================

  subroutine write_output()

    use global, only: time_init,time_run,tal,n_tallies,ana_kinf_mean,          &
   &                  ana_kinf_std,col_kinf_mean,col_kinf_std
    use hdf5

    ! local variables
    integer                        :: i            ! loop counter
    integer                        :: error        ! hdf5 error
    integer(HID_T)                 :: hdfile       ! hdf5 file
    integer(HID_T)                 :: dataspace_id ! dataspace identifier
    integer(HID_T)                 :: dataset_id   ! dataset identifier
    integer(HID_T)                 :: group_id     ! group id
    integer(HSIZE_T), dimension(1) :: dim1         ! vector for hdf5 dims
    integer(HSIZE_T), dimension(2) :: dim2         ! matrix for hdf5 dims
    character(11)                  :: talnum       ! tally number

    ! write results header
    write(*,'(/A,/,A,/)') "Results","-----------------------"

    ! write timing information
    write(*,100) "Initialization time",time_init%elapsed
    write(*,100) "Transport time",time_run%elapsed
    write(*,*)

    ! format for time write statements
100 format (1X,A,T35,"= ",ES11.4," seconds")

    ! write k-inf information
    write(*,101) 'k-inf (analog):',ana_kinf_mean,ana_kinf_std
    write(*,101) 'k-inf (coll):  ',col_kinf_mean,col_kinf_std
    write(*,*)

    ! format for kinf write statements
101 format(1X,A,2X,F7.5,1X,'+/-',1X,F7.5)

    ! open up output hdf5 file
    call h5fcreate_f("output.h5",H5F_ACC_TRUNC_F,hdfile,error)

    ! begin loop around tallies to write out
    do i = 1,n_tallies

      ! get tally number
      write (talnum, '(I11)') i
      talnum = adjustl(talnum) 

      ! open up a group
      call h5gcreate_f(hdfile,"tally_"//trim(talnum),group_id,error)

      ! write mean
      dim2 = (/size(tal(i)%mean,1),size(tal(i)%mean,2)/)
      call h5screate_simple_f(2,dim2,dataspace_id,error)
      call h5dcreate_f(hdfile,"tally_"//trim(talnum)//"/mean",H5T_NATIVE_DOUBLE&
     &                ,dataspace_id,dataset_id,error)
      call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,tal(i)%mean,dim2,error)
      call h5sclose_f(dataspace_id,error)
      call h5dclose_f(dataset_id,error)

      ! write standard deviation 
      dim2 = (/size(tal(i)%std,1),size(tal(i)%std,2)/)
      call h5screate_simple_f(2,dim2,dataspace_id,error)
      call h5dcreate_f(hdfile,"tally_"//trim(talnum)//"/std",H5T_NATIVE_DOUBLE &
     &                ,dataspace_id,dataset_id,error)
      call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,tal(i)%std,dim2,error)
      call h5sclose_f(dataspace_id,error)
      call h5dclose_f(dataset_id,error)

      ! only write energy edges if a user tally
      if (.not.tal(i)%flux_tally) then

        ! write tally data to file
        dim1 = (/size(tal(i)%E)/)
        call h5screate_simple_f(1,dim1,dataspace_id,error)
        call h5dcreate_f(hdfile,"tally_"//trim(talnum)//"/E",H5T_NATIVE_DOUBLE,&
       &                 dataspace_id,dataset_id,error)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,tal(i)%E,dim1,error)
        call h5sclose_f(dataspace_id,error)
        call h5dclose_f(dataset_id,error)

      end if

      ! close the group
      call h5gclose_f(group_id,error)

    end do

    ! close the file
    call h5fclose_f(hdfile,error)

  end subroutine write_output

!===============================================================================
! GET_TODAY
!> @brief calculates information about date/time of run
!===============================================================================

  subroutine get_today(today_date, today_time)

    character(10), intent(out) :: today_date
    character(8),  intent(out) :: today_time

    integer       :: val(8)
    character(8)  :: date_
    character(10) :: time_
    character(5)  :: zone

    call date_and_time(date_, time_, zone, val)
    ! val(1) = year (YYYY)
    ! val(2) = month (MM)
    ! val(3) = day (DD)
    ! val(4) = timezone
    ! val(5) = hours (HH)
    ! val(6) = minutes (MM)
    ! val(7) = seconds (SS)
    ! val(8) = milliseconds

    if (val(2) < 10) then
       if (val(3) < 10) then
          today_date = date_(6:6) // "/" // date_(8:8) // "/" // date_(1:4)
       else
          today_date = date_(6:6) // "/" // date_(7:8) // "/" // date_(1:4)
       end if
    else
       if (val(3) < 10) then
          today_date = date_(5:6) // "/" // date_(8:8) // "/" // date_(1:4)
       else
          today_date = date_(5:6) // "/" // date_(7:8) // "/" // date_(1:4)
       end if
    end if
    today_time = time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)

  end subroutine get_today

end module output
