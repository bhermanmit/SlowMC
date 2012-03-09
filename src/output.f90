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
  public :: print_heading

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
