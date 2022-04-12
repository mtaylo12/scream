module scream_abortutils

!-------------------------------------------------
!Utility to stop the model in case of
!catastrophic errors
!-------------------------------------------------
implicit none
private

!public subroutines
public :: endscreamrun

contains

  subroutine endscreamrun (msg)
    !-------------------------------------------------
    ! This subroutine will print the optional message
    ! received via optional arg "msg" and stop
    ! the simulation
    !-------------------------------------------------

    use shr_abort_mod, only: shr_abort => shr_abort_abort

    character(len=*), intent(in), optional :: msg

    if(present(msg)) then
      call shr_abort(msg,-1)
    else
      call shr_abort('ERROR: Aborting...',-1)
    endif

  end subroutine endscreamrun

end module scream_abortutils
