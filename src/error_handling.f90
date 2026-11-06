module error_handling
    implicit none
    private

    public :: failed_sanity_check, set_unsafe_mode
    public :: did_fail_any_sanity_check, reset_failed_check_counter

    logical, private :: unsafe_mode = .false.
    integer, private :: failed_check_counter = 0

contains

    !> Toggle unsafe mode for sanity checks.
    !!
    !! When .true., a failed check increments the counter instead of halting,
    !! allowing surface-by-surface iteration to continue. Call
    !! did_fail_any_sanity_check to detect failures.
    subroutine set_unsafe_mode(value)
        logical, intent(in) :: value
        unsafe_mode = value
    end subroutine set_unsafe_mode

    subroutine failed_sanity_check()
        if (unsafe_mode) then
            failed_check_counter = failed_check_counter + 1
            return
        else
            error stop
        end if
    end subroutine failed_sanity_check

    !> Reset the failed sanity check counter to zero.
    subroutine reset_failed_check_counter()
        failed_check_counter = 0
    end subroutine reset_failed_check_counter

    !> Return .true. if any sanity check has failed since the last reset.
    function did_fail_any_sanity_check()
        logical :: did_fail_any_sanity_check
        did_fail_any_sanity_check = (failed_check_counter > 0)
    end function did_fail_any_sanity_check

end module error_handling
