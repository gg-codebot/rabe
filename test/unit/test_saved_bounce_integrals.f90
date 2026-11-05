program test_saved_bounce_integrals
    ! Quick checks of reading back saved bounce integrals; resuming a real run
    ! is covered by the end-to-end test in test/restart.
    use constants, only: dp
    use grid_mod, only: integration_grid_t, fieldline_for_precession_t
    use grid_mod, only: set_integration_grids
    use saved_bounce_integrals, only: bounce_integrals_header_t
    use saved_bounce_integrals, only: write_bounce_integrals_header
    use saved_bounce_integrals, only: append_bounce_integrals, read_bounce_integrals

    implicit none

    character(len=*), parameter :: test_file = "test_saved_bounce_integrals.dat"
    integer, parameter :: n_fieldlines = 3
    integer, parameter :: n_grid = 7
    real(dp), parameter :: eta_c = 0.1_dp
    real(dp), parameter :: eta_t = 0.3_dp
    real(dp), parameter :: s_tor = 0.25_dp
    real(dp), parameter :: M_pol = -1.0_dp
    real(dp), parameter :: N_tor = 4.0_dp
    real(dp), parameter :: test_value = 1.0_dp/3.0_dp ! not exact in decimal
    character(len=*), parameter :: field_file = "wout_test.nc"
    character(len=*), parameter :: git_hash = "abc1234"
    character(len=*), parameter :: other_git_hash = "def5678"
    integer, parameter :: line_length = 163 ! 162 characters + newline

    type(bounce_integrals_header_t) :: header
    type(integration_grid_t) :: grid
    type(fieldline_for_precession_t), dimension(n_fieldlines) :: fieldlines
    logical, dimension(n_fieldlines) :: is_complete
    logical :: test_failed

    test_failed = .false.

    call set_integration_grids(eta_t, eta_c, grid, n_grid)
    grid%bounce_time_weighted = test_value
    grid%I_j = test_value
    grid%radial_drift_weighted = test_value
    grid%poloidal_drift_weighted = test_value
    fieldlines%grid = grid

    header = bounce_integrals_header_t(git_hash=git_hash, &
                                       field_file=field_file, &
                                       s_tor=s_tor, &
                                       M_pol=M_pol, &
                                       N_tor=N_tor, &
                                       n_grid=n_grid, &
                                       max_n_fieldlines=n_fieldlines, &
                                       n_fieldlines=n_fieldlines, &
                                       eta_c=eta_c, &
                                       eta_t=eta_t)

    ! cut inside the last line: the last fieldline must be dropped
    call write_all(test_file, header, fieldlines)
    call remove_trailing_bytes(test_file, line_length/2)
    call read_bounce_integrals(test_file, header, fieldlines, is_complete)
    call check_complete("cut inside line", is_complete, [.true., .true., .false.])

    ! cut between lines: every line is complete, but the last fieldline
    ! lacks its last line and must be dropped
    call write_all(test_file, header, fieldlines)
    call remove_trailing_bytes(test_file, line_length)
    fieldlines%grid = grid
    fieldlines(1)%grid%I_j = 0.0_dp
    fieldlines(2)%grid%I_j = 0.0_dp
    call read_bounce_integrals(test_file, header, fieldlines, is_complete)
    call check_complete("cut between lines", is_complete, [.true., .true., .false.])

    ! values of complete fieldlines must come back bit-exact
    if (any(fieldlines(1)%grid%I_j /= test_value) .or. &
        any(fieldlines(2)%grid%I_j /= test_value)) then
        print *, "test_saved_bounce_integrals failed: values not bit-exact"
        test_failed = .true.
    end if

    ! a different git hash only warns
    header%git_hash = other_git_hash
    call read_bounce_integrals(test_file, header, fieldlines, is_complete)
    call check_complete("other git hash", is_complete, [.true., .true., .false.])

    if (test_failed) error stop

contains

    subroutine write_all(filename, header, fieldlines)
        character(len=*), intent(in) :: filename
        type(bounce_integrals_header_t), intent(in) :: header
        type(fieldline_for_precession_t), dimension(:), intent(in) :: fieldlines
        integer :: idx

        call write_bounce_integrals_header(filename, header)
        do idx = 1, size(fieldlines)
            call append_bounce_integrals(filename, idx, fieldlines(idx)%grid)
        end do
    end subroutine write_all

    subroutine check_complete(label, is_complete, expected_complete)
        character(len=*), intent(in) :: label
        logical, dimension(:), intent(in) :: is_complete, expected_complete

        if (any(is_complete .neqv. expected_complete)) then
            print *, "-------------------------------------------------------------"
            print *, "test_saved_bounce_integrals failed: ", label
            print *, "complete fieldlines found:    ", is_complete
            print *, "complete fieldlines expected: ", expected_complete
            test_failed = .true.
        end if
    end subroutine check_complete

    subroutine remove_trailing_bytes(filename, n_bytes)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n_bytes

        character(len=:), allocatable :: content
        integer :: file_unit, file_size

        open (newunit=file_unit, file=filename, access='stream', &
              form='unformatted', status='old', action='read')
        inquire (unit=file_unit, size=file_size)
        allocate (character(len=file_size) :: content)
        read (file_unit) content
        close (file_unit)

        open (newunit=file_unit, file=filename, access='stream', &
              form='unformatted', status='replace', action='write')
        write (file_unit) content(1:file_size - n_bytes)
        close (file_unit)
    end subroutine remove_trailing_bytes

end program test_saved_bounce_integrals
