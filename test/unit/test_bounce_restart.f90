program test_bounce_restart
    use constants, only: dp
    use grid_mod, only: integration_grid_t, fieldline_for_precession_t
    use grid_mod, only: set_integration_grids
    use bounce_restart, only: restart_header_t, write_restart_header
    use bounce_restart, only: append_restart_fieldline, read_restart_file

    implicit none

    character(len=*), parameter :: test_file = "test_bounce_restart.dat"
    integer, parameter :: n_fieldlines = 5
    integer, parameter :: n_grid = 7
    integer, parameter :: line_length = 163 ! 162 characters + newline

    type(restart_header_t) :: header
    type(integration_grid_t) :: grid
    type(fieldline_for_precession_t), dimension(n_fieldlines) :: expected
    type(fieldline_for_precession_t), dimension(n_fieldlines) :: fieldlines
    logical, dimension(n_fieldlines) :: is_complete
    logical :: test_failed

    test_failed = .false.

    call set_integration_grids(0.3_dp, 0.1_dp, grid, n_grid)
    expected%grid = grid
    call fill_with_test_values(expected)
    header = make_header("abc1234")

    ! round trip of all fieldlines must be bit-exact
    call write_fieldlines(test_file, header, expected, [1, 2, 3, 4, 5])
    call read_into_fresh(test_file, header, fieldlines, is_complete)
    call check_complete("round trip", is_complete, [.true., .true., .true., &
                                                     .true., .true.])
    call check_values("round trip", fieldlines, expected, is_complete)

    ! reading rewrites the file, reading it again must give the same
    call read_into_fresh(test_file, header, fieldlines, is_complete)
    call check_complete("second read", is_complete, [.true., .true., .true., &
                                                      .true., .true.])
    call check_values("second read", fieldlines, expected, is_complete)

    ! a different git hash only warns
    call read_into_fresh(test_file, make_header("def5678"), fieldlines, is_complete)
    call check_complete("other git hash", is_complete, [.true., .true., .true., &
                                                         .true., .true.])

    ! missing fieldlines are reported as incomplete
    call write_fieldlines(test_file, header, expected, [1, 2, 4])
    call read_into_fresh(test_file, header, fieldlines, is_complete)
    call check_complete("missing fieldlines", is_complete, [.true., .true., .false., &
                                                             .true., .false.])
    call check_values("missing fieldlines", fieldlines, expected, is_complete)

    ! cut inside the last line: last fieldline is dropped
    call write_fieldlines(test_file, header, expected, [1, 2, 3])
    call remove_trailing_bytes(test_file, 50)
    call read_into_fresh(test_file, header, fieldlines, is_complete)
    call check_complete("cut inside line", is_complete, [.true., .true., .false., &
                                                          .false., .false.])
    call check_values("cut inside line", fieldlines, expected, is_complete)
    call check_resume("cut inside line", test_file, header, expected)

    ! cut between lines: last fieldline lacks a line and is dropped
    call write_fieldlines(test_file, header, expected, [1, 2, 3])
    call remove_trailing_bytes(test_file, line_length)
    call read_into_fresh(test_file, header, fieldlines, is_complete)
    call check_complete("cut between lines", is_complete, [.true., .true., .false., &
                                                            .false., .false.])
    call check_resume("cut between lines", test_file, header, expected)

    ! only the final newline missing: all data is there
    call write_fieldlines(test_file, header, expected, [1, 2, 3])
    call remove_trailing_bytes(test_file, 1)
    call read_into_fresh(test_file, header, fieldlines, is_complete)
    call check_complete("missing newline", is_complete, [.true., .true., .true., &
                                                          .false., .false.])
    call check_values("missing newline", fieldlines, expected, is_complete)
    call check_resume("missing newline", test_file, header, expected)

    ! cut inside the header: start from scratch
    call write_fieldlines(test_file, header, expected, [1])
    call remove_trailing_bytes(test_file, n_grid*line_length + 20)
    call read_into_fresh(test_file, header, fieldlines, is_complete)
    call check_complete("cut inside header", is_complete, [.false., .false., .false., &
                                                            .false., .false.])
    call check_resume("cut inside header", test_file, header, expected)

    call delete_file(test_file)
    if (test_failed) error stop

contains

    function make_header(git_hash) result(header)
        character(len=*), intent(in) :: git_hash
        type(restart_header_t) :: header

        header = restart_header_t(git_hash=git_hash, &
                                  field_file="wout_test.nc", &
                                  s_tor=0.25_dp, &
                                  M_pol=-1.0_dp, &
                                  N_tor=4.0_dp, &
                                  n_grid=n_grid, &
                                  max_n_fieldlines=n_fieldlines, &
                                  n_fieldlines=n_fieldlines, &
                                  eta_c=0.1_dp, &
                                  eta_t=0.3_dp)
    end function make_header

    subroutine fill_with_test_values(fieldlines)
        ! values whose decimal representation is not exact, plus extremes
        type(fieldline_for_precession_t), dimension(:), intent(inout) :: fieldlines
        integer :: idx, i_grid
        real(dp) :: sign_factor

        do idx = 1, size(fieldlines)
            do i_grid = 1, fieldlines(idx)%grid%n
                sign_factor = (-1.0_dp)**(idx + i_grid)
                fieldlines(idx)%grid%bounce_time_weighted(i_grid) = &
                    sign_factor*real(idx, dp)/3.0_dp*10.0_dp**(i_grid - 4)
                fieldlines(idx)%grid%I_j(i_grid) = acos(-1.0_dp)*real(i_grid, dp)/7.0_dp
                fieldlines(idx)%grid%radial_drift_weighted(i_grid) = &
                    sign_factor*tiny(1.0_dp)*real(idx + i_grid, dp)
                fieldlines(idx)%grid%poloidal_drift_weighted(i_grid) = &
                    -huge(1.0_dp)/real(idx*i_grid, dp)
            end do
        end do
    end subroutine fill_with_test_values

    subroutine write_fieldlines(filename, header, fieldlines, ids)
        character(len=*), intent(in) :: filename
        type(restart_header_t), intent(in) :: header
        type(fieldline_for_precession_t), dimension(:), intent(in) :: fieldlines
        integer, dimension(:), intent(in) :: ids
        integer :: i

        call write_restart_header(filename, header)
        do i = 1, size(ids)
            call append_restart_fieldline(filename, ids(i), fieldlines(ids(i))%grid)
        end do
    end subroutine write_fieldlines

    subroutine read_into_fresh(filename, header, fieldlines, is_complete)
        character(len=*), intent(in) :: filename
        type(restart_header_t), intent(in) :: header
        type(fieldline_for_precession_t), dimension(:), intent(inout) :: fieldlines
        logical, dimension(:), intent(out) :: is_complete

        fieldlines%grid = grid
        call read_restart_file(filename, header, fieldlines, is_complete)
    end subroutine read_into_fresh

    subroutine check_resume(label, filename, header, expected)
        ! appending the missing fieldlines after a read must give a complete file
        character(len=*), intent(in) :: label, filename
        type(restart_header_t), intent(in) :: header
        type(fieldline_for_precession_t), dimension(:), intent(in) :: expected

        type(fieldline_for_precession_t), dimension(size(expected)) :: fieldlines
        logical, dimension(size(expected)) :: is_complete
        integer :: idx

        call read_into_fresh(filename, header, fieldlines, is_complete)
        do idx = 1, size(expected)
            if (.not. is_complete(idx)) then
                call append_restart_fieldline(filename, idx, expected(idx)%grid)
            end if
        end do
        call read_into_fresh(filename, header, fieldlines, is_complete)
        call check_complete(label//" resumed", is_complete, &
                            spread(.true., 1, size(expected)))
        call check_values(label//" resumed", fieldlines, expected, is_complete)
    end subroutine check_resume

    subroutine check_complete(label, is_complete, expected_complete)
        character(len=*), intent(in) :: label
        logical, dimension(:), intent(in) :: is_complete, expected_complete

        if (any(is_complete .neqv. expected_complete)) then
            print *, "-------------------------------------------------------------"
            print *, "test_bounce_restart failed: ", label
            print *, "complete fieldlines found:    ", is_complete
            print *, "complete fieldlines expected: ", expected_complete
            test_failed = .true.
        end if
    end subroutine check_complete

    subroutine check_values(label, fieldlines, expected, is_complete)
        character(len=*), intent(in) :: label
        type(fieldline_for_precession_t), dimension(:), intent(in) :: fieldlines
        type(fieldline_for_precession_t), dimension(:), intent(in) :: expected
        logical, dimension(:), intent(in) :: is_complete
        integer :: idx

        do idx = 1, size(expected)
            if (.not. is_complete(idx)) cycle
            associate (found => fieldlines(idx)%grid, wanted => expected(idx)%grid)
                if (any(found%bounce_time_weighted /= wanted%bounce_time_weighted) &
                    .or. any(found%I_j /= wanted%I_j) &
                    .or. any(found%radial_drift_weighted &
                             /= wanted%radial_drift_weighted) &
                    .or. any(found%poloidal_drift_weighted &
                             /= wanted%poloidal_drift_weighted)) then
                    print *, "---------------------------------------------------------"
                    print *, "test_bounce_restart failed: ", label
                    print *, "values of fieldline ", idx, " not bit-exact"
                    test_failed = .true.
                end if
            end associate
        end do
    end subroutine check_values

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

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: file_unit

        open (newunit=file_unit, file=filename, status='old')
        close (file_unit, status='delete')
    end subroutine delete_file

end program test_bounce_restart
