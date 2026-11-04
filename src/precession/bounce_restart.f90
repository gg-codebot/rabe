module bounce_restart
    ! Restart file for the per-fieldline bounce integrals of electric_rabe.
    !
    ! Plain text: a header of "# key = value" lines describing the inputs the
    ! integrals depend on, followed by one fixed-width line per grid point:
    !   fieldline  t  eta  bounce_time_weighted  I_j
    !   radial_drift_weighted  poloidal_drift_weighted
    ! Fieldlines are appended as soon as they are computed. A fieldline counts
    ! as complete only if all grid%n of its lines are present and have full
    ! width, so a file cut off by a killed run is detected on reading.
    use constants, only: dp
    use grid_mod, only: integration_grid_t, fieldline_for_precession_t
    implicit none
    private

    public :: restart_header_t
    public :: write_restart_header
    public :: append_restart_fieldline
    public :: read_restart_file

    integer, parameter :: format_version = 1
    character(len=*), parameter :: line_format = '(I6, 6(1X, ES25.17E3))'
    integer, parameter :: line_width = 6 + 6*26
    character(len=*), parameter :: real_format = '(ES25.17E3)'
    integer, parameter :: max_line_length = 1024
    integer, parameter :: n_header_lines = 13
    integer, parameter :: git_hash_line = 3

    type :: restart_header_t
        character(len=:), allocatable :: git_hash
        character(len=:), allocatable :: field_file
        real(dp) :: s_tor
        real(dp) :: M_pol
        real(dp) :: N_tor
        integer :: n_grid
        integer :: max_n_fieldlines
        integer :: n_fieldlines
        real(dp) :: eta_c
        real(dp) :: eta_t
    end type restart_header_t

    interface
        function c_rename(old_name, new_name) bind(C, name="rename") result(ierr)
            use, intrinsic :: iso_c_binding, only: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: old_name
            character(kind=c_char), dimension(*), intent(in) :: new_name
            integer(c_int) :: ierr
        end function c_rename
    end interface

contains

    subroutine write_restart_header(filename, header)
        character(len=*), intent(in) :: filename
        type(restart_header_t), intent(in) :: header

        character(len=max_line_length), dimension(n_header_lines) :: lines
        integer :: file_unit, i

        lines = header_lines(header)
        open (newunit=file_unit, file=filename, status='replace', action='write')
        do i = 1, n_header_lines
            write (file_unit, '(A)') trim(lines(i))
        end do
        close (file_unit)
    end subroutine write_restart_header

    subroutine append_restart_fieldline(filename, fieldline_id, grid)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: fieldline_id
        type(integration_grid_t), intent(in) :: grid

        integer :: file_unit, i_grid

        open (newunit=file_unit, file=filename, status='old', &
              position='append', action='write')
        do i_grid = 1, grid%n
            write (file_unit, line_format) fieldline_id, &
                grid%t(i_grid), &
                grid%eta(i_grid), &
                grid%bounce_time_weighted(i_grid), &
                grid%I_j(i_grid), &
                grid%radial_drift_weighted(i_grid), &
                grid%poloidal_drift_weighted(i_grid)
        end do
        close (file_unit)
    end subroutine append_restart_fieldline

    subroutine read_restart_file(filename, header, fieldlines, is_complete)
        ! Loads all complete fieldlines from filename into fieldlines(:)%grid and
        ! rewrites the file with only those, so new fieldlines can be appended.
        ! Stops if the header does not match the current inputs.
        character(len=*), intent(in) :: filename
        type(restart_header_t), intent(in) :: header
        type(fieldline_for_precession_t), dimension(:), intent(inout) :: fieldlines
        logical, dimension(:), intent(out) :: is_complete

        character(len=max_line_length), dimension(n_header_lines) :: expected
        character(len=max_line_length) :: line
        real(dp), dimension(6) :: values
        integer :: file_unit, ios, i, fieldline_id, current_id, n_lines

        is_complete = .false.
        expected = header_lines(header)

        open (newunit=file_unit, file=filename, status='old', action='read')
        do i = 1, n_header_lines
            read (file_unit, '(A)', iostat=ios) line
            if (ios == 0 .and. trim(line) == trim(expected(i))) cycle
            if (ios /= 0 .or. is_cut_off(file_unit, line, expected(i))) then
                close (file_unit)
                print *, "Warning: incomplete header in ", filename, &
                    ", starting from scratch"
                call rewrite_restart_file(filename, header, fieldlines, is_complete)
                return
            end if
            if (i == git_hash_line) then
                print *, "Warning: restart file written by a different version"
                print *, "  file:    ", trim(line)
                print *, "  current: ", trim(expected(i))
                cycle
            end if
            print *, "Error: restart file ", filename, &
                " does not match the current input"
            print *, "  file:    ", trim(line)
            print *, "  current: ", trim(expected(i))
            print *, "Remove or rename it to start from scratch."
            error stop
        end do

        current_id = 0
        n_lines = 0
        do
            read (file_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) /= line_width) exit
            read (line, line_format, iostat=ios) fieldline_id, values
            if (ios /= 0) exit
            if (fieldline_id < 1 .or. fieldline_id > size(fieldlines)) exit

            if (fieldline_id /= current_id) then
                call mark_if_complete(current_id, n_lines, fieldlines, is_complete)
                current_id = fieldline_id
                n_lines = 0
                is_complete(current_id) = .false.
            end if
            n_lines = n_lines + 1
            if (n_lines > fieldlines(current_id)%grid%n) exit

            associate (grid => fieldlines(current_id)%grid)
                if (values(1) /= grid%t(n_lines) .or. &
                    values(2) /= grid%eta(n_lines)) then
                    print *, "Error: integration grid in restart file ", filename, &
                        " differs from the current one at fieldline ", current_id
                    error stop
                end if
                grid%bounce_time_weighted(n_lines) = values(3)
                grid%I_j(n_lines) = values(4)
                grid%radial_drift_weighted(n_lines) = values(5)
                grid%poloidal_drift_weighted(n_lines) = values(6)
            end associate
        end do
        call mark_if_complete(current_id, n_lines, fieldlines, is_complete)
        close (file_unit)

        call rewrite_restart_file(filename, header, fieldlines, is_complete)
    end subroutine read_restart_file

    function is_cut_off(file_unit, line, expected) result(cut_off)
        ! true if line is the beginning of expected and the last line of the file
        integer, intent(in) :: file_unit
        character(len=*), intent(in) :: line, expected
        logical :: cut_off

        character(len=1) :: next_line
        integer :: ios

        cut_off = .false.
        if (len_trim(line) >= len_trim(expected)) return
        if (index(expected, trim(line)) /= 1) return
        read (file_unit, '(A)', iostat=ios) next_line
        cut_off = is_iostat_end(ios)
    end function is_cut_off

    subroutine mark_if_complete(fieldline_id, n_lines, fieldlines, is_complete)
        integer, intent(in) :: fieldline_id, n_lines
        type(fieldline_for_precession_t), dimension(:), intent(in) :: fieldlines
        logical, dimension(:), intent(inout) :: is_complete

        if (fieldline_id < 1) return
        is_complete(fieldline_id) = (n_lines == fieldlines(fieldline_id)%grid%n)
    end subroutine mark_if_complete

    subroutine rewrite_restart_file(filename, header, fieldlines, is_complete)
        ! Writes the complete fieldlines to a temporary file first and renames it,
        ! so a crash during the rewrite leaves the original file intact.
        use, intrinsic :: iso_c_binding, only: c_null_char
        character(len=*), intent(in) :: filename
        type(restart_header_t), intent(in) :: header
        type(fieldline_for_precession_t), dimension(:), intent(in) :: fieldlines
        logical, dimension(:), intent(in) :: is_complete

        character(len=:), allocatable :: tmp_filename
        integer :: idx, file_unit

        tmp_filename = filename//".tmp"
        call write_restart_header(tmp_filename, header)
        do idx = 1, size(fieldlines)
            if (is_complete(idx)) then
                call append_restart_fieldline(tmp_filename, idx, fieldlines(idx)%grid)
            end if
        end do

        if (c_rename(tmp_filename//c_null_char, filename//c_null_char) /= 0) then
            ! rename does not replace an existing file on all platforms
            open (newunit=file_unit, file=filename, status='old')
            close (file_unit, status='delete')
            if (c_rename(tmp_filename//c_null_char, filename//c_null_char) /= 0) then
                print *, "Error: could not rename ", tmp_filename, " to ", filename
                error stop
            end if
        end if
    end subroutine rewrite_restart_file

    function header_lines(header) result(lines)
        type(restart_header_t), intent(in) :: header
        character(len=max_line_length), dimension(n_header_lines) :: lines

        lines(1) = "# electric_rabe restart file"
        lines(2) = "# format_version = "//int_to_str(format_version)
        lines(git_hash_line) = "# git_hash = "//header%git_hash
        lines(4) = "# field_file = "//trim(header%field_file)
        lines(5) = "# s_tor = "//real_to_str(header%s_tor)
        lines(6) = "# n_grid = "//int_to_str(header%n_grid)
        lines(7) = "# max_n_fieldlines = "//int_to_str(header%max_n_fieldlines)
        lines(8) = "# M_pol = "//real_to_str(header%M_pol)
        lines(9) = "# N_tor = "//real_to_str(header%N_tor)
        lines(10) = "# n_fieldlines = "//int_to_str(header%n_fieldlines)
        lines(11) = "# eta_c = "//real_to_str(header%eta_c)
        lines(12) = "# eta_t = "//real_to_str(header%eta_t)
        lines(13) = "# columns: fieldline  t  eta  bounce_time_weighted  I_j" &
                    //"  radial_drift_weighted  poloidal_drift_weighted"
    end function header_lines

    function int_to_str(value) result(str)
        integer, intent(in) :: value
        character(len=:), allocatable :: str
        character(len=32) :: buffer

        write (buffer, '(I0)') value
        str = trim(buffer)
    end function int_to_str

    function real_to_str(value) result(str)
        real(dp), intent(in) :: value
        character(len=:), allocatable :: str
        character(len=32) :: buffer

        write (buffer, real_format) value
        str = trim(adjustl(buffer))
    end function real_to_str

end module bounce_restart
