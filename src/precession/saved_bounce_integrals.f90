module saved_bounce_integrals
    ! Per-fieldline bounce integrals of electric_rabe saved as plain text:
    ! "# key = value" header lines with the inputs they depend on, then one
    ! fixed-width line per fieldline and grid point. Fieldlines are appended as
    ! soon as they are computed. On reading, a fieldline counts only if all of
    ! its lines are present at full width, so a file cut off by a killed run is
    ! safe to resume from. Only a clean build of the same commit reuses a file.
    use constants, only: dp
    use grid_mod, only: integration_grid_t, fieldline_for_precession_t
    implicit none
    private

    public :: bounce_integrals_header_t
    public :: write_bounce_integrals_header
    public :: append_bounce_integrals
    public :: read_bounce_integrals

    character(len=*), parameter :: line_format = '(I6, 6(1X, ES25.17E3))'
    integer, parameter :: line_width = 6 + 6*26
    integer, parameter :: n_header_lines = 12
    integer, parameter :: git_hash_line = 2

    type :: bounce_integrals_header_t
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
    end type bounce_integrals_header_t

contains

    subroutine write_bounce_integrals_header(filename, header)
        character(len=*), intent(in) :: filename
        type(bounce_integrals_header_t), intent(in) :: header

        character(len=1024), dimension(n_header_lines) :: lines
        integer :: file_unit, i

        lines = header_lines(header)
        open (newunit=file_unit, file=filename, status='replace', action='write')
        do i = 1, n_header_lines
            write (file_unit, '(A)') trim(lines(i))
        end do
        close (file_unit)
    end subroutine write_bounce_integrals_header

    subroutine append_bounce_integrals(filename, fieldline_id, grid)
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
    end subroutine append_bounce_integrals

    subroutine read_bounce_integrals(filename, header, fieldlines, is_complete)
        ! Loads all complete fieldlines into fieldlines(:)%grid and writes the
        ! file back with only those, so new fieldlines append to a clean end.
        ! Stops if the header does not match the current input and git hash,
        ! or if the file was written by a dirty build.
        character(len=*), intent(in) :: filename
        type(bounce_integrals_header_t), intent(in) :: header
        type(fieldline_for_precession_t), dimension(:), intent(inout) :: fieldlines
        logical, dimension(:), intent(out) :: is_complete

        character(len=1024), dimension(n_header_lines) :: expected
        character(len=1024) :: line
        real(dp), dimension(6) :: values
        integer :: file_unit, ios, i, id, line_id, n_lines

        expected = header_lines(header)
        open (newunit=file_unit, file=filename, status='old', action='read')
        do i = 1, n_header_lines
            read (file_unit, '(A)', iostat=ios) line
            if (ios /= 0) line = ""
            if (i == git_hash_line .and. index(line, "-dirty") > 0) then
                print *, "Error: ", filename, " was written by a dirty build"
                error stop
            end if
            if (line == expected(i)) cycle
            print *, "Error: ", filename, " does not match the current input"
            print *, "  file:    ", trim(line)
            print *, "  current: ", trim(expected(i))
            error stop
        end do

        is_complete = .false.
        id = 0
        n_lines = 0
        do
            read (file_unit, '(A)', iostat=ios) line
            if (ios /= 0 .or. len_trim(line) /= line_width) exit
            read (line, line_format) line_id, values
            if (line_id /= id) then
                if (id > 0) is_complete(id) = n_lines == fieldlines(id)%grid%n
                id = line_id
                n_lines = 0
            end if
            n_lines = n_lines + 1
            fieldlines(id)%grid%bounce_time_weighted(n_lines) = values(3)
            fieldlines(id)%grid%I_j(n_lines) = values(4)
            fieldlines(id)%grid%radial_drift_weighted(n_lines) = values(5)
            fieldlines(id)%grid%poloidal_drift_weighted(n_lines) = values(6)
        end do
        if (id > 0) is_complete(id) = n_lines == fieldlines(id)%grid%n
        close (file_unit)

        call write_bounce_integrals_header(filename, header)
        do id = 1, size(fieldlines)
            if (is_complete(id)) then
                call append_bounce_integrals(filename, id, fieldlines(id)%grid)
            end if
        end do
    end subroutine read_bounce_integrals

    function header_lines(header) result(lines)
        type(bounce_integrals_header_t), intent(in) :: header
        character(len=1024), dimension(n_header_lines) :: lines

        lines(1) = "# electric_rabe saved bounce integrals"
        lines(git_hash_line) = "# git_hash = "//header%git_hash
        lines(3) = "# field_file = "//header%field_file
        write (lines(4), '(A, ES26.17E3)') "# s_tor =", header%s_tor
        write (lines(5), '(A, I0)') "# n_grid = ", header%n_grid
        write (lines(6), '(A, I0)') "# max_n_fieldlines = ", header%max_n_fieldlines
        write (lines(7), '(A, ES26.17E3)') "# M_pol =", header%M_pol
        write (lines(8), '(A, ES26.17E3)') "# N_tor =", header%N_tor
        write (lines(9), '(A, I0)') "# n_fieldlines = ", header%n_fieldlines
        write (lines(10), '(A, ES26.17E3)') "# eta_c =", header%eta_c
        write (lines(11), '(A, ES26.17E3)') "# eta_t =", header%eta_t
        lines(12) = "# columns: fieldline  t  eta  bounce_time_weighted  I_j" &
                    //"  radial_drift_weighted  poloidal_drift_weighted"
    end function header_lines

end module saved_bounce_integrals
