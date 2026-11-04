program test_bounce_restart_mismatch
    ! Reading a restart file written for different inputs must stop.
    ! The test passes if the expected error message is printed.
    use constants, only: dp
    use grid_mod, only: integration_grid_t, fieldline_for_precession_t
    use grid_mod, only: set_integration_grids
    use bounce_restart, only: restart_header_t, write_restart_header
    use bounce_restart, only: read_restart_file

    implicit none

    character(len=*), parameter :: test_file = "test_bounce_restart_mismatch.dat"
    integer, parameter :: n_fieldlines = 3

    type(restart_header_t) :: header
    type(integration_grid_t) :: grid
    type(fieldline_for_precession_t), dimension(n_fieldlines) :: fieldlines
    logical, dimension(n_fieldlines) :: is_complete

    call set_integration_grids(0.3_dp, 0.1_dp, grid, 7)
    fieldlines%grid = grid

    header = restart_header_t(git_hash="abc1234", &
                              field_file="wout_test.nc", &
                              s_tor=0.25_dp, &
                              M_pol=-1.0_dp, &
                              N_tor=4.0_dp, &
                              n_grid=7, &
                              max_n_fieldlines=n_fieldlines, &
                              n_fieldlines=n_fieldlines, &
                              eta_c=0.1_dp, &
                              eta_t=0.3_dp)
    call write_restart_header(test_file, header)

    header%s_tor = 0.5_dp
    call read_restart_file(test_file, header, fieldlines, is_complete)

    print *, "test_bounce_restart_mismatch failed: mismatch not detected"
end program test_bounce_restart_mismatch
