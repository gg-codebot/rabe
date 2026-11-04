program test_saved_bounce_integrals_mismatch
    ! Reading saved bounce integrals written for different inputs must stop.
    ! The test passes if the expected error message is printed.
    use constants, only: dp
    use grid_mod, only: integration_grid_t, fieldline_for_precession_t
    use grid_mod, only: set_integration_grids
    use saved_bounce_integrals, only: bounce_integrals_header_t
    use saved_bounce_integrals, only: write_bounce_integrals_header
    use saved_bounce_integrals, only: read_bounce_integrals

    implicit none

    character(len=*), parameter :: test_file = &
                                   "test_saved_bounce_integrals_mismatch.dat"
    integer, parameter :: n_fieldlines = 3
    integer, parameter :: n_grid = 7
    real(dp), parameter :: eta_c = 0.1_dp
    real(dp), parameter :: eta_t = 0.3_dp
    real(dp), parameter :: s_tor = 0.25_dp
    real(dp), parameter :: M_pol = -1.0_dp
    real(dp), parameter :: N_tor = 4.0_dp
    character(len=*), parameter :: field_file = "wout_test.nc"
    real(dp), parameter :: other_s_tor = 0.5_dp
    character(len=*), parameter :: git_hash = "abc1234"

    type(bounce_integrals_header_t) :: header
    type(integration_grid_t) :: grid
    type(fieldline_for_precession_t), dimension(n_fieldlines) :: fieldlines
    logical, dimension(n_fieldlines) :: is_complete

    call set_integration_grids(eta_t, eta_c, grid, n_grid)
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
    call write_bounce_integrals_header(test_file, header)

    header%s_tor = other_s_tor
    call read_bounce_integrals(test_file, header, fieldlines, is_complete)

    print *, "test_saved_bounce_integrals_mismatch failed: mismatch not detected"
end program test_saved_bounce_integrals_mismatch
