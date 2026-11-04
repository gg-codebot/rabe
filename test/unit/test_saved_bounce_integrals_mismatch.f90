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

    type(bounce_integrals_header_t) :: header
    type(integration_grid_t) :: grid
    type(fieldline_for_precession_t), dimension(n_fieldlines) :: fieldlines
    logical, dimension(n_fieldlines) :: is_complete

    call set_integration_grids(0.3_dp, 0.1_dp, grid, 7)
    fieldlines%grid = grid

    header = bounce_integrals_header_t(git_hash="abc1234", &
                                       field_file="wout_test.nc", &
                                       s_tor=0.25_dp, &
                                       M_pol=-1.0_dp, &
                                       N_tor=4.0_dp, &
                                       n_grid=7, &
                                       max_n_fieldlines=n_fieldlines, &
                                       n_fieldlines=n_fieldlines, &
                                       eta_c=0.1_dp, &
                                       eta_t=0.3_dp)
    call write_bounce_integrals_header(test_file, header)

    header%s_tor = 0.5_dp
    call read_bounce_integrals(test_file, header, fieldlines, is_complete)

    print *, "test_saved_bounce_integrals_mismatch failed: mismatch not detected"
end program test_saved_bounce_integrals_mismatch
