program test_read_namelist
    use constants, only: dp
    use utils, only: not_same
    use read_file, only: read_namelist
    use read_file, only: field_file, &
                         field_type, &
                         M_pol, &
                         N_tor, &
                         s_tor, &
                         s_tor_min, &
                         s_tor_max, &
                         n_s_tor, &
                         sign_sqrtg, &
                         max_n_fieldlines, &
                         should_calc_shaing_callen, &
                         n_eta

    implicit none

    integer :: this
    real(dp), parameter :: abstol = 0.0_dp, reltol = 1e-14
    character(len=*), parameter :: test_file = "rabe.in"

    logical :: test_failed

    character(len=*), parameter :: test_field_file = "test.nc"
    character(len=*), parameter :: test_field_type = "booz_xform"
    real(dp), parameter :: test_M_pol = -1.0_dp
    real(dp), parameter :: test_N_tor = 2.0_dp
    real(dp), dimension(:), allocatable :: test_s_tor
    real(dp), parameter :: test_sign_sqrtg = -1.0_dp
    integer, parameter :: test_max_n_fieldlines = 10
    logical, parameter :: test_should_calc_shaing_callen = .true.
    integer, parameter :: test_n_eta = 11

    test_failed = .false.

    if (allocated(test_s_tor)) deallocate (test_s_tor)
    allocate (test_s_tor, source=[0.25_dp, 0.5_dp])
    call write_test_file(test_file, &
                         test_field_file=test_field_file, &
                         test_field_type=test_field_type, &
                         test_M_pol=test_M_pol, &
                         test_N_tor=test_N_tor, &
                         test_s_tor=test_s_tor, &
                         test_sign_sqrtg=test_sign_sqrtg, &
                         test_max_n_fieldlines=test_max_n_fieldlines, &
                         test_should_calc_shaing_callen= &
                         test_should_calc_shaing_callen, &
                         test_n_eta=test_n_eta)
    call read_namelist(test_file)

    if (field_file /= test_field_file) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: bc_filename"
        print *, "found: ", field_file
        print *, "expected: ", test_field_file
        test_failed = .true.
    end if
    if (field_type /= test_field_type) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: field_type"
        print *, "found: ", trim(field_type)
        print *, "expected: ", test_field_type
        test_failed = .true.
    end if
    if (not_same(M_pol, &
                 test_M_pol, &
                 reltol_in=reltol, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: M_pol"
        print *, "found: ", M_pol
        print *, "expected: ", test_M_pol
        test_failed = .true.
    end if
    if (not_same(N_tor, &
                 test_N_tor, &
                 reltol_in=reltol, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: N_tor"
        print *, "found: ", N_tor
        print *, "expected: ", test_N_tor
        test_failed = .true.
    end if
    if (not_same(s_tor, &
                 test_s_tor, &
                 reltol_in=reltol, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: s_tor"
        print *, "found: ", s_tor
        print *, "expected: ", test_s_tor
        test_failed = .true.
    end if
    if (not_same(sign_sqrtg, &
                 test_sign_sqrtg, &
                 reltol_in=reltol, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: sign_sqrtg"
        print *, "found: ", sign_sqrtg
        print *, "expected: ", test_sign_sqrtg
        test_failed = .true.
    end if
    if (max_n_fieldlines /= test_max_n_fieldlines) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: max_n_fieldlines"
        print *, "found: ", max_n_fieldlines
        print *, "expected: ", test_max_n_fieldlines
        test_failed = .true.
    end if
    if (should_calc_shaing_callen .neqv. test_should_calc_shaing_callen) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: should_calc_shaing_callen"
        print *, "found: ", should_calc_shaing_callen
        print *, "expected: ", test_should_calc_shaing_callen
        test_failed = .true.
    end if
    if (n_eta /= test_n_eta) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: n_eta"
        print *, "found: ", n_eta
        print *, "expected: ", test_n_eta
        test_failed = .true.
    end if

    call remove_test_file(test_file)

    if (allocated(test_s_tor)) deallocate (test_s_tor)
    allocate (test_s_tor, source=[0.314_dp])
    call write_test_file(test_file, &
                         test_field_file=test_field_file, &
                         test_M_pol=test_M_pol, &
                         test_N_tor=test_N_tor, &
                         test_s_tor=test_s_tor, &
                         test_sign_sqrtg=test_sign_sqrtg, &
                         test_max_n_fieldlines=test_max_n_fieldlines, &
                         test_should_calc_shaing_callen= &
                         test_should_calc_shaing_callen, &
                         test_n_eta=test_n_eta)
    call read_namelist(test_file)

    if (not_same(s_tor, &
                 test_s_tor, &
                 reltol_in=reltol, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: s_tor"
        print *, "found: ", s_tor
        print *, "expected: ", test_s_tor
        test_failed = .true.
    end if

    call remove_test_file(test_file)

    ! Test range specification: s_tor_min/s_tor_max/n_s_tor
    call write_test_file(test_file, &
                         test_field_file=test_field_file, &
                         test_M_pol=test_M_pol, &
                         test_N_tor=test_N_tor, &
                         test_s_tor_min=0.2_dp, &
                         test_s_tor_max=0.8_dp, &
                         test_n_s_tor=4, &
                         test_sign_sqrtg=test_sign_sqrtg, &
                         test_max_n_fieldlines=test_max_n_fieldlines, &
                         test_should_calc_shaing_callen= &
                         test_should_calc_shaing_callen, &
                         test_n_eta=test_n_eta)
    call read_namelist(test_file)

    if (allocated(test_s_tor)) deallocate (test_s_tor)
    allocate (test_s_tor, source=[0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp])
    if (size(s_tor) /= size(test_s_tor)) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: s_tor size from range"
        print *, "found: ", size(s_tor)
        print *, "expected: ", size(test_s_tor)
        test_failed = .true.
    else if (not_same(s_tor, &
                      test_s_tor, &
                      reltol_in=reltol, &
                      abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_read_namelist failed: s_tor values from range"
        print *, "found: ", s_tor
        print *, "expected: ", test_s_tor
        test_failed = .true.
    end if

    call remove_test_file(test_file)

    if (test_failed) error stop

contains

    subroutine write_test_file(filename, &
                               test_field_file, &
                               test_M_pol, &
                               test_N_tor, &
                               test_s_tor, &
                               test_s_tor_min, &
                               test_s_tor_max, &
                               test_n_s_tor, &
                               test_sign_sqrtg, &
                               test_max_n_fieldlines, &
                               test_should_calc_shaing_callen, &
                               test_n_eta, &
                               test_field_type)

        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: test_field_file
        real(dp), intent(in), optional :: test_M_pol
        real(dp), intent(in), optional :: test_N_tor
        real(dp), intent(in), dimension(:), optional :: test_s_tor
        real(dp), intent(in), optional :: test_s_tor_min
        real(dp), intent(in), optional :: test_s_tor_max
        integer, intent(in), optional :: test_n_s_tor
        real(dp), intent(in), optional :: test_sign_sqrtg
        integer, intent(in), optional :: test_max_n_fieldlines
        logical, intent(in), optional :: test_should_calc_shaing_callen
        integer, intent(in), optional :: test_n_eta
        character(len=*), intent(in), optional :: test_field_type

        integer :: unit

        open (newunit=unit, file=filename, status="replace")
        write (unit, "(A)") "&rabe_config"
        if (present(test_field_file)) then
            write (unit, "(A,A,A)") "field_file = '", test_field_file, "',"
        end if
        if (present(test_field_type)) then
            write (unit, "(A,A,A)") "field_type = '", test_field_type, "',"
        end if
        if (present(test_M_pol)) then
            write (unit, "(A,F4.1,A)") "M_pol = ", test_M_pol, ","
        end if
        if (present(test_N_tor)) then
            write (unit, "(A,F4.1,A)") "N_tor = ", test_N_tor, ","
        end if
        if (present(test_s_tor)) then
            write (unit, "(A)", advance="no") "s_tor = "
            write (unit, "(*(G0.3,1X))", advance="no") test_s_tor
            write (unit, "(A)") ","
        end if
        if (present(test_s_tor_min)) then
            write (unit, "(A,G0.3,A)") "s_tor_min = ", test_s_tor_min, ","
        end if
        if (present(test_s_tor_max)) then
            write (unit, "(A,G0.3,A)") "s_tor_max = ", test_s_tor_max, ","
        end if
        if (present(test_n_s_tor)) then
            write (unit, "(A,I6,A)") "n_s_tor = ", test_n_s_tor, ","
        end if
        if (present(test_sign_sqrtg)) then
            write (unit, "(A,F4.1,A)") "sign_sqrtg = ", test_sign_sqrtg, ","
        end if
        if (present(test_max_n_fieldlines)) then
            write (unit, "(A,I6,A)") "max_n_fieldlines = ", test_max_n_fieldlines, ","
        end if
        if (present(test_should_calc_shaing_callen)) then
            write (unit, "(A,L,A)") "should_calc_shaing_callen = ", &
                test_should_calc_shaing_callen, ","
        end if
        if (present(test_n_eta)) then
            write (unit, "(A,I6,A)") "n_eta = ", test_n_eta
        end if
        write (unit, "(A)") "/"
        close (unit)

    end subroutine write_test_file

    subroutine remove_test_file(filename)
        character(len=*), intent(in) :: filename

        logical :: file_exists
        integer :: unit

        inquire (file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, "error: file not found: ", filename
            stop
        end if

        open (newunit=unit, file=filename, status="old")
        close (unit, status="delete")

    end subroutine remove_test_file

end program test_read_namelist
