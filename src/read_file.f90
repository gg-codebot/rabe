module read_file
    use constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none

    character(len=100), public, protected :: field_file
    character(len=20), public, protected :: field_type
    real(dp), public, protected :: M_pol
    real(dp), public, protected :: N_tor
    real(dp), dimension(:), allocatable, public, protected :: s_tor
    real(dp), public, protected :: ds_dr ![1/cm]
    real(dp), public, protected :: sign_sqrtg
    integer, public, protected :: max_n_fieldlines
    logical, public, protected :: should_calc_shaing_callen
    integer, public, protected :: n_eta
    logical, public, protected :: unsafe_mode

    real(dp), public, protected :: s_tor_min
    real(dp), public, protected :: s_tor_max
    integer, public, protected :: n_s_tor

    namelist /rabe_config/ &
        field_file, &
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
        n_eta, &
        unsafe_mode

contains

    subroutine read_namelist(filename)
        use utils, only: linspace
        implicit none
        character(len=*), intent(in) :: filename
        integer :: ios, unit
        logical :: file_exists

        real(dp), dimension(:), allocatable :: s_tor_temp
        integer, parameter :: n_stor_max = 1000
        integer :: n_list
        real(dp) :: nan_value
        logical :: has_list, has_range

        inquire (file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, "Error in read_namelist: file not found: ", filename
            error stop
        end if

        open (newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print *, "Error in read_namelist: cannot open file: ", filename
            error stop
        end if

        nan_value = ieee_value(nan_value, ieee_quiet_nan)

        ! Default values
        field_type = 'vmec_nc'
        unsafe_mode = .false.
        should_calc_shaing_callen = .false.
        n_eta = 100
        s_tor_min = nan_value
        s_tor_max = nan_value
        n_s_tor = 0
        if (allocated(s_tor)) deallocate (s_tor)
        allocate (s_tor(n_stor_max))
        s_tor = nan_value

        read (unit, nml=rabe_config, iostat=ios)
        if (ios /= 0) then
            print *, "Error in read_namelist:"
            print *, "iostat = ", ios
            error stop
        end if
        close (unit)

        n_list = count(.not. ieee_is_nan(s_tor))
        has_list = n_list > 0
        has_range = (.not. ieee_is_nan(s_tor_min)) .and. &
                    (.not. ieee_is_nan(s_tor_max)) .and. &
                    (n_s_tor > 0)

        if (has_list .and. has_range) then
            print *, "Error in read_namelist: both s_tor list and s_tor range " &
                //"(s_tor_min, s_tor_max, n_s_tor) are specified. Use one or the other."
            error stop
        end if

        if (.not. has_list .and. .not. has_range) then
            print *, "Error in read_namelist: no s_tor values provided. " &
                //"Set either s_tor or (s_tor_min, s_tor_max, n_s_tor)."
            error stop
        end if

        if (has_range) then
            deallocate (s_tor)
            allocate (s_tor(n_s_tor))
            call linspace(s_tor_min, s_tor_max, n_s_tor, s_tor)
        else
            if (n_list >= n_stor_max) then
                print *, "Error in read_namelist: too many s_tor values ", &
                    "provided in namelist"
                error stop
            end if
            if (allocated(s_tor_temp)) deallocate (s_tor_temp)
            allocate (s_tor_temp(n_list))
            s_tor_temp = pack(s_tor,.not. ieee_is_nan(s_tor))
            deallocate (s_tor)
            allocate (s_tor(n_list))
            s_tor = s_tor_temp
            deallocate (s_tor_temp)
        end if

        call check_if_valid_namelist()

    end subroutine read_namelist

    subroutine check_if_valid_namelist()
        use make_fieldline, only: is_not_integer
        use utils, only: not_same

        real(dp), parameter :: tol = 1e-15
        logical :: is_valid

        is_valid = .true.

        if (len(trim(field_file)) == 0) then
            print *, "field_file is empty!"
            is_valid = .false.
        end if
        if (len(trim(field_type)) == 0 .or. &
            (trim(field_type) /= 'vmec_nc' .and. &
             trim(field_type) /= 'booz_xform' .and. &
             trim(field_type) /= 'chartmap')) then
            print *, "field_type must be 'vmec_nc', 'booz_xform', or 'chartmap'"
            is_valid = .false.
        end if
        if (ieee_is_nan(M_pol)) then
            print *, "M_pol is NaN!"
            is_valid = .false.
        end if
        if (ieee_is_nan(N_tor)) then
            print *, "N_tor is NaN!"
            is_valid = .false.
        end if
        if (any(ieee_is_nan(s_tor))) then
            print *, "s_tor is NaN!"
            is_valid = .false.
        end if
        if (ieee_is_nan(sign_sqrtg)) then
            print *, "sign_sqrtg is NaN!"
            is_valid = .false.
        end if

        if (is_not_integer(M_pol, tol)) then
            print *, "M_pol must be integer"
            is_valid = .false.
        end if
        if (is_not_integer(N_tor, tol)) then
            print *, "N_tor must be integer"
            is_valid = .false.
        end if
        if (not_same(sign_sqrtg, 1.0_dp, reltol_in=0.0_dp, abstol_in=tol) .and. &
            not_same(sign_sqrtg, -1.0_dp, reltol_in=0.0_dp, abstol_in=tol)) then
            print *, "sign_sqrtg must be +-1"
            is_valid = .false.
        end if
        if (any(s_tor <= 0.0_dp) .or. any(s_tor >= 1.0_dp)) then
            print *, "s_tor must be in ]0.0, 1.0["
            is_valid = .false.
        end if
        if (n_s_tor > 0) then
            if (s_tor_min >= s_tor_max) then
                print *, "s_tor_min must be less than s_tor_max"
                is_valid = .false.
            end if
            if (n_s_tor < 2) then
                print *, "n_s_tor must be at least 2"
                is_valid = .false.
            end if
        end if
        if (max_n_fieldlines <= 1) then
            print *, "max_n_fieldlines must be bigger than 1"
            is_valid = .false.
        end if
        if (should_calc_shaing_callen .and. n_eta <= 3) then
            print *, "n_eta must be bigger than 3"
        end if

        if (.not. is_valid) error stop

    end subroutine check_if_valid_namelist

end module read_file
