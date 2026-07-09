program main
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use constants, only: dp
    use netcdf_mod, only: netcdf_t
    use read_file, only: read_namelist, &
                         field_file, &
                         field_type, &
                         M_pol, &
                         N_tor, &
                         s_tor, &
                         sign_sqrtg, &
                         max_n_fieldlines, &
                         should_calc_shaing_callen, &
                         n_eta, &
                         unsafe_mode
    use boozer_field, only: boozer_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_of_fieldlines
    use coefficients, only: calc_nu_star_crit, &
                            calc_finite_boundary_layer_correction, &
                            calc_gradient_scaling_factor_r_eff, &
                            calc_offset_coefficients
    use shaing_callen_mod, only: calc_lambda_LC, get_non_omnigenous_remainder
    use error_handling, only: set_unsafe_mode
    use error_handling, only: reset_failed_check_counter, did_fail_any_sanity_check
    use git_version, only: git_hash

    implicit none

    character(len=*), parameter :: input_file = "rabe.in"
    character(len=*), parameter :: output_file = "rabe.nc"
    character(len=*), parameter :: dat_file = "rabe.dat"

    type(boozer_field_t) :: field

    integer :: n_stor
    integer :: this

    real(dp) :: R ![m]
    real(dp) :: dr_dAtheta ![rad/Tm/]
    real(dp) :: iota
    real(dp) :: nfp

    type(flock_of_fieldlines_t) :: flock
    logical :: too_strong_violation

    real(dp), dimension(:), allocatable :: Lambda_A, Lambda_B
    real(dp), dimension(:), allocatable :: nu_star_crit
    real(dp), dimension(:), allocatable :: Lambda_S
    integer, dimension(:), allocatable :: split_maxima
    real(dp) :: helical_factor

    real(dp) :: trapped_fraction
    real(dp), dimension(:), allocatable :: lambda_LC, remainder

    real(dp) :: nan_value

    type(netcdf_t) :: nc_output
    character(len=*), parameter :: dim_name = "surface"
    character(len=1024) :: description

    call read_namelist(input_file)
    call set_unsafe_mode(unsafe_mode)
    nan_value = ieee_value(nan_value, ieee_quiet_nan)

    n_stor = size(s_tor)
    allocate (Lambda_A(n_stor))
    allocate (Lambda_B(n_stor))
    allocate (nu_star_crit(n_stor))
    allocate (Lambda_S(n_stor))
    allocate (split_maxima(n_stor))
    if (should_calc_shaing_callen) then
        allocate (lambda_LC(n_stor))
        allocate (remainder(n_stor))
    end if

    select case (trim(field_type))
    case ('vmec_nc')
        call field%init_from_vmec(field_file, grid_refinement=6)
    case ('booz_xform')
        call field%init_from_boozmn(field_file)
    case ('chartmap')
        call field%init_from_chartmap(field_file)
    case default
        error stop "main: unknown field_type"
    end select
    do this = 1, n_stor
        call reset_failed_check_counter()
        call field%fix_to_surface(s_tor(this))
        call field%get_iota(s_tor(this), iota)
        nfp = field%nfp

        call make_flock_of_fieldlines(flock, max_n_fieldlines, iota, &
                                      field, M_pol, N_tor, nfp, &
                                      split_maxima(this))

        dr_dAtheta = calc_gradient_scaling_factor_r_eff(flock, field%psi_tor_edge, &
                                                        nint(sign_sqrtg))
        R = field%R
        call calc_offset_coefficients(flock, R, dr_dAtheta, &
                                      Lambda_A(this), Lambda_B(this))
        nu_star_crit(this) = calc_nu_star_crit(flock, R)
        Lambda_S(this) = calc_finite_boundary_layer_correction(flock, field, &
                                                               R, &
                                                               dr_dAtheta)

        print *, "s_tor: ", s_tor(this)
        if (should_calc_shaing_callen) then
            lambda_LC(this) = calc_lambda_LC(flock, field, n_eta, dr_dAtheta)
            remainder(this) = get_non_omnigenous_remainder(flock, field, &
                                                           n_eta, dr_dAtheta)
            print *, "omnigenous lambda_LC_bB: ", lambda_LC(this)
            print *, "non-omnigneous remainder: ", remainder(this)
        end if
        print *, "1/sqrt(nu_star) factor: ", Lambda_A(this)
        print *, "1/nu_star factor: ", Lambda_B(this)
        print *, "Lambda_S: ", Lambda_S(this)

        if (allocated(flock%fieldlines)) deallocate (flock%fieldlines)

        if (unsafe_mode .and. did_fail_any_sanity_check()) then
            lambda_A(this) = nan_value
            lambda_B(this) = nan_value
            nu_star_crit(this) = nan_value
            Lambda_S(this) = nan_value
            if (should_calc_shaing_callen) then
                lambda_LC(this) = nan_value
                remainder(this) = nan_value
            end if
        end if

    end do

    call nc_output%create(output_file)
    call nc_output%add_global_attribute("title", &
                                        "asymptotic bootstrap coefficient lambda_bB")
    call nc_output%add_global_attribute("definition", &
                          "lambda^{off}_bB = Lambda_A/sqrt(nu_star) + Lambda_B/nu_star")
    call nc_output%add_global_attribute("git_hash", git_hash)
    call nc_output%def_dim(dim_name, n_stor)
    call nc_output%add_real_1d("s_tor", dim_name)
    call nc_output%add_attr("s_tor", "long_name", &
                            "normalized toroidal flux label")
    call nc_output%add_attr("s_tor", "unit", &
                            "[1]")
    call nc_output%add_real_1d("Lambda_A", dim_name)
    call nc_output%add_attr("Lambda_A", "long_name", &
                            "1/sqrt(nu_star) factor")
    call nc_output%add_attr("Lambda_A", "unit", &
                            "[1]")
    call nc_output%add_real_1d("Lambda_B", dim_name)
    call nc_output%add_attr("Lambda_B", "long_name", &
                            "1/nu_star factor")
    call nc_output%add_attr("Lambda_B", "unit", &
                            "[1]")
    call nc_output%add_real_1d("nu_star_crit", dim_name)
    call nc_output%add_attr("nu_star_crit", "long_name", &
                          "lower collisionality limit for validity of asymptotic model")
    call nc_output%add_attr("nu_star_crit", "unit", &
                            "[1]")
    call nc_output%add_real_1d("Lambda_S", dim_name)
    call nc_output%add_attr("Lambda_S", "long_name", &
                      "sqrt(nu_star) factor accounting for finite boundary layer width")
    call nc_output%add_attr("Lambda_S", "unit", &
                            "[1]")
    call nc_output%add_int_1d("split_maxima", dim_name)
    call nc_output%add_attr("split_maxima", "long_name", &
                            "1 if violation of omnigeneity is too strong, 0 otherwise")
    call nc_output%add_real("R")
    call nc_output%add_attr("R", "long_name", &
                            "major radius")
    call nc_output%add_attr("R", "unit", "[m]")
    write (description, "(A,A,A)") "defines the reference length scale to convert ", &
        "to dimensionless quantities i.e. defines coeffients in respect to ", &
        "nu_star = pi*R/(2*mean_free_path) = pi*R*deflection_frequency/particle_speed"
    call nc_output%add_attr("R", "definition", description)
    if (should_calc_shaing_callen) then
        call nc_output%add_real_1d("lambda_LC_bB", dim_name)
        call nc_output%add_attr("lambda_LC_bB", "long_name", &
                                "omnigenous Shaing-Callen coefficient")
        call nc_output%add_attr("lambda_LC_bB", "unit", "[1]")
        call nc_output%add_real_1d("remainder", dim_name)
        call nc_output%add_attr("remainder", "long_name", &
                                "non-omnigenous remainder of Shaing-Callen coefficient")
        call nc_output%add_attr("remainder", "unit", "[1]")

        call nc_output%write_real_1d("lambda_LC_bB", lambda_LC)
        call nc_output%write_real_1d("remainder", remainder)
    end if
    call nc_output%write_real_1d("Lambda_A", Lambda_A)
    call nc_output%write_real_1d("Lambda_B", Lambda_B)
    call nc_output%write_real_1d("nu_star_crit", nu_star_crit)
    call nc_output%write_real_1d("Lambda_S", Lambda_S)
    call nc_output%write_int_1d("split_maxima", split_maxima)
    call nc_output%write_real_1d("s_tor", s_tor)
    call nc_output%write_real("R", field%R)
    call nc_output%close()

    if (should_calc_shaing_callen) then
        call write_dat_output(dat_file, git_hash, field%R, n_stor, &
                              s_tor, Lambda_A, Lambda_B, &
                              nu_star_crit, Lambda_S, split_maxima, &
                              lambda_LC=lambda_LC, remainder=remainder)
    else
        call write_dat_output(dat_file, git_hash, field%R, n_stor, &
                              s_tor, Lambda_A, Lambda_B, &
                              nu_star_crit, Lambda_S, split_maxima)
    end if

    if (allocated(Lambda_A)) deallocate (Lambda_A)
    if (allocated(Lambda_B)) deallocate (Lambda_B)
    if (allocated(nu_star_crit)) deallocate (nu_star_crit)
    if (allocated(Lambda_S)) deallocate (Lambda_S)
    if (allocated(lambda_LC)) deallocate (lambda_LC)
    if (allocated(remainder)) deallocate (remainder)

contains

    subroutine write_dat_output(filename, git_hash, R, n, s_tor, &
                                Lambda_A, Lambda_B, nu_star_crit, &
                                Lambda_S, split_maxima, &
                                lambda_LC, remainder)
        use constants, only: dp
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: git_hash
        real(dp), intent(in) :: R
        integer, intent(in) :: n
        real(dp), intent(in) :: s_tor(n), Lambda_A(n), Lambda_B(n)
        real(dp), intent(in) :: nu_star_crit(n), Lambda_S(n)
        integer, intent(in) :: split_maxima(n)
        real(dp), optional, intent(in) :: lambda_LC(n), remainder(n)

        integer :: u, i
        logical :: with_sc

        with_sc = present(lambda_LC)

        open (newunit=u, file=filename, status="replace", action="write")
        write (u, "(A)") "# asymptotic bootstrap coefficient lambda_bB"
        write (u, "(A)") &
            "# lambda^{off}_bB = Lambda_A/sqrt(nu_star) + Lambda_B/nu_star"
        write (u, "(A,A)") "# git_hash: ", trim(git_hash)
        write (u, "(A,ES23.15)") "# R [m]: ", R
        if (with_sc) then
            write (u, "(A)") &
                "# s_tor  Lambda_A  Lambda_B  nu_star_crit" &
                //"  Lambda_S  split_maxima  lambda_LC_bB  remainder"
        else
            write (u, "(A)") &
                "# s_tor  Lambda_A  Lambda_B  nu_star_crit" &
                //"  Lambda_S  split_maxima"
        end if
        do i = 1, n
            if (with_sc) then
                write (u, "(ES23.15,4(1X,ES23.15),1X,I2,2(1X,ES23.15))") &
                    s_tor(i), Lambda_A(i), Lambda_B(i), &
                    nu_star_crit(i), Lambda_S(i), split_maxima(i), &
                    lambda_LC(i), remainder(i)
            else
                write (u, "(ES23.15,4(1X,ES23.15),1X,I2)") &
                    s_tor(i), Lambda_A(i), Lambda_B(i), &
                    nu_star_crit(i), Lambda_S(i), split_maxima(i)
            end if
        end do
        close (u)
    end subroutine write_dat_output

end program main
