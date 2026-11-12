program electric_rabe
    use constants, only: dp, pi, machine_eps
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use boozer_field, only: boozer_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_of_fieldlines

    use precession, only: bounce_integrals_t
    use precession, only: make_fieldlines_for_precession
    use precession, only: assemble_bounce_integrals
    use precession, only: compute_flux_modes
    use precession, only: use_dkes_electric_drift
    use precession, only: get_boundary_bounce_length
    use precession, only: get_maximum_valid_precession
    use grid_mod, only: fieldline_for_precession_t
    use grid_mod, only: compute_bounce_integrals
    use saved_bounce_integrals, only: bounce_integrals_header_t
    use saved_bounce_integrals, only: write_bounce_integrals_header
    use saved_bounce_integrals, only: append_bounce_integrals, read_bounce_integrals
    use field_instance, only: initialize_field_instance
    use fieldline_labels, only: modes_t, calc_rms_mode_number

    use surface_average_mod, only: surface_average_t, calc_surface_averages
    use coefficients, only: calc_nu_star_crit
    use netcdf_mod, only: netcdf_t
    use fieldline_labels, only: fieldline_modes_t, fourier_transform_over_label
    use git_version, only: git_hash

    use read_file, only: read_namelist
    use read_file, only: field_file, &
                         field_type, &
                         n_rho, n_theta, n_zeta, &
                         M_pol, &
                         N_tor, &
                         s_tor, &
                         sign_sqrtg, &
                         max_n_fieldlines
    use read_precession, only: read_precession_namelist
    use read_precession, only: Omega_hat, &
                               nu_star, &
                               l_c_input => l_c, &
                               n_eta, &
                               electric_drift

    implicit none

    character(len=*), parameter :: input_file = "rabe.in"
    character(len=*), parameter :: output_file = "electric_rabe.nc"
    character(len=*), parameter :: bounce_integrals_file = "saved_bounce_integrals.dat"

    type(boozer_field_t) :: field

    integer :: n_stor
    integer :: this

    real(dp) :: R ![m]
    type(surface_average_t) :: average
    real(dp) :: iota
    real(dp) :: B_theta_covariant, B_phi_covariant
    real(dp) :: nfp

    integer :: n_fieldlines
    type(flock_of_fieldlines_t) :: flock

    real(dp) :: covariant_factor
    real(dp), dimension(:), allocatable :: nu_star_crit
    real(dp), dimension(:), allocatable :: l_c_crit
    ! collisionality grid as collision length l_c [m] and as nu_star with the
    ! code's R; the input gives one of them, the other follows from
    ! l_c = pi*R/(2*nu_star)
    real(dp), dimension(:), allocatable :: l_c_grid, nu_star_grid
    integer, dimension(:), allocatable :: err_flag

    type(fieldline_for_precession_t), dimension(:), allocatable :: fieldlines_precession
    type(bounce_integrals_t) :: bounce_integrals

    real(dp), dimension(:), allocatable :: flux_mode
    type(modes_t) :: g_off_modes

    integer :: n_modes
    integer :: idx
    logical :: file_exists
    type(bounce_integrals_header_t) :: bounce_integrals_header
    logical, dimension(:), allocatable :: is_complete
    integer :: n_Omega_hat
    integer :: idx_Omega
    integer :: n_nu
    integer :: id_nu
    real(dp) :: l_c
    real(dp), dimension(:, :), allocatable :: lambda_off
    real(dp), dimension(:, :, :), allocatable :: lambda_off_modes
    ! rms mode number of g_off over the field line label, the poloidal
    ! precession limit per collisionality and the lowest of them
    real(dp), dimension(:), allocatable :: g_off_m_rms
    real(dp), dimension(:), allocatable :: Omega_hat_limits
    real(dp), dimension(:), allocatable :: Omega_hat_lowest_limit

    type(netcdf_t) :: nc_output
    character(len=*), parameter :: dim_name = "surface"
    character(len=1024) :: description
    character(len=1) :: radial_label

    call read_namelist(input_file)
    call read_precession_namelist(input_file)

    n_stor = size(s_tor)
    if (n_stor > 1) then
        error stop "Error: only one surface allowed!"
    end if
    allocate (nu_star_crit(n_stor))
    allocate (l_c_crit(n_stor))
    allocate (err_flag(n_stor))
    allocate (Omega_hat_lowest_limit(n_stor))
    n_Omega_hat = size(Omega_hat)
    n_nu = max(size(nu_star), size(l_c_input))
    allocate (l_c_grid(n_nu), nu_star_grid(n_nu))
    allocate (lambda_off(n_nu, n_Omega_hat))
    allocate (g_off_m_rms(n_nu), Omega_hat_limits(n_nu))

    select case (trim(field_type))
    case ('vmec_nc')
        call field%init_from_vmec(field_file, grid_refinement=6)
    case ('booz_xform')
        call field%init_from_boozmn(field_file, n_rho, n_theta, n_zeta)
    case ('chartmap')
        call field%init_from_chartmap(field_file)
    case default
        error stop "electric_rabe: unknown field_type"
    end select
    do this = 1, n_stor
        call field%fix_to_surface(s_tor(this))
        call field%get_iota(s_tor(this), iota)
        call field%get_covariant_components(B_theta_covariant, &
                                            B_phi_covariant)
        nfp = field%nfp

        call make_flock_of_fieldlines(flock, max_n_fieldlines, iota, &
                                      field, M_pol, N_tor, nfp, &
                                      err_flag(this))
        iota = flock%iota
        n_fieldlines = size(flock%fieldlines)

        covariant_factor = (B_phi_covariant + B_theta_covariant*iota)
        call calc_surface_averages(flock, average)
        R = field%R
        if (size(l_c_input) > 0) then
            l_c_grid = l_c_input
            nu_star_grid = pi*R/(2.0_dp*l_c_grid)
        else
            nu_star_grid = nu_star
            l_c_grid = pi*R/(2.0_dp*nu_star_grid)
        end if

! --------------------------------------------------------------------------
        call make_fieldlines_for_precession(field, flock, fieldlines_precession, n_eta)
        call initialize_field_instance(field)

        bounce_integrals_header = &
            bounce_integrals_header_t(git_hash=git_hash, &
                                      field_file=trim(field_file), &
                                      s_tor=s_tor(this), &
                                      M_pol=M_pol, &
                                      N_tor=N_tor, &
                                      n_grid=n_eta, &
                                      max_n_fieldlines=max_n_fieldlines, &
                                      n_fieldlines=n_fieldlines, &
                                      eta_c=fieldlines_precession(1)%grid%eta(n_eta), &
                                      eta_t=fieldlines_precession(1)%grid%eta(1))
        allocate (is_complete(n_fieldlines))
        inquire (file=bounce_integrals_file, exist=file_exists)
        if (file_exists) then
            call read_bounce_integrals(bounce_integrals_file, &
                                       bounce_integrals_header, &
                                       fieldlines_precession, is_complete)
            print *, "Read ", count(is_complete), " of ", n_fieldlines, &
                " fieldlines from ", bounce_integrals_file
        else
            call write_bounce_integrals_header(bounce_integrals_file, &
                                               bounce_integrals_header)
            is_complete = .false.
        end if

        do idx = 1, n_fieldlines
            if (.not. is_complete(idx)) then
                call compute_bounce_integrals(field, &
                                              fieldlines_precession(idx), &
                                              s_tor(this), &
                                              fieldlines_precession(idx)%grid)
                call append_bounce_integrals(bounce_integrals_file, idx, &
                                             fieldlines_precession(idx)%grid)
            end if
        end do
        call assemble_bounce_integrals(fieldlines_precession, bounce_integrals)
        if (trim(electric_drift) == 'dkes') then
            call use_dkes_electric_drift(bounce_integrals, s_tor(this))
        end if
        print *, "E x B drift: ", trim(electric_drift)
        n_modes = bounce_integrals%n_modes

        allocate (lambda_off_modes(n_modes, n_nu, n_Omega_hat))

        ! g_off without precession does not depend on Omega_hat; the limit in
        ! the label xi_0 is converted to the poloidal convention of the input
        ! (xi_0 drifts by Omega_hat nfp/(M iota - N), see compute_flux_modes)
        do id_nu = 1, n_nu
            call get_g_modes_from_fieldlines(flock, l_c_grid(id_nu), g_off_modes, &
                                             covariant_factor)
            g_off_m_rms(id_nu) = calc_rms_mode_number(g_off_modes)
            Omega_hat_limits(id_nu) = get_maximum_valid_precession( &
                                      get_boundary_bounce_length(bounce_integrals), &
                                      g_off_m_rms(id_nu)) &
                                      *abs(bounce_integrals%M_pol*bounce_integrals%iota &
                                           - bounce_integrals%N_tor)/bounce_integrals%nfp
        end do
        Omega_hat_lowest_limit(this) = minval(Omega_hat_limits)
        print *, "lowest Omega_hat limit: ", Omega_hat_lowest_limit(this)

        do idx_Omega = 1, n_Omega_hat
            print *, "idx_Omega: ", idx_Omega, " / ", n_Omega_hat
            do id_nu = 1, n_nu
                l_c = l_c_grid(id_nu)
                call compute_flux_modes(bounce_integrals, l_c, Omega_hat(idx_Omega), &
                                        flux_mode)
                call get_g_modes_from_fieldlines(flock, l_c, g_off_modes, &
                                                 covariant_factor)
        lambda_off_modes(:, id_nu, idx_Omega) = pi*g_off_modes%sin_coeffs*flux_mode/average%normalization
               lambda_off(id_nu, idx_Omega) = sum(lambda_off_modes(:, id_nu, idx_Omega))

            end do
        end do
        lambda_off = lambda_off*0.75_dp/covariant_factor
        lambda_off_modes = lambda_off_modes*0.75_dp/covariant_factor
        ! per effective radius r if <|grad s|> is known, otherwise per s
        if (ieee_is_nan(average%nabla_s) .or. abs(average%nabla_s) < machine_eps) then
            radial_label = "s"
        else
            radial_label = "r"
            lambda_off = lambda_off/average%nabla_s
            lambda_off_modes = lambda_off_modes/average%nabla_s
        end if
        print *, "lambda_off given with respect to ", radial_label

        deallocate (fieldlines_precession)
        deallocate (is_complete)

! --------------------------------------------------------------------------

        nu_star_crit(this) = calc_nu_star_crit(flock, R)
        l_c_crit(this) = pi*R/(2.0_dp*nu_star_crit(this))

        print *, "s_tor: ", s_tor(this)

    end do

    call nc_output%create(output_file)
    call nc_output%add_global_attribute("title", &
                                        "asymptotic bootstrap coefficient lambda_bB")
    call nc_output%add_global_attribute("definition", &
                                        "lambda^{off}_bB = ?")
    call nc_output%add_global_attribute("git_hash", git_hash)
    call nc_output%def_dim(dim_name, n_stor)
    call nc_output%add_real_1d("s_tor", dim_name)
    call nc_output%add_attr("s_tor", "long_name", &
                            "normalized toroidal flux label")
    call nc_output%add_attr("s_tor", "unit", &
                            "[1]")

    call nc_output%def_dim("collisionality", n_nu)
    call nc_output%add_real_1d("nu_star", "collisionality")
    call nc_output%add_attr("nu_star", "long_name", &
                            "collisionality pi*R/(2*l_c) with the major radius R of this file")
    call nc_output%add_real_1d("l_c", "collisionality")
    call nc_output%add_attr("l_c", "long_name", &
                            "collision length l_c = pi*R/(2*nu_star); independent of the choice of R")
    call nc_output%add_attr("l_c", "unit", "[m]")
    call nc_output%add_real_1d("g_off_m_rms", "collisionality")
    call nc_output%add_attr("g_off_m_rms", "long_name", &
                            "rms mode number of g_off over the field line label")
    call nc_output%add_attr("g_off_m_rms", "unit", "[1]")
    call nc_output%add_real_1d("Omega_hat_limits", "collisionality")
    call nc_output%add_attr("Omega_hat_limits", "long_name", &
                            "upper limit of the poloidal Omega_hat: drift of the "// &
                            "field line label during a bounce below 1/g_off_m_rms")
    call nc_output%add_attr("Omega_hat_limits", "unit", "[1/m]")
    call nc_output%add_real_1d("Omega_hat_lowest_limit", dim_name)
    call nc_output%add_attr("Omega_hat_lowest_limit", "long_name", &
                            "lowest Omega_hat_limits, valid for all collisionalities")
    call nc_output%add_attr("Omega_hat_lowest_limit", "unit", "[1/m]")
    call nc_output%def_dim("precession_frequency", n_Omega_hat)
    call nc_output%add_real_1d("Omega_hat", "precession_frequency")

    call nc_output%add_real_2d("lambda_off", "collisionality", "precession_frequency")
    call nc_output%add_attr("lambda_off", "long_name", &
                            "geometrical factor")
    call nc_output%add_attr("lambda_off", "unit", &
                            "[1]")
    call nc_output%add_attr("lambda_off", "radial_label", radial_label)
    call nc_output%add_attr("lambda_off", "electric_drift", trim(electric_drift))

    call nc_output%def_dim("label_mode", n_modes)
    call nc_output%add_real_1d("mode_numbers", "label_mode")
    call nc_output%add_attr("mode_numbers", "long_name", &
                            "Fourier mode numbers")
    call nc_output%add_real_3d("lambda_off_modes", &
                               "label_mode", "collisionality", &
                               "precession_frequency")
    call nc_output%add_attr("lambda_off_modes", "long_name", &
                            "per-mode contribution to lambda_off")
    call nc_output%add_attr("lambda_off_modes", "unit", &
                            "[1]")

    call nc_output%add_real_1d("nu_star_crit", dim_name)
    call nc_output%add_attr("nu_star_crit", "long_name", &
                          "lower collisionality limit for validity of asymptotic model")
    call nc_output%add_attr("nu_star_crit", "unit", &
                            "[1]")
    call nc_output%add_real_1d("l_c_crit", dim_name)
    call nc_output%add_attr("l_c_crit", "long_name", &
                            "upper collision length limit for validity of asymptotic model")
    call nc_output%add_attr("l_c_crit", "unit", "[m]")
    call nc_output%add_int_1d("err_flag", dim_name)
    call nc_output%add_attr("err_flag", "long_name", &
                            "1 if violation of omnigeneity is too strong, 0 otherwise")
    call nc_output%add_real("R")
    call nc_output%add_attr("R", "long_name", &
                            "major radius")
    call nc_output%add_attr("R", "unit", "[m]")
    write (description, "(A,A,A)") "defines the reference length scale to convert ", &
        "to dimensionless quantities i.e. defines coeffients in respect to ", &
        "nu_star = pi*R/(2*mean_free_path) = pi*R*deflection_frequency/particle_speed"
    call nc_output%add_attr("R", "definition", description)
    call nc_output%write_real_2d("lambda_off", lambda_off)
    call nc_output%write_real_1d("mode_numbers", g_off_modes%mode_numbers)
    call nc_output%write_real_3d("lambda_off_modes", lambda_off_modes)
    call nc_output%write_real_1d("Omega_hat", Omega_hat)
    call nc_output%write_real_1d("nu_star", nu_star_grid)
    call nc_output%write_real_1d("l_c", l_c_grid)
    call nc_output%write_real_1d("g_off_m_rms", g_off_m_rms)
    call nc_output%write_real_1d("Omega_hat_limits", Omega_hat_limits)
    call nc_output%write_real_1d("Omega_hat_lowest_limit", Omega_hat_lowest_limit)
    call nc_output%write_real_1d("nu_star_crit", nu_star_crit)
    call nc_output%write_real_1d("l_c_crit", l_c_crit)
    call nc_output%write_int_1d("err_flag", err_flag)
    call nc_output%write_real_1d("s_tor", s_tor)
    call nc_output%write_real("R", field%R)
    call nc_output%close()

    if (allocated(lambda_off_modes)) deallocate (lambda_off_modes)
    if (allocated(nu_star_crit)) deallocate (nu_star_crit)
    if (allocated(l_c_crit)) deallocate (l_c_crit)
    if (allocated(l_c_grid)) deallocate (l_c_grid)
    if (allocated(g_off_m_rms)) deallocate (g_off_m_rms)
    if (allocated(Omega_hat_limits)) deallocate (Omega_hat_limits)
    if (allocated(Omega_hat_lowest_limit)) deallocate (Omega_hat_lowest_limit)
    if (allocated(nu_star_grid)) deallocate (nu_star_grid)
    if (allocated(err_flag)) deallocate (err_flag)

contains

    subroutine get_g_modes_from_fieldlines(flock, l_c, g_off_modes, covariant_factor)
        use fit_functions, only: S_A, S_B
        use fieldline_labels, only: allocate_modes
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), intent(in) :: l_c
        real(dp), intent(in) :: covariant_factor
        type(modes_t), intent(out) :: g_off_modes

        type(fieldline_modes_t) :: modes
        type(surface_average_t) :: averages

        real(dp) :: prefactor_A

        integer :: max_mode

        call fourier_transform_over_label(flock, modes)
        call calc_surface_averages(flock, averages)

        max_mode = size(modes%delta_eta%cos_coeffs, dim=1)
        if (.not. allocated(g_off_modes%sin_coeffs)) then
            call allocate_modes(g_off_modes, max_mode)
        end if

        ! flock%I_ref already contains the covariant factor (857dc93)
        prefactor_A = 2.0_dp*sqrt(flock%eta_b*flock%I_ref/l_c)
        g_off_modes%sin_coeffs = prefactor_A*modes%delta_aspect_ratio%cos_coeffs* &
                                 S_A(flock%iota_p*modes%delta_aspect_ratio%mode_numbers)
        g_off_modes%sin_coeffs = g_off_modes%sin_coeffs + &
                                 modes%delta_eta%cos_coeffs* &
                                 S_B(flock%iota_p*modes%delta_eta%mode_numbers)
        g_off_modes%sin_coeffs = g_off_modes%sin_coeffs* &
                                 averages%B_squared/averages%lambda_b* &
                                 l_c*0.5_dp
    end subroutine get_g_modes_from_fieldlines
end program electric_rabe
