module read_precession
    use constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none

    real(dp), dimension(:), allocatable, public, protected :: Omega_hat
    ! collisionality grid: either nu_star or l_c = pi*R/(2*nu_star) [m] is given;
    ! the other one is allocated with size 0
    real(dp), dimension(:), allocatable, public, protected :: nu_star
    real(dp), dimension(:), allocatable, public, protected :: l_c
    integer, public, protected :: n_eta ! eta grid points of the bounce integrals
    ! E x B drift with 1/B^2 ('exact') or with 1/<B^2> as in DKES/MONKES ('dkes')
    character(len=8), public, protected :: electric_drift

    namelist /precession_config/ &
        Omega_hat, &
        nu_star, &
        l_c, &
        n_eta, &
        electric_drift

contains

    subroutine read_precession_namelist(filename)
        implicit none
        character(len=*), intent(in) :: filename
        integer :: ios, unit
        logical :: file_exists

        real(dp), dimension(:), allocatable :: Omega_hat_temp
        integer, parameter :: n_Omega_hat_max = 1000
        integer :: n_Omega_hat
        real(dp), dimension(:), allocatable :: nu_star_temp
        integer, parameter :: n_nu_star_max = 1000
        integer :: n_nu_star, n_l_c
        real(dp) :: nan_value

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

        ! Default values
        nan_value = ieee_value(nan_value, ieee_quiet_nan)
        if (allocated(Omega_hat)) deallocate (Omega_hat)
        allocate (Omega_hat(n_Omega_hat_max))
        Omega_hat = nan_value
        if (allocated(nu_star)) deallocate (nu_star)
        allocate (nu_star(n_nu_star_max))
        nu_star = nan_value
        if (allocated(l_c)) deallocate (l_c)
        allocate (l_c(n_nu_star_max))
        l_c = nan_value
        n_eta = 200
        electric_drift = 'exact'

        read (unit, nml=precession_config, iostat=ios)
        if (ios /= 0) then
            print *, "Error in read_precession_namelist:"
            print *, "iostat = ", ios
            error stop
        end if

        n_nu_star = count(.not. ieee_is_nan(nu_star))
        n_l_c = count(.not. ieee_is_nan(l_c))
        if (n_nu_star == 0 .and. n_l_c == 0) then
            print *, "Error in read_reprecession_namelist: ", &
                "neither nu_star nor l_c values provided in namelist"
            error stop
        end if
        if (n_nu_star > 0 .and. n_l_c > 0) then
            print *, "Error in read_reprecession_namelist: ", &
                "give either nu_star or l_c, not both"
            error stop
        end if
        if (max(n_nu_star, n_l_c) >= n_nu_star_max) then
            print *, "Error in read_reprecession_namelist: ", &
                "too many nu_star or l_c values provided in namelist"
            error stop
        end if
        if (any(l_c(1:n_l_c) <= 0.0_dp)) then
            print *, "Error in read_reprecession_namelist: l_c must be positive"
            error stop
        end if
        if (allocated(nu_star_temp)) deallocate (nu_star_temp)
        allocate (nu_star_temp(n_nu_star))
        nu_star_temp = pack(nu_star,.not. ieee_is_nan(nu_star))
        deallocate (nu_star)
        allocate (nu_star(n_nu_star))
        nu_star = nu_star_temp
        deallocate (nu_star_temp)
        allocate (nu_star_temp(n_l_c))
        nu_star_temp = pack(l_c,.not. ieee_is_nan(l_c))
        deallocate (l_c)
        allocate (l_c(n_l_c))
        l_c = nu_star_temp
        deallocate (nu_star_temp)

        n_Omega_hat = count(.not. ieee_is_nan(Omega_hat))
        if (n_Omega_hat >= n_Omega_hat_max) then
            print *, "Error in read_reprecession_namelist: ", &
                "too many Omega_hat values provided in namelist"
            error stop
        end if
        if (allocated(Omega_hat_temp)) deallocate (Omega_hat_temp)
        allocate (Omega_hat_temp(n_Omega_hat))
        Omega_hat_temp = pack(Omega_hat,.not. ieee_is_nan(Omega_hat))
        deallocate (Omega_hat)
        allocate (Omega_hat(n_Omega_hat))
        Omega_hat = Omega_hat_temp
        deallocate (Omega_hat_temp)

        if (trim(electric_drift) /= 'exact' .and. trim(electric_drift) /= 'dkes') then
            print *, "Error in read_precession_namelist: electric_drift must be ", &
                "'exact' or 'dkes'"
            error stop
        end if

        close (unit)
    end subroutine read_precession_namelist

end module read_precession
