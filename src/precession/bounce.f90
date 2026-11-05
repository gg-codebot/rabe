module bounce
    use constants, only: dp, pi
    use params, only: ntimstep, nplagr, nder, npl_half, dtaumin, relerr, ntau

    implicit none

contains

    subroutine trace_orbit_till_bounce(z_start, z_end, orbit_buffer)
        use plag_coeff_sub, only: plag_coeff
        use orbit, only: orbit_timestep_axis
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        integer, parameter :: ndim = 8
        real(dp), dimension(ndim), intent(in) :: z_start
        real(dp), dimension(ndim), intent(out) :: z_end
        real(dp), dimension(:, :), intent(inout), optional :: orbit_buffer

        integer :: ierr
        real(dp), dimension(ndim) :: z
        integer :: it, ktau
        integer(8) :: kt

        integer, dimension(:), allocatable :: ipoi
        real(dp), dimension(:), allocatable :: xp
        real(dp), dimension(:, :), allocatable :: coef, orb_sten
        real(dp), parameter :: zerolam = 0.0_dp
        real(dp) :: alam_prev
        integer :: itip
        logical :: did_bounce

        did_bounce = .false.
        z_end = 0d0

        z = z_start

        if (present(orbit_buffer)) orbit_buffer = ieee_value(1.0_dp, ieee_quiet_nan)

        !$omp critical
        if (.not. allocated(ipoi)) then
            allocate (ipoi(nplagr))
            allocate (coef(0:nder, nplagr))
            allocate (orb_sten(ndim, nplagr))
            allocate (xp(nplagr))
        end if
        !$omp end critical
        do it = 1, nplagr
            ipoi(it) = it
        end do

        kt = 0

        itip = npl_half + 1
        alam_prev = z(5)

        ! End initialize tip detector
        !--------------------------------
        !
        TIMELOOP: do ktau = 1, ntau*(ntimstep - 1)
            call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr)
            if (ierr .ne. 0) then
                error stop "Orbit lost during integration!"
            end if
            kt = kt + 1
            if (present(orbit_buffer)) then
                if (kt <= size(orbit_buffer, 1)) then
                    orbit_buffer(kt, :) = z
                end if
            end if
            if (kt .le. nplagr) then          !<=first nplagr points to initialize stencil
                orb_sten(:, kt) = z
            else                          !<=normal case, shift stencil
                orb_sten(:, ipoi(1)) = z
                ipoi = cshift(ipoi, 1)
            end if

            ! Tip detection and interpolation
            if (alam_prev .lt. 0.d0 .and. z(5) .gt. 0.d0) itip = 0   !<=tip has been passed
            itip = itip + 1
            alam_prev = z(5)
            if (kt .gt. nplagr) then          !<=use only initialized stencil
                if (itip .eq. npl_half) then   !<=stencil around tip is complete, interpolate
                    xp = orb_sten(5, ipoi)

                    call plag_coeff(nplagr, nder, zerolam, xp, coef)

                    z = matmul(orb_sten(:, ipoi), coef(0, :))
                    z(2) = modulo(z(2), 2.0_dp*pi)
                    z(3) = modulo(z(3), 2.0_dp*pi)
                    if (present(orbit_buffer)) then
                        orbit_buffer(kt, :) = z
                    end if

                    did_bounce = .true.
                    exit TIMELOOP
                end if
            end if
            ! End tip detection and interpolation
        end do TIMELOOP

        !$omp critical
        z_end = z
        !$omp end critical
        if (.not. did_bounce) then
            print *, "Warning: trace_orbit_till_bounce reached end of time loop without detecting bounce."
        end if
    end subroutine trace_orbit_till_bounce

end module bounce
