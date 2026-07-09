module bc_file
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    implicit none

contains

    subroutine write_field_B_mod_to_bc(field, s_tors, m_max, n_max, filename)
        type(boozer_field_t), intent(inout) :: field
        real(dp), dimension(:), intent(in) :: s_tors
        integer, intent(in) :: m_max, n_max
        character(len=*), intent(in) :: filename

        integer :: is, i_m, i_n, idx
        integer :: number_of_modes
        integer, parameter :: grid_multiplier = 9
        integer :: n_surf
        real(dp) :: s_tor
        integer :: nu, nv
        real(dp) :: iota, B_theta_covariant, B_phi_covariant
        real(dp) :: Jpol_over_nper, Itor
        integer :: iunit
        integer, dimension(:), allocatable :: m, n
        real(dp), dimension(:), allocatable :: bmnc
        real(dp), dimension(:), allocatable :: rmnc, zmns, vmns
        real(dp) :: dummy_minor_radius
        real(dp), parameter :: magnetic_to_current = -0.5_dp*1e7

        number_of_modes = (m_max + 1)*(2*n_max + 1) - n_max
        allocate (m(number_of_modes))
        allocate (n(number_of_modes))
        allocate (bmnc(number_of_modes))
        n_surf = size(s_tors)
        nu = grid_multiplier*m_max
        nv = grid_multiplier*n_max + 1
        idx = 0
        do i_m = 0, m_max
            do i_n = -n_max, n_max
                if (i_m == 0 .and. i_n < 0) cycle
                idx = idx + 1
                m(idx) = i_m
                n(idx) = i_n
            end do
        end do
        n = n*nint(field%nfp)

        allocate (rmnc(number_of_modes))
        allocate (zmns(number_of_modes))
        allocate (vmns(number_of_modes))

        rmnc = 0.0_dp
        zmns = 0.0_dp
        vmns = 0.0_dp
        rmnc(1) = field%R
        dummy_minor_radius = 0.5_dp*field%R
        rmnc(2*(n_max + 1)) = dummy_minor_radius
        zmns(2*(n_max + 1) + 1) = dummy_minor_radius

        open (newunit=iunit, file=trim(filename), status='replace', action='write')

        ! --- Global header ---
        write (iunit, '(A)') 'CC Boozer-coordinate data file'
        write (iunit, '(A)') 'CC Version:'
        write (iunit, '(A)') 'CC Author:'
        write (iunit, '(A)') 'CC shot:    0'
        write (iunit, '(A)') ' m0b   n0b  nsurf  nper    flux [Tm^2]'// &
            '        a [m]          R [m]'
        write (iunit, '(2I6, 2I6, 3ES16.6)') &
            m_max, n_max, n_surf, int(field%nfp), field%psi_tor_edge, 0.0_dp, field%R

        ! --- Per-surface blocks ---
        do is = 1, n_surf
            print *, "Writing surface ", is, " of ", n_surf
            s_tor = s_tors(is)
            call field%get_iota(s_tor, iota)
            call field%fix_to_surface(s_tor)
            call field%get_covariant_components(B_theta_covariant, B_phi_covariant)
            Jpol_over_nper = B_phi_covariant/field%nfp*magnetic_to_current
            Itor = B_theta_covariant*magnetic_to_current

            ! Surface header
            write (iunit, '(A)') '        s               iota'// &
                '           Jpol/nper          Itor            pprime'// &
                '         sqrt g(0,0)'
            write (iunit, '(A)') '                                          [A]'// &
                '           [A]             [Pa]         (dV/ds)/nper'
            write (iunit, '(6ES16.8)') s_tor, iota, &
                Jpol_over_nper, Itor, 0.0_dp, 0.0_dp

            ! Mode table header
            write (iunit, '(A)') '    m    n      rmnc [m]'// &
                '         zmns [m]         vmns [ ]         bmnc [T]'

            call compute_2D_fourier(field, m, n, nu=nu, nv=nv, bmnc=bmnc, &
                                    nfp=field%nfp)
            do idx = 1, size(m)
                write (iunit, '(2I6, 4ES16.8)') &
               m(idx), n(idx)/int(field%nfp), rmnc(idx), zmns(idx), vmns(idx), bmnc(idx)
            end do

        end do

        close (iunit)

    end subroutine write_field_B_mod_to_bc

    subroutine compute_2D_fourier(field, m, n, nu, nv, bmnc, nfp)
        use constants, only: pi
        use field_base, only: field_t
        class(field_t), intent(in) :: field

        integer, intent(in) :: nu, nv
        integer, dimension(:), intent(in) :: m, n
        real(dp), dimension(:), intent(out) :: bmnc
        real(dp), intent(in) :: nfp

        integer :: i, j, k, idx
        real(dp) :: du, dv
        real(dp) :: prefac
        real(dp) :: u(nu*nv), v(nu*nv), fval(nu*nv)

        du = 2.0_dp*pi/real(nu, kind=dp)
        dv = 2.0_dp*pi/real(nv, kind=dp)/nfp
        prefac = 2.0_dp/(2.0_dp*pi)**2*du*dv*nfp

        idx = 0
        do j = 1, nv
            do i = 1, nu
                idx = idx + 1
                u(idx) = real(i - 1, kind=dp)*du
                v(idx) = real(j - 1, kind=dp)*dv
                call field%compute_B_mod(u(idx), v(idx), fval(idx))
            end do
        end do

        do k = 1, size(m)
           bmnc(k) = prefac*sum(cos(real(m(k), kind=dp)*u - real(n(k), kind=dp)*v)*fval)
            if (m(k) == 0 .and. n(k) == 0) bmnc(k) = bmnc(k)/2.0_dp
        end do

    end subroutine compute_2D_fourier

    subroutine delete_bc_file(filename)
        character(len=*), intent(in) :: filename

        integer :: iunit
        logical :: exists
        inquire (file=filename, exist=exists)
        if (exists) then
            open (newunit=iunit, file=filename, status='old')
            close (iunit, status='delete')
            inquire (file=filename, exist=exists)
            if (exists) then
                error stop "Failed to delete bc file: "//trim(filename)
            end if
        end if
    end subroutine delete_bc_file

end module bc_file
