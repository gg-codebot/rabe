program test_fft
    use constants, only: dp
    use utils, only: not_same
    use fourier_field, only: ifft_modes_to_B
    use horner_fourier_field, only: ift_modes_to_B

    implicit none

    real(dp), parameter :: reltol = 1e-14
    integer, parameter :: n_m = 6, n_n = 3, n_total = n_m*n_n
    integer :: m(n_total), n(n_total)
    real(dp) :: B_modes(n_total)
    integer :: m_idx, n_idx, k_mode

    real(dp), dimension(:, :), allocatable :: ifft_B, ift_B
    real(dp), dimension(:), allocatable :: theta, phi
    real(dp) :: max_rel_err_B
    integer, parameter :: n_grid = 257
    integer, parameter :: n_fft = n_grid - 1

    logical :: test_failed

    test_failed = .false.

    k_mode = 0
    do m_idx = 1, n_m
        do n_idx = 1, n_n
            k_mode = k_mode + 1
            m(k_mode) = m_idx - 1
            n(k_mode) = n_idx - 2         ! n = -1, 0, 1
            if (m(k_mode) == 0 .and. n(k_mode) == 0) then
                B_modes(k_mode) = 1.0_dp
            else
                B_modes(k_mode) = 0.01_dp
            end if
        end do
    end do

    allocate (ifft_B(n_fft, n_fft))
    call ifft_modes_to_B(m, n, B_modes, n_fft, ifft_B)
    allocate (ift_B(n_grid, n_grid))
    allocate (theta(n_grid), phi(n_grid))
    call ift_modes_to_B(m, n, B_modes, ift_B, theta, phi, nfp=1, n_grid=n_grid)

    if (not_same(ifft_B, ift_B(1:n_fft, 1:n_fft), &
                 reltol_in=reltol, abstol_in=0.0_dp)) then
        print *, "Test failed: ifft_modes_to_B and ift_modes_to_B do not match."
        print *, "Max relative error:", &
            maxval(abs(ifft_B - ift_B(1:n_fft, 1:n_fft))/abs(ifft_B))
        print *, "bigger than tolerance:", reltol
        test_failed = .true.
    end if

    if (test_failed) error stop "test_fft failed!"

end program test_fft
