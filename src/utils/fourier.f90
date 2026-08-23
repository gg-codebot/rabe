module fourier
    use constants, only: dp, pi
    implicit none
    private

    public :: real_ft
    public :: check_is_equidistant, check_has_correct_endpoints

contains

    subroutine real_ft(x, f, f_cos, f_sin)
        real(dp), dimension(:), intent(in) :: x, f
        real(dp), dimension(size(f)/2 + 1), intent(out) :: f_cos, f_sin

        integer :: N, k, j
        real(dp) :: dx
        complex(dp) :: phase_step, phase, sum
        complex(dp), dimension(size(f)/2 + 1) :: f_exp
        complex(dp), parameter :: Icmplx = (0.0_dp, 1.0_dp)

        call check_is_equidistant(x)
        call check_has_correct_endpoints(x)

        N = size(f)
        dx = x(2) - x(1)

        do k = 0, N/2
            phase_step = exp(-Icmplx*dx*real(k, kind=dp))
            phase = (1.0_dp, 0.0_dp)*exp(-Icmplx*x(1)*real(k, kind=dp))
            sum = (0.0_dp, 0.0_dp)
            do j = 0, N - 1
                sum = sum + f(j + 1)*phase
                phase = phase*phase_step
            end do
            f_exp(k + 1) = sum
        end do

        f_cos(1) = real(f_exp(1))
        f_sin(1) = 0.0_dp
        do k = 1, N/2
            f_cos(k + 1) = 2.0_dp*real(f_exp(k + 1))
            f_sin(k + 1) = -2.0_dp*aimag(f_exp(k + 1))
        end do

        if (mod(N, 2) == 0) then
            f_cos(N/2 + 1) = f_cos(N/2 + 1)/2.0_dp
            f_sin(N/2 + 1) = 0.0_dp
        end if

        f_cos = f_cos/real(N, kind=dp)
        f_sin = f_sin/real(N, kind=dp)

    end subroutine real_ft

    subroutine check_is_equidistant(x)
        real(dp), dimension(:), intent(in) :: x
        logical :: is_equidistant

        real(dp), parameter :: retol = 1e-13, abstol = 1e-14
        real(dp), dimension(size(x) - 1) :: dx
        real(dp) :: tol
        integer :: N

        N = size(x)

        dx = x(2:N) - x(1:N - 1)
        is_equidistant = all(abs(dx - dx(1)) < retol*dx(1) + abstol)

        if (.not. is_equidistant) then
            print *, "Input x has to be equidistant!"
            print *, "violation by ", maxval(abs(dx - dx(1)))
            error stop
        end if
    end subroutine check_is_equidistant

    subroutine check_has_correct_endpoints(x)
        real(dp), dimension(:), intent(in) :: x
        logical :: has_correct_range

        real(dp), parameter :: tol = 1e-13
        real(dp) :: correct_range, range
        integer :: N

        N = size(x)
        range = x(N) - x(1)
        correct_range = 2.0_dp*pi*(1 - 1/real(N, kind=dp))
        has_correct_range = abs(correct_range - range) < tol*correct_range
        if (.not. has_correct_range) then
            print *, "Input x has wrong endpoints!"
            print *, "actual: ", x(1), x(N)
            print *, "required: ", x(1), x(1) + correct_range
            error stop
        end if
    end subroutine check_has_correct_endpoints

end module fourier
