program test_plot
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace
    implicit none

    type(myplot) :: plt

    integer, parameter :: n_points = 100
    real(dp), dimension(n_points) :: x, f

    call linspace(0.0_dp, 2.0_dp*pi, n_points, x)
    f = sin(x)

    call plt%initialize()
    call plt%add_plot(x, f, label="sin(x)", linestyle="--")
    call plt%save("sin")

end program test_plot
