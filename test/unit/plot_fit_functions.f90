program plot_fit_functions
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace
    use fit_functions, only: S_A, S_B
    implicit none

    type(myplot) :: plt

    integer, parameter :: n_points = 100
    real(dp), dimension(n_points) :: angle

    call linspace(-2.0_dp*pi, 2.0_dp*pi, n_points, angle)

    call plt%initialize(xlabel="angle [$\pi$]", &
                        ylabel="value [$\pi$]", &
                        legend=.true.)
    call plt%add_plot(angle/pi, &
                      S_A(angle)/pi, &
                      label="$S_A$", &
                      linestyle="r-")
    call plt%add_plot(angle/pi, &
                      S_B(angle)/pi, &
                      label="$S_B$", &
                      linestyle="b--")
    call plt%save("fit_functions")

end program plot_fit_functions
