<!-- When adding new input files or changing the myplot API, update this README accordingly. -->

# Plot Tests

## CMake pattern

```cmake
add_executable(my_plot.x my_plot.f90)
target_link_libraries(my_plot.x myplot_lib rabe_lib helpers_lib)
add_test(MyPlot my_plot.x)
set_tests_properties(MyPlot PROPERTIES LABELS "plot")
```

Link `myplot_lib` always. Add `rabe_lib` and `helpers_lib` when using field types.

## Input files

The CMake build copies files from `test/integration/input/` and `test/integration/vmec/input/` into `build/test/plot/input/`. Reference them as `"input/<filename>"` in Fortran.

Available input files:
- **NEO `.bc` files**: `quasi_helical.bc`, `helical_anti.bc`, `helical_anti_small_aspect.bc`, `poloidal_anti_minuspert.bc`, `single_mode_m_2_n_minus4.bc`
- **VMEC `.nc` files**: `wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc`, `wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc`
- **VMEC-derived `.bc` files**: `landreman_paul_qh.bc`, `landreman_paul_qa.bc`

## myplot API (`use myplot_module, only: myplot`)

```fortran
type(myplot) :: plt

! Initialize figure
call plt%initialize(xlabel, ylabel, legend, figsize, title, fontsize)
!   xlabel, ylabel: character, optional (LaTeX via raw strings)
!   legend: logical, optional
!   figsize: integer(2), optional, default [10,8]
!   title: character, optional
!   fontsize: integer, optional, default 30

! Add a line plot
call plt%add_plot(x, f, label, linestyle, linewidth, markersize, &
                  xscale, yscale, xlim, ylim)
!   x, f: real(dp)(:), required
!   label, linestyle: character, required (e.g. "-", "--", "o")
!   linewidth, markersize: integer, optional
!   xscale, yscale: character, optional (e.g. "log")
!   xlim, ylim: real(dp)(2), optional

! Add a contour plot
call plt%add_contour(x, y, f, levels, colorbar, filled, cmap)
!   x: real(dp)(:), y: real(dp)(:), f: real(dp)(:,:), required
!   levels: integer, optional (number of contour levels)
!   colorbar: logical, optional
!   filled: logical, optional
!   cmap: character, optional

! Add a line colored by a third variable
call plt%add_colored_line(x, y, color, cmap, clabel, linewidth)

! Show the figure
call plt%show()
```

## Minimal contour example

```fortran
program my_plot
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace
    implicit none

    integer, parameter :: n = 50
    real(dp), dimension(n) :: x, y
    real(dp), dimension(n, n) :: f
    type(myplot) :: plt
    integer :: i, j

    call linspace(0.0_dp, 2.0_dp*pi, n, x)
    call linspace(0.0_dp, 2.0_dp*pi, n, y)
    do j = 1, n
        do i = 1, n
            f(i, j) = sin(x(i))*cos(y(j))
        end do
    end do

    call plt%initialize(xlabel="x", ylabel="y", title="Example")
    call plt%add_contour(x, y, f, levels=20, colorbar=.true., filled=.true.)
    call plt%show()
end program my_plot
```
