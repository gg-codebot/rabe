module myplot_module
    use pyplot_module, only: pyplot
    use constants, only: dp

    implicit none

    type :: myplot
        type(pyplot) :: plt
    contains
        procedure :: initialize
        procedure :: add_plot
        procedure :: add_colored_line
        procedure :: add_contour
        procedure :: show
    end type myplot

contains

    subroutine initialize(self, xlabel, ylabel, legend, figsize, title, fontsize)
        class(myplot), intent(out) :: self

        character(len=*), intent(in), optional :: xlabel, ylabel
        logical, intent(in), optional :: legend
        integer, dimension(2), intent(in), optional :: figsize
        character(len=*), intent(in), optional :: title
        integer, intent(in), optional :: fontsize

        integer :: fontsize_internal
        integer, dimension(2) :: figsize_internal
        logical, parameter :: tight_layout = .true.

        if (present(fontsize)) then
            fontsize_internal = fontsize
        else
            fontsize_internal = 30
        end if

        if (present(figsize)) then
            figsize_internal = figsize
        else
            figsize_internal = [10, 8]
        end if

        call self%plt%initialize(grid=.true., &
                                 xlabel=xlabel, &
                                 ylabel=ylabel, &
                                 legend=legend, &
                                 figsize=figsize_internal, &
                                 title=title, &
                                 raw_strings=.true., &
                                 font_size=fontsize_internal, &
                                 legend_fontsize=fontsize_internal, &
                                 axes_labelsize=fontsize_internal, &
                                 xtick_labelsize=fontsize_internal, &
                                 ytick_labelsize=fontsize_internal, &
                                 tight_layout=tight_layout)
    end subroutine initialize

    subroutine add_plot(self, x, f, label, linestyle, linewidth, markersize, &
                        xscale, yscale, xlim, ylim)
        class(myplot), intent(inout) :: self
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(in) :: f
        character(len=*), intent(in) :: label, linestyle

        integer, intent(in), optional :: linewidth
        integer, intent(in), optional :: markersize
        character(len=*), intent(in), optional :: xscale, yscale
        real(dp), dimension(2), optional :: xlim, ylim

        integer :: my_linewidth

        if (present(linewidth)) then
            my_linewidth = linewidth
        else
            my_linewidth = 1
        end if
        call self%plt%add_plot(x, f, label, linestyle, &
                               linewidth=my_linewidth, &
                               markersize=markersize, &
                               xscale=xscale, &
                               yscale=yscale, &
                               xlim=xlim, &
                               ylim=ylim)
    end subroutine add_plot

    subroutine add_colored_line(self, x, y, color, cmap, clabel, linewidth)
        class(myplot), intent(inout) :: self
        real(dp), dimension(:), intent(in) :: x, y, color
        character(len=*), intent(in), optional :: cmap, clabel
        integer, intent(in), optional :: linewidth

        character(len=:), allocatable :: cmap_str
        integer :: lw

        if (present(cmap)) then
            cmap_str = cmap
        else
            cmap_str = "viridis"
        end if
        if (present(linewidth)) then
            lw = linewidth
        else
            lw = 1
        end if

        call self%plt%add_str( &
            'from matplotlib.collections import LineCollection')
        call self%plt%add_str( &
            '_x = np.array('//array_to_string(x)//')')
        call self%plt%add_str( &
            '_y = np.array('//array_to_string(y)//')')
        call self%plt%add_str( &
            '_c = np.array('//array_to_string(color)//')')
        call self%plt%add_str( &
            '_pts = np.array([_x, _y]).T.reshape(-1, 1, 2)')
        call self%plt%add_str( &
            '_seg = np.concatenate([_pts[:-1], _pts[1:]], axis=1)')
        call self%plt%add_str( &
            '_mask = ~np.isnan(_c[:-1]) ' &
            //'& ~np.isnan(_x[:-1]) & ~np.isnan(_x[1:]) ' &
            //'& ~np.isnan(_y[:-1]) & ~np.isnan(_y[1:])')
        call self%plt%add_str( &
            '_lc = LineCollection(_seg[_mask], cmap="' &
            //trim(cmap_str)//'")')
        call self%plt%add_str('_lc.set_array(_c[:-1][_mask])')
        write (cmap_str, '(I0)') lw
        call self%plt%add_str( &
            '_lc.set_linewidth('//trim(cmap_str)//')')
        call self%plt%add_str('ax.add_collection(_lc)')
        call self%plt%add_str('ax.autoscale_view()')
        if (present(clabel)) then
            call self%plt%add_str( &
                '_cb = fig.colorbar(_lc, ax=ax)')
            call self%plt%add_str( &
                '_cb.set_label(r"'//trim(clabel)//'")')
        end if

    end subroutine add_colored_line

    function array_to_string(arr) result(str)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        real(dp), dimension(:), intent(in) :: arr
        character(len=:), allocatable :: str

        character(len=25) :: val_str
        character(len=:), allocatable :: buf
        integer :: i, p, l

        allocate (character(len=26*size(arr) + 2) :: buf)
        buf(1:1) = '['
        p = 1
        do i = 1, size(arr)
            if (i > 1) then
                buf(p + 1:p + 1) = ','
                p = p + 1
            end if
            if (ieee_is_nan(arr(i))) then
                buf(p + 1:p + 6) = 'np.nan'
                p = p + 6
            else
                write (val_str, '(ES23.15E3)') arr(i)
                val_str = adjustl(val_str)
                l = len_trim(val_str)
                buf(p + 1:p + l) = val_str(1:l)
                p = p + l
            end if
        end do
        buf(p + 1:p + 1) = ']'
        str = buf(1:p + 1)
    end function array_to_string

    subroutine add_contour(self, x, y, f, levels, colorbar, filled, cmap)
        use utils, only: linspace
        class(myplot), intent(inout) :: self
        real(dp), dimension(:), intent(in) :: x, y
        real(dp), dimension(:, :), intent(in) :: f
        integer, optional :: levels
        logical, optional :: colorbar
        logical, optional :: filled
        character(len=*), intent(in), optional :: cmap

        real(dp), dimension(:), allocatable :: levels_pyplot

        if (present(levels)) then
            allocate (levels_pyplot(levels))
            call linspace(minval(f), maxval(f), levels, levels_pyplot)
            call self%plt%add_contour(x, y, f, linestyle="-", &
                                      levels=levels_pyplot, &
                                      colorbar=colorbar, &
                                      filled=filled, &
                                      cmap=cmap)
        else
            call self%plt%add_contour(x, y, f, linestyle="-", &
                                      colorbar=colorbar, &
                                      filled=filled, &
                                      cmap=cmap)
        end if
    end subroutine add_contour

    subroutine show(self)
        class(myplot), intent(inout) :: self

        call self%plt%showfig()
    end subroutine show

end module myplot_module
