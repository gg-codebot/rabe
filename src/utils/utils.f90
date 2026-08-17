module utils
    use constants, only: dp

    implicit none

    interface not_same
        module procedure not_same_scalar
        module procedure not_same_array
        module procedure not_same_matrix
    end interface

contains
    subroutine linspace(a, b, n, x, include_endpoint)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n
        real(dp), dimension(:), intent(out) :: x
        logical, intent(in), optional :: include_endpoint

        real(dp) :: dx
        integer :: i

        if (n < 2) then
            print *, "Error in linspace: n must be at least 2."
            error stop
        end if
        if (size(x) < n) then
            print *, "Error in linspace: output array size must be at least n."
            error stop
        end if

        if (present(include_endpoint)) then
            if (include_endpoint) then
                dx = (b - a)/(n - 1)
            else
                dx = (b - a)/n
            end if
        else
            dx = (b - a)/(n - 1)
        end if

        do i = 1, n
            x(i) = a + (i - 1)*dx
        end do
    end subroutine linspace

    function not_same_scalar(scalar_1, scalar_2, reltol_in, abstol_in)
        use ieee_arithmetic, only: ieee_is_nan
        real(dp), intent(in) :: scalar_1, scalar_2
        real(dp), intent(in), optional :: reltol_in, abstol_in
        logical :: not_same_scalar

        real(dp) :: reltol, abstol

        if (present(reltol_in)) then
            reltol = reltol_in
        else
            reltol = 0.0_dp
        end if
        if (present(abstol_in)) then
            abstol = abstol_in
        else
            abstol = reltol
        end if

        if (ieee_is_nan(scalar_1) .or. ieee_is_nan(scalar_2)) then
            not_same_scalar = .true.
        else
            not_same_scalar = abs(scalar_1 - scalar_2) > (reltol*abs(scalar_1) + abstol)
        end if
    end function not_same_scalar

    function not_same_array(array_1, array_2, reltol_in, abstol_in)
        real(dp), dimension(:), intent(in) :: array_1, array_2
        real(dp), intent(in), optional :: reltol_in, abstol_in
        logical :: not_same_array

        real(dp) :: reltol, abstol

        if (present(reltol_in)) then
            reltol = reltol_in
        else
            reltol = 0.0_dp
        end if
        if (present(abstol_in)) then
            abstol = abstol_in
        else
            abstol = reltol
        end if

        if (contains_nan(array_1) .or. contains_nan(array_2)) then
            not_same_array = .true.
        else
            not_same_array = any(abs(array_1 - array_2) > &
                                 (reltol*abs(array_1) + abstol))
        end if
    end function not_same_array

    function contains_nan(array)
        use ieee_arithmetic, only: ieee_is_nan
        implicit none
        real(dp), intent(in) :: array(:)
        logical :: contains_nan
        integer :: i

        contains_nan = .false.
        do i = 1, size(array)
            if (ieee_is_nan(array(i))) then
                contains_nan = .true.
                return
            end if
        end do
    end function contains_nan

    function not_same_matrix(matrix_1, matrix_2, reltol_in, abstol_in)
        real(dp), dimension(:, :), intent(in) :: matrix_1, matrix_2
        real(dp), intent(in), optional :: reltol_in, abstol_in
        logical :: not_same_matrix

        real(dp) :: reltol, abstol
        integer :: rows_1, cols_1, col
        integer :: rows_2, cols_2

        if (present(reltol_in)) then
            reltol = reltol_in
        else
            reltol = 0.0_dp
        end if
        if (present(abstol_in)) then
            abstol = abstol_in
        else
            abstol = reltol
        end if

        not_same_matrix = .false.

        rows_1 = size(matrix_1, 1)
        cols_1 = size(matrix_1, 2)
        rows_2 = size(matrix_2, 1)
        cols_2 = size(matrix_2, 2)
        if (rows_1 /= rows_2 .or. cols_1 /= cols_2) then
            print *, "Error in not_same_matrix: matrices must have the same dimensions."
            print *, "Matrix 1 dimensions: ", rows_1, "x", cols_1
            print *, "Matrix 2 dimensions: ", rows_2, "x", cols_2
            error stop
        end if
        do col = 1, cols_1
            if (not_same_array(matrix_1(:, col), matrix_2(:, col), reltol, abstol)) then
                not_same_matrix = .true.
                return
            end if
        end do
    end function not_same_matrix

end module utils
