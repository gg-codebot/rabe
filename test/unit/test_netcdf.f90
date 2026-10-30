program test_netcdf
    use constants, only: dp
    use netcdf_mod, only: netcdf_t
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-12_dp

    type(netcdf_t) :: nc_out, nc_in
    character(len=*), parameter :: test_file = "test_output.nc"
    real(dp), parameter :: test_factor_a = 1.23456789_dp
    real(dp), parameter :: test_factor_b = 9.87654321_dp
    character(len=*), parameter :: dim_name = "radius"
    integer, parameter :: n_dim = 3
    real(dp), parameter :: test_factor_a_array(n_dim) = [1.23_dp, 4.56_dp, 7.89_dp]
    real(dp), parameter :: test_factor_b_array(n_dim) = [98.7_dp, 6.54_dp, 32.1_dp]
    integer, parameter :: test_int_array(n_dim) = [1, 0, -102]
    character(len=*), parameter :: dim_name2 = "theta"
    integer, parameter :: n_dim2 = 2
    real(dp), parameter :: test_2d(n_dim, n_dim2) = reshape( &
                           [1.1_dp, 2.2_dp, 3.3_dp, 4.4_dp, 5.5_dp, 6.6_dp], &
                           [n_dim, n_dim2])
    character(len=*), parameter :: dim_name3 = "phi"
    integer, parameter :: n_dim3 = 4
    real(dp), parameter :: test_3d(n_dim, n_dim2, n_dim3) = reshape( &
                           [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, &
                            7._dp, 8._dp, 9._dp, 10._dp, 11._dp, 12._dp, &
                            13._dp, 14._dp, 15._dp, 16._dp, 17._dp, 18._dp, &
                            19._dp, 20._dp, 21._dp, 22._dp, 23._dp, 24._dp], &
                           [n_dim, n_dim2, n_dim3])
    character(len=*), parameter :: title = "RABE Bootstrap Current Analysis Results"
    character(len=*), parameter :: long_name_a = "1/sqrt(nu_star) factor"
    character(len=*), parameter :: long_name_b = "1/nu_star factor"
    character(len=*), parameter :: long_name_flags = "status flags"

    real(dp) :: read_factor_a, read_factor_b
    real(dp), dimension(n_dim) :: read_factor_a_array
    real(dp), dimension(n_dim) :: read_factor_b_array
    character(len=100) :: read_title, read_long_name_a, read_long_name_b
    character(len=100) :: read_long_name_flags

    integer, dimension(n_dim) :: read_int_array
    real(dp), dimension(n_dim, n_dim2) :: read_2d
    real(dp), dimension(n_dim, n_dim2, n_dim3) :: read_3d

    logical :: file_exists
    integer :: unit, iostat

    call nc_out%create(test_file)
    call nc_out%add_global_attribute("title", title)
    call nc_out%add_real("off_factor_a")
    call nc_out%add_attr("off_factor_a", "long_name", long_name_a)
    call nc_out%add_real("off_factor_b")
    call nc_out%add_attr("off_factor_b", "long_name", long_name_b)
    call nc_out%write_real("off_factor_a", test_factor_a)
    call nc_out%write_real("off_factor_b", test_factor_b)
    call nc_out%close()

    inquire (file=test_file, exist=file_exists)
    if (.not. file_exists) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: file does not exist after creation"
        error stop
    end if

    call nc_in%open(test_file)
    call nc_in%read_real("off_factor_a", read_factor_a)
    call nc_in%read_real("off_factor_b", read_factor_b)
    call nc_in%read_global_attribute("title", read_title)
    call nc_in%read_attr("off_factor_a", "long_name", read_long_name_a)
    call nc_in%read_attr("off_factor_b", "long_name", read_long_name_b)
    call nc_in%close()

    if (not_same(test_factor_a, read_factor_a, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_a mismatch"
        print *, "found: ", read_factor_a
        print *, "expected: ", test_factor_a
        error stop
    end if

    if (not_same(test_factor_b, read_factor_b, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_b mismatch"
        print *, "found: ", read_factor_b
        print *, "expected: ", test_factor_b
        error stop
    end if

    if (trim(read_title) /= title) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: global title mismatch"
        print *, "found: ", trim(read_title)
        print *, "expected: RABE Bootstrap Current Analysis Results"
        error stop
    end if

    if (trim(read_long_name_a) /= long_name_a) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_a long_name mismatch"
        print *, "found: ", trim(read_long_name_a)
        print *, "expected: ", trim(long_name_a)
        error stop
    end if

    if (trim(read_long_name_b) /= long_name_b) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_b long_name mismatch"
        print *, "found: ", trim(read_long_name_b)
        print *, "expected: ", trim(long_name_b)
        error stop
    end if

    inquire (file=test_file, exist=file_exists)
    if (file_exists) then
        open (newunit=unit, file=test_file, iostat=iostat)
        if (iostat == 0) then
            close (unit, status="delete")
        end if
    else
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: file does not exist after reading"
        error stop
    end if

    call nc_out%create(test_file)
    call nc_out%add_global_attribute("title", title)
    call nc_out%def_dim(dim_name, n_dim)
    call nc_out%add_real_1d("off_factor_a", dim_name)
    call nc_out%add_attr("off_factor_a", "long_name", long_name_a)
    call nc_out%add_real_1d("off_factor_b", dim_name)
    call nc_out%add_attr("off_factor_b", "long_name", long_name_b)
    call nc_out%add_int_1d("flags", dim_name)
    call nc_out%add_attr("flags", "long_name", long_name_flags)
    call nc_out%write_real_1d("off_factor_a", test_factor_a_array)
    call nc_out%write_real_1d("off_factor_b", test_factor_b_array)
    call nc_out%write_int_1d("flags", test_int_array)
    call nc_out%close()

    inquire (file=test_file, exist=file_exists)
    if (.not. file_exists) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: file does not exist after creation"
        error stop
    end if

    call nc_in%open(test_file)
    call nc_in%read_real_1d("off_factor_a", read_factor_a_array)
    call nc_in%read_real_1d("off_factor_b", read_factor_b_array)
    call nc_in%read_int_1d("flags", read_int_array)
    call nc_in%read_global_attribute("title", read_title)
    call nc_in%read_attr("off_factor_a", "long_name", read_long_name_a)
    call nc_in%read_attr("off_factor_b", "long_name", read_long_name_b)
    call nc_in%read_attr("flags", "long_name", read_long_name_flags)
    call nc_in%close()

    if (not_same(test_factor_a_array, read_factor_a_array, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_a_array mismatch"
        print *, "found: ", read_factor_a_array
        print *, "expected: ", test_factor_a_array
        error stop
    end if

    if (not_same(test_factor_b_array, read_factor_b_array, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_b mismatch"
        print *, "found: ", read_factor_b_array
        print *, "expected: ", test_factor_b_array
        error stop
    end if

    if (trim(read_title) /= title) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: global title mismatch"
        print *, "found: ", trim(read_title)
        print *, "expected: RABE Bootstrap Current Analysis Results"
        error stop
    end if

    if (trim(read_long_name_a) /= long_name_a) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_a long_name mismatch"
        print *, "found: ", trim(read_long_name_a)
        print *, "expected: ", trim(long_name_a)
        error stop
    end if

    if (trim(read_long_name_b) /= long_name_b) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_b long_name mismatch"
        print *, "found: ", trim(read_long_name_b)
        print *, "expected: ", trim(long_name_b)
        error stop
    end if

    if (any(read_int_array /= test_int_array)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: int array mismatch"
        print *, "found: ", read_int_array
        print *, "expected: ", test_int_array
        error stop
    end if

    if (trim(read_long_name_flags) /= long_name_flags) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: flags long_name mismatch"
        print *, "found: ", trim(read_long_name_flags)
        print *, "expected: ", trim(long_name_flags)
        error stop
    end if

    inquire (file=test_file, exist=file_exists)
    if (file_exists) then
        open (newunit=unit, file=test_file, iostat=iostat)
        if (iostat == 0) then
            close (unit, status="delete")
        end if
    else
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: file does not exist after reading"
        error stop
    end if

    ! --- 2D array round-trip test ---
    call nc_out%create(test_file)
    call nc_out%def_dim(dim_name, n_dim)
    call nc_out%def_dim(dim_name2, n_dim2)
    call nc_out%add_real_2d("matrix", dim_name, dim_name2)
    call nc_out%write_real_2d("matrix", test_2d)
    call nc_out%close()

    call nc_in%open(test_file)
    call nc_in%read_real_2d("matrix", read_2d)
    call nc_in%close()

    if (any(abs(test_2d - read_2d) > reltol*abs(test_2d))) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: 2D array mismatch"
        print *, "found: ", read_2d
        print *, "expected: ", test_2d
        error stop
    end if

    inquire (file=test_file, exist=file_exists)
    if (file_exists) then
        open (newunit=unit, file=test_file, iostat=iostat)
        if (iostat == 0) then
            close (unit, status="delete")
        end if
    end if

    ! --- 3D array round-trip test ---
    call nc_out%create(test_file)
    call nc_out%def_dim(dim_name, n_dim)
    call nc_out%def_dim(dim_name2, n_dim2)
    call nc_out%def_dim(dim_name3, n_dim3)
    call nc_out%add_real_3d("tensor", dim_name, dim_name2, dim_name3)
    call nc_out%write_real_3d("tensor", test_3d)
    call nc_out%close()

    call nc_in%open(test_file)
    call nc_in%read_real_3d("tensor", read_3d)
    call nc_in%close()

    if (any(abs(test_3d - read_3d) > reltol*abs(test_3d))) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: 3D array mismatch"
        print *, "found: ", read_3d
        print *, "expected: ", test_3d
        error stop
    end if

    inquire (file=test_file, exist=file_exists)
    if (file_exists) then
        open (newunit=unit, file=test_file, iostat=iostat)
        if (iostat == 0) then
            close (unit, status="delete")
        end if
    end if

end program test_netcdf
