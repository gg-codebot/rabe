module netcdf_mod
    use constants, only: dp

    implicit none
    private

    public :: netcdf_t

    integer, parameter :: max_vars = 100
    integer, parameter :: max_dims = 20

    type :: var_info_t
        character(len=100) :: name
        integer :: var_id
    end type var_info_t

    type :: dim_info_t
        character(len=100) :: name
        integer :: dim_id = -1
    end type dim_info_t

    type :: netcdf_t
        integer :: ncid = -1
        logical :: is_open = .false.
        logical :: in_define_mode = .false.
        integer :: n_vars = 0
        type(var_info_t) :: vars(max_vars)
        integer :: n_dims = 0
        type(dim_info_t) :: dims(max_dims)
    contains
        procedure :: create => netcdf_create
        procedure :: open => netcdf_open
        procedure :: add_global_attribute => netcdf_add_global_attr
        procedure :: read_global_attribute => netcdf_read_global_attr
        procedure :: def_dim => netcdf_def_dim
        procedure :: add_real => netcdf_add_real
        procedure :: add_real_1d => netcdf_add_real_1d
        procedure :: add_attr => netcdf_add_attr
        procedure :: read_attr => netcdf_read_attr
        procedure :: add_int_1d => netcdf_add_int_1d
        procedure :: end_define => netcdf_end_define
        procedure :: write_real => netcdf_write_real
        procedure :: write_real_1d => netcdf_write_real_1d
        procedure :: write_int_1d => netcdf_write_int_1d
        procedure :: read_real => netcdf_read_real
        procedure :: read_real_1d => netcdf_read_real_1d
        procedure :: read_int_1d => netcdf_read_int_1d
        procedure :: close => netcdf_close
        final :: netcdf_final
    end type netcdf_t

contains

    subroutine netcdf_create(this, filename)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: status

        if (this%is_open) then
            call this%close()
        end if

        status = nf90_create(filename, ior(NF90_CLOBBER, NF90_CLASSIC_MODEL), this%ncid)
        call check_netcdf_status(status, "creating file: "//filename)

        this%is_open = .true.
        this%in_define_mode = .true.
        this%n_vars = 0
    end subroutine netcdf_create

    subroutine netcdf_add_global_attr(this, attr_name, attr_value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: attr_name
        character(len=*), intent(in) :: attr_value
        integer :: status

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            error stop "Not in define mode"
        end if

        status = nf90_put_att(this%ncid, NF90_GLOBAL, attr_name, attr_value)
        call check_netcdf_status(status, "setting global attribute: " &
                                 //attr_name)
    end subroutine netcdf_add_global_attr

    subroutine netcdf_def_dim(this, dim_name, dim_length)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: dim_name
        integer, intent(in) :: dim_length
        integer :: status, dim_id

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            error stop "Not in define mode"
        end if

        if (this%n_dims >= max_dims) then
            error stop "Maximum number of dimensions exceeded"
        end if

        status = nf90_def_dim(this%ncid, dim_name, dim_length, dim_id)
        call check_netcdf_status(status, "defining dimension "//dim_name)

        this%n_dims = this%n_dims + 1
        this%dims(this%n_dims)%name = dim_name
        this%dims(this%n_dims)%dim_id = dim_id
    end subroutine

    subroutine netcdf_add_real(this, var_name)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        integer :: status, var_id

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            error stop "Not in define mode"
        end if

        if (this%n_vars >= max_vars) then
            error stop "Maximum number of variables exceeded"
        end if

        status = nf90_def_var(this%ncid, var_name, NF90_DOUBLE, varid=var_id)
        call check_netcdf_status(status, "defining variable: "//var_name)

        this%n_vars = this%n_vars + 1
        this%vars(this%n_vars)%name = var_name
        this%vars(this%n_vars)%var_id = var_id
    end subroutine netcdf_add_real

    subroutine netcdf_add_real_1d(this, var_name, dim_name)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        character(len=*), intent(in) :: dim_name
        integer :: status, var_id

        integer :: i, j
        integer :: dim_id

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            error stop "Not in define mode"
        end if

        if (this%n_vars >= max_vars) then
            error stop "Maximum number of variables exceeded"
        end if

        dim_id = -1
        do j = 1, this%n_dims
            if (trim(this%dims(j)%name) == trim(dim_name)) then
                dim_id = this%dims(j)%dim_id
                exit
            end if
        end do
        if (dim_id == -1) then
            error stop "Dimension not found: "//dim_name
        end if

        status = nf90_def_var(this%ncid, var_name, NF90_DOUBLE, [dim_id], var_id)
        call check_netcdf_status(status, "defining variable: "//var_name)

        this%n_vars = this%n_vars + 1
        this%vars(this%n_vars)%name = var_name
        this%vars(this%n_vars)%var_id = var_id
    end subroutine netcdf_add_real_1d

    subroutine netcdf_add_int_1d(this, var_name, dim_name)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        character(len=*), intent(in) :: dim_name
        integer :: status, var_id

        integer :: j
        integer :: dim_id

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            error stop "Not in define mode"
        end if

        if (this%n_vars >= max_vars) then
            error stop "Maximum number of variables exceeded"
        end if

        dim_id = -1
        do j = 1, this%n_dims
            if (trim(this%dims(j)%name) == trim(dim_name)) then
                dim_id = this%dims(j)%dim_id
                exit
            end if
        end do
        if (dim_id == -1) then
            error stop "Dimension not found: "//dim_name
        end if

        status = nf90_def_var(this%ncid, var_name, NF90_INT, &
                              [dim_id], var_id)
        call check_netcdf_status(status, &
                                 "defining variable: "//var_name)

        this%n_vars = this%n_vars + 1
        this%vars(this%n_vars)%name = var_name
        this%vars(this%n_vars)%var_id = var_id
    end subroutine netcdf_add_int_1d

    subroutine netcdf_add_attr(this, var_name, attr_name, &
                               attr_value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        character(len=*), intent(in) :: attr_name
        character(len=*), intent(in) :: attr_value
        integer :: status, var_id, i

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            error stop "Not in define mode"
        end if

        var_id = -1
        do i = 1, this%n_vars
            if (trim(this%vars(i)%name) == trim(var_name)) then
                var_id = this%vars(i)%var_id
                exit
            end if
        end do

        if (var_id == -1) then
            error stop "Variable not found: "//var_name
        end if

        status = nf90_put_att(this%ncid, var_id, attr_name, attr_value)
        call check_netcdf_status(status, "setting attribute "//attr_name &
                                 //" for variable: "//var_name)
    end subroutine netcdf_add_attr

    subroutine netcdf_end_define(this)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        integer :: status

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            return
        end if

        status = nf90_enddef(this%ncid)
        call check_netcdf_status(status, "ending definition mode")

        this%in_define_mode = .false.
    end subroutine netcdf_end_define

    subroutine netcdf_write_real(this, var_name, value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        real(dp), intent(in) :: value
        integer :: status, var_id, i

        if (.not. this%is_open) then
            error stop "NetCDF file not open for writing"
        end if

        if (this%in_define_mode) then
            call this%end_define()
        end if

        var_id = -1
        do i = 1, this%n_vars
            if (trim(this%vars(i)%name) == trim(var_name)) then
                var_id = this%vars(i)%var_id
                exit
            end if
        end do

        if (var_id == -1) then
            error stop "Variable not found: "//var_name
        end if

        status = nf90_put_var(this%ncid, var_id, value)
        call check_netcdf_status(status, "writing variable: "//var_name)
    end subroutine netcdf_write_real

    subroutine netcdf_write_real_1d(this, var_name, value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        real(dp), dimension(:), intent(in) :: value
        integer :: status, var_id, i

        if (.not. this%is_open) then
            error stop "NetCDF file not open for writing"
        end if

        if (this%in_define_mode) then
            call this%end_define()
        end if

        var_id = -1
        do i = 1, this%n_vars
            if (trim(this%vars(i)%name) == trim(var_name)) then
                var_id = this%vars(i)%var_id
                exit
            end if
        end do

        if (var_id == -1) then
            error stop "Variable not found: "//var_name
        end if

        status = nf90_put_var(this%ncid, var_id, value)
        call check_netcdf_status(status, "writing 1D array variable: "//var_name)
    end subroutine netcdf_write_real_1d

    subroutine netcdf_write_int_1d(this, var_name, value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        integer, dimension(:), intent(in) :: value
        integer :: status, var_id, i

        if (.not. this%is_open) then
            error stop "NetCDF file not open for writing"
        end if

        if (this%in_define_mode) then
            call this%end_define()
        end if

        var_id = -1
        do i = 1, this%n_vars
            if (trim(this%vars(i)%name) == trim(var_name)) then
                var_id = this%vars(i)%var_id
                exit
            end if
        end do

        if (var_id == -1) then
            error stop "Variable not found: "//var_name
        end if

        status = nf90_put_var(this%ncid, var_id, value)
        call check_netcdf_status(status, &
                                 "writing integer array: "//var_name)
    end subroutine netcdf_write_int_1d

    subroutine netcdf_close(this)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        integer :: status

        if (this%is_open) then
            if (this%in_define_mode) then
                call this%end_define()
            end if

            status = nf90_close(this%ncid)
            call check_netcdf_status(status, "closing NetCDF file")
            this%is_open = .false.
            this%in_define_mode = .false.
            this%ncid = -1
            this%n_vars = 0
        end if
    end subroutine netcdf_close

    subroutine netcdf_open(this, filename)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: status

        if (this%is_open) then
            call this%close()
        end if

        status = nf90_open(filename, NF90_NOWRITE, this%ncid)
        call check_netcdf_status(status, "opening file: "//filename)

        this%is_open = .true.
        this%in_define_mode = .false.
        this%n_vars = 0
    end subroutine netcdf_open

    subroutine netcdf_read_real(this, var_name, value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        real(dp), intent(out) :: value
        integer :: status, var_id

        if (.not. this%is_open) then
            error stop "NetCDF file not open for reading"
        end if

        status = nf90_inq_varid(this%ncid, var_name, var_id)
        call check_netcdf_status(status, "finding variable: "//var_name)

        status = nf90_get_var(this%ncid, var_id, value)
        call check_netcdf_status(status, "reading variable: "//var_name)
    end subroutine netcdf_read_real

    subroutine netcdf_read_real_1d(this, var_name, value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        real(dp), dimension(:), intent(out) :: value
        integer :: status, var_id

        if (.not. this%is_open) then
            error stop "NetCDF file not open for reading"
        end if

        status = nf90_inq_varid(this%ncid, var_name, var_id)
        call check_netcdf_status(status, "finding variable: "//var_name)

        status = nf90_get_var(this%ncid, var_id, value)
        call check_netcdf_status(status, "reading variable: "//var_name)
    end subroutine netcdf_read_real_1d

    subroutine netcdf_read_int_1d(this, var_name, value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        integer, dimension(:), intent(out) :: value
        integer :: status, var_id

        if (.not. this%is_open) then
            error stop "NetCDF file not open for reading"
        end if

        status = nf90_inq_varid(this%ncid, var_name, var_id)
        call check_netcdf_status(status, &
                                 "finding variable: "//var_name)

        status = nf90_get_var(this%ncid, var_id, value)
        call check_netcdf_status(status, &
                                 "reading integer array: "//var_name)
    end subroutine netcdf_read_int_1d

    subroutine netcdf_read_global_attr(this, attr_name, attr_value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: attr_name
        character(len=*), intent(out) :: attr_value
        integer :: status

        if (.not. this%is_open) then
            error stop "NetCDF file not open for reading"
        end if

        status = nf90_get_att(this%ncid, NF90_GLOBAL, attr_name, attr_value)
        call check_netcdf_status(status, "reading global attribute: " &
                                 //attr_name)
    end subroutine netcdf_read_global_attr

    subroutine netcdf_read_attr(this, var_name, attr_name, attr_value)
        use netcdf
        class(netcdf_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        character(len=*), intent(in) :: attr_name
        character(len=*), intent(out) :: attr_value
        integer :: status, var_id

        if (.not. this%is_open) then
            error stop "NetCDF file not open for reading"
        end if

        status = nf90_inq_varid(this%ncid, var_name, var_id)
        call check_netcdf_status(status, "finding variable: "//var_name)

        status = nf90_get_att(this%ncid, var_id, attr_name, attr_value)
        call check_netcdf_status(status, "reading attribute "//attr_name &
                                 //" for variable: "//var_name)
    end subroutine netcdf_read_attr

    subroutine netcdf_final(this)
        type(netcdf_t), intent(inout) :: this

        if (this%is_open) then
            call this%close()
        end if
    end subroutine netcdf_final

    subroutine check_netcdf_status(status, operation)
        use netcdf
        integer, intent(in) :: status
        character(len=*), intent(in) :: operation

        if (status /= NF90_NOERR) then
            print *, "NetCDF error during ", operation, ":"
            print *, trim(nf90_strerror(status))
            error stop "NetCDF operation failed"
        end if
    end subroutine check_netcdf_status

end module netcdf_mod
