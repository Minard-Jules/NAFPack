submodule(NAFPack_meshgrid) NAFPack_meshgrid_integer_3D

    implicit none(type, external)

    interface
        module subroutine compute_meshgrid_integer_3D_i8( &
            x_vector, y_vector, z_vector, &
            X, Y, Z, &
            Nx, Ny, Nz, &
            loop_method)
            integer(i8), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            integer(i8), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            integer(isp), intent(in) :: Nx, Ny, Nz
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_integer_3D_i8
    end interface

    interface
        module subroutine compute_meshgrid_integer_3D_i16( &
            x_vector, y_vector, z_vector, &
            X, Y, Z, &
            Nx, Ny, Nz, &
            loop_method)
            integer(i16), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            integer(i16), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            integer(isp), intent(in) :: Nx, Ny, Nz
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_integer_3D_i16
    end interface

    interface
        module subroutine compute_meshgrid_integer_3D_isp( &
            x_vector, y_vector, z_vector, &
            X, Y, Z, &
            Nx, Ny, Nz, &
            loop_method)
            integer(isp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            integer(isp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            integer(isp), intent(in) :: Nx, Ny, Nz
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_integer_3D_isp
    end interface

    interface
        module subroutine compute_meshgrid_integer_3D_idp( &
            x_vector, y_vector, z_vector, &
            X, Y, Z, &
            Nx, Ny, Nz, &
            loop_method)
            integer(idp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            integer(idp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            integer(isp), intent(in) :: Nx, Ny, Nz
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_integer_3D_idp
    end interface

contains

    module subroutine meshgrid_integer_i8_3D( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        indexing, &
        strict_mode, &
        loop_method)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(i8), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: Nx, Ny, Nz

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(x_vector, 1)
        Ny = size(y_vector, 1)
        Nz = size(z_vector, 1)
        if (use_ij_indexing) then
            call compute_meshgrid_integer_3D_i8(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_integer_3D_i8(y_vector, x_vector, z_vector, Y, X, Z, Ny, Nx, Nz, loop_method_used)
        end if

    end subroutine meshgrid_integer_i8_3D

    module subroutine meshgrid_integer_i16_3D( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        indexing, &
        strict_mode, &
        loop_method)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(i16), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: Nx, Ny, Nz

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(x_vector, 1)
        Ny = size(y_vector, 1)
        Nz = size(z_vector, 1)
        if (use_ij_indexing) then
            call compute_meshgrid_integer_3D_i16(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_integer_3D_i16(y_vector, x_vector, z_vector, Y, X, Z, Ny, Nx, Nz, loop_method_used)
        end if

    end subroutine meshgrid_integer_i16_3D

    module subroutine meshgrid_integer_isp_3D( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        indexing, &
        strict_mode, &
        loop_method)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: Nx, Ny, Nz

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(x_vector, 1)
        Ny = size(y_vector, 1)
        Nz = size(z_vector, 1)
        if (use_ij_indexing) then
            call compute_meshgrid_integer_3D_isp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_integer_3D_isp(y_vector, x_vector, z_vector, Y, X, Z, Ny, Nx, Nz, loop_method_used)
        end if

    end subroutine meshgrid_integer_isp_3D

    module subroutine meshgrid_integer_idp_3D( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        indexing, &
        strict_mode, &
        loop_method)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(idp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: Nx, Ny, Nz

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(x_vector, 1)
        Ny = size(y_vector, 1)
        Nz = size(z_vector, 1)
        if (use_ij_indexing) then
            call compute_meshgrid_integer_3D_idp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_integer_3D_idp(y_vector, x_vector, z_vector, Y, X, Z, Ny, Nx, Nz, loop_method_used)
        end if

    end subroutine meshgrid_integer_idp_3D

end submodule NAFPack_meshgrid_integer_3D
