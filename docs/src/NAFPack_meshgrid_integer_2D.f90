submodule(NAFPack_meshgrid) NAFPack_meshgrid_integer_2D

    implicit none(type, external)

    interface
        module subroutine compute_meshgrid_integer_2D_i8( &
            x_vector, y_vector, &
            X, Y, &
            Nx, Ny, &
            loop_method)
            integer(i8), dimension(:), intent(in) :: x_vector, y_vector
            integer(i8), dimension(:, :), allocatable, intent(out) :: X, Y
            integer(isp), intent(in) :: Nx, Ny
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_integer_2D_i8
    end interface

    interface
        module subroutine compute_meshgrid_integer_2D_i16( &
            x_vector, y_vector, &
            X, Y, &
            Nx, Ny, &
            loop_method)
            integer(i16), dimension(:), intent(in) :: x_vector, y_vector
            integer(i16), dimension(:, :), allocatable, intent(out) :: X, Y
            integer(isp), intent(in) :: Nx, Ny
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_integer_2D_i16
    end interface

    interface
        module subroutine compute_meshgrid_integer_2D_isp( &
            x_vector, y_vector, &
            X, Y, &
            Nx, Ny, &
            loop_method)
            integer(isp), dimension(:), intent(in) :: x_vector, y_vector
            integer(isp), dimension(:, :), allocatable, intent(out) :: X, Y
            integer(isp), intent(in) :: Nx, Ny
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_integer_2D_isp
    end interface

    interface
        module subroutine compute_meshgrid_integer_2D_idp( &
            x_vector, y_vector, &
            X, Y, &
            Nx, Ny, &
            loop_method)
            integer(idp), dimension(:), intent(in) :: x_vector, y_vector
            integer(idp), dimension(:, :), allocatable, intent(out) :: X, Y
            integer(isp), intent(in) :: Nx, Ny
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_integer_2D_idp
    end interface

contains

    module subroutine meshgrid_integer_i8_2D( &
        x_vector, y_vector, &
        X, Y, &
        indexing, &
        strict_mode, &
        loop_method)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector
        integer(i8), dimension(:, :), allocatable, intent(out) :: X, Y
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: Nx, Ny

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(x_vector, 1)
        Ny = size(y_vector, 1)
        if (use_ij_indexing) then
            call compute_meshgrid_integer_2D_i8(x_vector, y_vector, X, Y, Nx, Ny, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_integer_2D_i8(y_vector, x_vector, Y, X, Ny, Nx, loop_method_used)
        end if

    end subroutine meshgrid_integer_i8_2D

    module subroutine meshgrid_integer_i16_2D( &
        x_vector, y_vector, &
        X, Y, &
        indexing, &
        strict_mode, &
        loop_method)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector
        integer(i16), dimension(:, :), allocatable, intent(out) :: X, Y
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: Nx, Ny

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(x_vector, 1)
        Ny = size(y_vector, 1)
        if (use_ij_indexing) then
            call compute_meshgrid_integer_2D_i16(x_vector, y_vector, X, Y, Nx, Ny, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_integer_2D_i16(y_vector, x_vector, Y, X, Ny, Nx, loop_method_used)
        end if

    end subroutine meshgrid_integer_i16_2D

    module subroutine meshgrid_integer_isp_2D( &
        x_vector, y_vector, &
        X, Y, &
        indexing, &
        strict_mode, &
        loop_method)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), dimension(:, :), allocatable, intent(out) :: X, Y
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: Nx, Ny

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(x_vector, 1)
        Ny = size(y_vector, 1)
        if (use_ij_indexing) then
            call compute_meshgrid_integer_2D_isp(x_vector, y_vector, X, Y, Nx, Ny, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_integer_2D_isp(y_vector, x_vector, Y, X, Ny, Nx, loop_method_used)
        end if

    end subroutine meshgrid_integer_isp_2D

    module subroutine meshgrid_integer_idp_2D( &
        x_vector, y_vector, &
        X, Y, &
        indexing, &
        strict_mode, &
        loop_method)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector
        integer(idp), dimension(:, :), allocatable, intent(out) :: X, Y
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: Nx, Ny

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(x_vector, 1)
        Ny = size(y_vector, 1)
        if (use_ij_indexing) then
            call compute_meshgrid_integer_2D_idp(x_vector, y_vector, X, Y, Nx, Ny, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_integer_2D_idp(y_vector, x_vector, Y, X, Ny, Nx, loop_method_used)
        end if

    end subroutine meshgrid_integer_idp_2D

end submodule NAFPack_meshgrid_integer_2D
