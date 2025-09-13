submodule(NAFPack_meshgrid) NAFPack_meshgrid_real_2D

    implicit none(type, external)

    interface
        module subroutine compute_meshgrid_real_2D_sp( &
            x_vector, y_vector, &
            X, Y, &
            Nx, Ny, &
            loop_method)
            real(sp), dimension(:), intent(in) :: x_vector, y_vector
            real(sp), dimension(:, :), allocatable, intent(out) :: X, Y
            integer(isp), intent(in) :: Nx, Ny
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_real_2D_sp
    end interface

    interface
        module subroutine compute_meshgrid_real_2D_dp( &
            x_vector, y_vector, &
            X, Y, &
            Nx, Ny, &
            loop_method)
            real(dp), dimension(:), intent(in) :: x_vector, y_vector
            real(dp), dimension(:, :), allocatable, intent(out) :: X, Y
            integer(isp), intent(in) :: Nx, Ny
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_real_2D_dp
    end interface

    interface
        module subroutine compute_meshgrid_real_2D_qp( &
            x_vector, y_vector, &
            X, Y, &
            Nx, Ny, &
            loop_method)
            real(qp), dimension(:), intent(in) :: x_vector, y_vector
            real(qp), dimension(:, :), allocatable, intent(out) :: X, Y
            integer(isp), intent(in) :: Nx, Ny
            type(LoopMethod), intent(in) :: loop_method
        end subroutine compute_meshgrid_real_2D_qp
    end interface

contains

    module subroutine meshgrid_real_sp_2D( &
        x_vector, y_vector, &
        X, Y, &
        indexing, &
        strict_mode, &
        loop_method)
        real(sp), dimension(:), intent(in) :: x_vector, y_vector
        real(sp), dimension(:, :), allocatable, intent(out) :: X, Y
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
            call compute_meshgrid_real_2D_sp(x_vector, y_vector, X, Y, Nx, Ny, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_real_2D_sp(y_vector, x_vector, Y, X, Ny, Nx, loop_method_used)
        end if

    end subroutine meshgrid_real_sp_2D

    module subroutine meshgrid_real_dp_2D( &
        x_vector, y_vector, &
        X, Y, &
        indexing, &
        strict_mode, &
        loop_method)
        real(dp), dimension(:), intent(in) :: x_vector, y_vector
        real(dp), dimension(:, :), allocatable, intent(out) :: X, Y
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
            call compute_meshgrid_real_2D_dp(x_vector, y_vector, X, Y, Nx, Ny, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_real_2D_dp(y_vector, x_vector, Y, X, Ny, Nx, loop_method_used)
        end if

    end subroutine meshgrid_real_dp_2D

    module subroutine meshgrid_real_qp_2D( &
        x_vector, y_vector, &
        X, Y, &
        indexing, &
        strict_mode, &
        loop_method)
        real(qp), dimension(:), intent(in) :: x_vector, y_vector
        real(qp), dimension(:, :), allocatable, intent(out) :: X, Y
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
            call compute_meshgrid_real_2D_qp(x_vector, y_vector, X, Y, Nx, Ny, loop_method_used)
        else if (use_xy_indexing) then
            call compute_meshgrid_real_2D_qp(y_vector, x_vector, Y, X, Ny, Nx, loop_method_used)
        end if

    end subroutine meshgrid_real_qp_2D

end submodule NAFPack_meshgrid_real_2D
