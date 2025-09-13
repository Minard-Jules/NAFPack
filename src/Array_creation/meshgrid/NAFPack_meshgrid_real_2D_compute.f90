submodule(NAFPack_meshgrid:NAFPack_meshgrid_real_2D) NAFPack_meshgrid_real_2D_compute

    implicit none(type, external)

contains

    !=================================================================================
    ! Compute the meshgrid of two real vectors in 32-bit (single precision)
    !=================================================================================

    module subroutine compute_meshgrid_real_2D_sp( &
        x_vector, y_vector, &
        X, Y, &
        Nx, Ny, &
        loop_method)
        real(sp), dimension(:), intent(in) :: x_vector, y_vector
        real(sp), dimension(:, :), allocatable, intent(out) :: X, Y
        integer(isp), intent(in) :: Nx, Ny
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx), Y(Ny, Nx))

        if (loop_method%use_do_classic) then
            call compute_do_classic_real_2D_sp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_real_2D_sp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_real_2D_sp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_real_2D_sp(x_vector, y_vector, X, Y, Nx, Ny, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_real_2D_sp

    pure subroutine compute_do_classic_real_2D_sp(x_vector, y_vector, X, Y, Nx, Ny)
        real(sp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        real(sp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do

        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_classic_real_2D_sp

    pure subroutine compute_do_vectorized_real_2D_sp(x_vector, y_vector, X, Y, Nx, Ny)
        real(sp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        real(sp), dimension(:, :), intent(out) :: X, Y

        X = spread(x_vector, 1, Ny)
        Y = spread(y_vector, 2, Nx)

    end subroutine compute_do_vectorized_real_2D_sp

    pure subroutine compute_do_concurrent_real_2D_sp(x_vector, y_vector, X, Y, Nx, Ny)
        real(sp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        real(sp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do concurrent(i=1:Ny)
            X(i, :) = x_vector(:)
        end do

        do concurrent(j=1:Nx)
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_concurrent_real_2D_sp

    subroutine compute_openmp_real_2D_sp(x_vector, y_vector, X, Y, Nx, Ny, threads)
        real(sp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny, threads
        real(sp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        !$omp parallel default(none) private(i, j) &
        !$omp& shared(X, Y, x_vector, y_vector, Nx, Ny) &
        !$omp& num_threads(threads)

        !$omp do
        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do
        !$omp end do

        !$omp do
        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do
        !$omp end do

        !$omp end parallel

    end subroutine compute_openmp_real_2D_sp

    !=================================================================================
    ! Compute the meshgrid of two real vectors in 64-bit (double precision)
    !=================================================================================

    module subroutine compute_meshgrid_real_2D_dp( &
        x_vector, y_vector, &
        X, Y, &
        Nx, Ny, &
        loop_method)
        real(dp), dimension(:), intent(in) :: x_vector, y_vector
        real(dp), dimension(:, :), allocatable, intent(out) :: X, Y
        integer(isp), intent(in) :: Nx, Ny
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx), Y(Ny, Nx))

        if (loop_method%use_do_classic) then
            call compute_do_classic_real_2D_dp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_real_2D_dp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_real_2D_dp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_real_2D_dp(x_vector, y_vector, X, Y, Nx, Ny, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_real_2D_dp

    pure subroutine compute_do_classic_real_2D_dp(x_vector, y_vector, X, Y, Nx, Ny)
        real(dp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        real(dp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do

        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_classic_real_2D_dp

    pure subroutine compute_do_vectorized_real_2D_dp(x_vector, y_vector, X, Y, Nx, Ny)
        real(dp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        real(dp), dimension(:, :), intent(out) :: X, Y

        X = spread(x_vector, 1, Ny)
        Y = spread(y_vector, 2, Nx)

    end subroutine compute_do_vectorized_real_2D_dp

    pure subroutine compute_do_concurrent_real_2D_dp(x_vector, y_vector, X, Y, Nx, Ny)
        real(dp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        real(dp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do concurrent(i=1:Ny)
            X(i, :) = x_vector(:)
        end do

        do concurrent(j=1:Nx)
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_concurrent_real_2D_dp

    subroutine compute_openmp_real_2D_dp(x_vector, y_vector, X, Y, Nx, Ny, threads)
        real(dp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny, threads
        real(dp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        !$omp parallel default(none) private(i, j) &
        !$omp& shared(X, Y, x_vector, y_vector, Nx, Ny) &
        !$omp& num_threads(threads)

        !$omp do
        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do
        !$omp end do

        !$omp do
        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do
        !$omp end do

        !$omp end parallel

    end subroutine compute_openmp_real_2D_dp

    !=================================================================================
    ! Compute the meshgrid of two real vectors in 128-bit (quadruple precision)
    !=================================================================================

    module subroutine compute_meshgrid_real_2D_qp( &
        x_vector, y_vector, &
        X, Y, &
        Nx, Ny, &
        loop_method)
        real(qp), dimension(:), intent(in) :: x_vector, y_vector
        real(qp), dimension(:, :), allocatable, intent(out) :: X, Y
        integer(isp), intent(in) :: Nx, Ny
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx), Y(Ny, Nx))

        if (loop_method%use_do_classic) then
            call compute_do_classic_real_2D_qp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_real_2D_qp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_real_2D_qp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_real_2D_qp(x_vector, y_vector, X, Y, Nx, Ny, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_real_2D_qp

    pure subroutine compute_do_classic_real_2D_qp(x_vector, y_vector, X, Y, Nx, Ny)
        real(qp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        real(qp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do

        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_classic_real_2D_qp

    pure subroutine compute_do_vectorized_real_2D_qp(x_vector, y_vector, X, Y, Nx, Ny)
        real(qp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        real(qp), dimension(:, :), intent(out) :: X, Y

        X = spread(x_vector, 1, Ny)
        Y = spread(y_vector, 2, Nx)

    end subroutine compute_do_vectorized_real_2D_qp

    pure subroutine compute_do_concurrent_real_2D_qp(x_vector, y_vector, X, Y, Nx, Ny)
        real(qp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        real(qp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do concurrent(i=1:Ny)
            X(i, :) = x_vector(:)
        end do

        do concurrent(j=1:Nx)
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_concurrent_real_2D_qp

    subroutine compute_openmp_real_2D_qp(x_vector, y_vector, X, Y, Nx, Ny, threads)
        real(qp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny, threads
        real(qp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        !$omp parallel default(none) private(i, j) &
        !$omp& shared(X, Y, x_vector, y_vector, Nx, Ny) &
        !$omp& num_threads(threads)

        !$omp do
        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do
        !$omp end do

        !$omp do
        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do
        !$omp end do

        !$omp end parallel

    end subroutine compute_openmp_real_2D_qp

end submodule NAFPack_meshgrid_real_2D_compute
