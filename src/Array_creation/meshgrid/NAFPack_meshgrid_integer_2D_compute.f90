submodule(NAFPack_meshgrid:NAFPack_meshgrid_integer_2D) NAFPack_meshgrid_integer_2D_compute

    implicit none(type, external)

contains

    !=================================================================================
    ! Compute the meshgrid of two integer vectors in 8-bit
    !=================================================================================

    module subroutine compute_meshgrid_integer_2D_i8( &
        x_vector, y_vector, &
        X, Y, &
        Nx, Ny, &
        loop_method)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector
        integer(i8), dimension(:, :), allocatable, intent(out) :: X, Y
        integer(isp), intent(in) :: Nx, Ny
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx), Y(Ny, Nx))

        if (loop_method%use_do_classic) then
            call compute_do_classic_integer_2D_i8(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_integer_2D_i8(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_integer_2D_i8(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_integer_2D_i8(x_vector, y_vector, X, Y, Nx, Ny, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_integer_2D_i8

    pure subroutine compute_do_classic_integer_2D_i8(x_vector, y_vector, X, Y, Nx, Ny)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(i8), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do

        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_classic_integer_2D_i8

    pure subroutine compute_do_vectorized_integer_2D_i8(x_vector, y_vector, X, Y, Nx, Ny)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(i8), dimension(:, :), intent(out) :: X, Y

        X = spread(x_vector, 1, Ny)
        Y = spread(y_vector, 2, Nx)

    end subroutine compute_do_vectorized_integer_2D_i8

    pure subroutine compute_do_concurrent_integer_2D_i8(x_vector, y_vector, X, Y, Nx, Ny)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(i8), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do concurrent(i=1:Ny)
            X(i, :) = x_vector(:)
        end do

        do concurrent(j=1:Nx)
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_concurrent_integer_2D_i8

    subroutine compute_openmp_integer_2D_i8(x_vector, y_vector, X, Y, Nx, Ny, threads)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny, threads
        integer(i8), dimension(:, :), intent(out) :: X, Y
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

    end subroutine compute_openmp_integer_2D_i8

    !=================================================================================
    ! Compute the meshgrid of two integer vectors in 16-bit
    !=================================================================================

    module subroutine compute_meshgrid_integer_2D_i16( &
        x_vector, y_vector, &
        X, Y, &
        Nx, Ny, &
        loop_method)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector
        integer(i16), dimension(:, :), allocatable, intent(out) :: X, Y
        integer(isp), intent(in) :: Nx, Ny
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx), Y(Ny, Nx))

        if (loop_method%use_do_classic) then
            call compute_do_classic_integer_2D_i16(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_integer_2D_i16(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_integer_2D_i16(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_integer_2D_i16(x_vector, y_vector, X, Y, Nx, Ny, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_integer_2D_i16

    pure subroutine compute_do_classic_integer_2D_i16(x_vector, y_vector, X, Y, Nx, Ny)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(i16), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do

        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_classic_integer_2D_i16

    pure subroutine compute_do_vectorized_integer_2D_i16(x_vector, y_vector, X, Y, Nx, Ny)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(i16), dimension(:, :), intent(out) :: X, Y

        X = spread(x_vector, 1, Ny)
        Y = spread(y_vector, 2, Nx)

    end subroutine compute_do_vectorized_integer_2D_i16

    pure subroutine compute_do_concurrent_integer_2D_i16(x_vector, y_vector, X, Y, Nx, Ny)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(i16), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do concurrent(i=1:Ny)
            X(i, :) = x_vector(:)
        end do

        do concurrent(j=1:Nx)
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_concurrent_integer_2D_i16

    subroutine compute_openmp_integer_2D_i16(x_vector, y_vector, X, Y, Nx, Ny, threads)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny, threads
        integer(i16), dimension(:, :), intent(out) :: X, Y
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

    end subroutine compute_openmp_integer_2D_i16

    !=================================================================================
    ! Compute the meshgrid of two integer vectors in 32-bit (Single precision)
    !=================================================================================

    module subroutine compute_meshgrid_integer_2D_isp( &
        x_vector, y_vector, &
        X, Y, &
        Nx, Ny, &
        loop_method)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), dimension(:, :), allocatable, intent(out) :: X, Y
        integer(isp), intent(in) :: Nx, Ny
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx), Y(Ny, Nx))

        if (loop_method%use_do_classic) then
            call compute_do_classic_integer_2D_isp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_integer_2D_isp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_integer_2D_isp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_integer_2D_isp(x_vector, y_vector, X, Y, Nx, Ny, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_integer_2D_isp

    pure subroutine compute_do_classic_integer_2D_isp(x_vector, y_vector, X, Y, Nx, Ny)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(isp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do

        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_classic_integer_2D_isp

    pure subroutine compute_do_vectorized_integer_2D_isp(x_vector, y_vector, X, Y, Nx, Ny)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(isp), dimension(:, :), intent(out) :: X, Y

        X = spread(x_vector, 1, Ny)
        Y = spread(y_vector, 2, Nx)

    end subroutine compute_do_vectorized_integer_2D_isp

    pure subroutine compute_do_concurrent_integer_2D_isp(x_vector, y_vector, X, Y, Nx, Ny)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(isp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do concurrent(i=1:Ny)
            X(i, :) = x_vector(:)
        end do

        do concurrent(j=1:Nx)
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_concurrent_integer_2D_isp

    subroutine compute_openmp_integer_2D_isp(x_vector, y_vector, X, Y, Nx, Ny, threads)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny, threads
        integer(isp), dimension(:, :), intent(out) :: X, Y
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

    end subroutine compute_openmp_integer_2D_isp

    !=================================================================================
    ! Compute the meshgrid of two integer vectors in 64-bit (double precision)
    !=================================================================================

    module subroutine compute_meshgrid_integer_2D_idp( &
        x_vector, y_vector, &
        X, Y, &
        Nx, Ny, &
        loop_method)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector
        integer(idp), dimension(:, :), allocatable, intent(out) :: X, Y
        integer(isp), intent(in) :: Nx, Ny
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx), Y(Ny, Nx))

        if (loop_method%use_do_classic) then
            call compute_do_classic_integer_2D_idp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_integer_2D_idp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_integer_2D_idp(x_vector, y_vector, X, Y, Nx, Ny)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_integer_2D_idp(x_vector, y_vector, X, Y, Nx, Ny, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_integer_2D_idp

    pure subroutine compute_do_classic_integer_2D_idp(x_vector, y_vector, X, Y, Nx, Ny)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(idp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do i = 1, Ny
            X(i, :) = x_vector(:)
        end do

        do j = 1, Nx
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_classic_integer_2D_idp

    pure subroutine compute_do_vectorized_integer_2D_idp(x_vector, y_vector, X, Y, Nx, Ny)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(idp), dimension(:, :), intent(out) :: X, Y

        X = spread(x_vector, 1, Ny)
        Y = spread(y_vector, 2, Nx)

    end subroutine compute_do_vectorized_integer_2D_idp

    pure subroutine compute_do_concurrent_integer_2D_idp(x_vector, y_vector, X, Y, Nx, Ny)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny
        integer(idp), dimension(:, :), intent(out) :: X, Y
        integer(isp) :: i, j

        do concurrent(i=1:Ny)
            X(i, :) = x_vector(:)
        end do

        do concurrent(j=1:Nx)
            Y(:, j) = y_vector(:)
        end do

    end subroutine compute_do_concurrent_integer_2D_idp

    subroutine compute_openmp_integer_2D_idp(x_vector, y_vector, X, Y, Nx, Ny, threads)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), intent(in) :: Nx, Ny, threads
        integer(idp), dimension(:, :), intent(out) :: X, Y
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

    end subroutine compute_openmp_integer_2D_idp

end submodule NAFPack_meshgrid_integer_2D_compute
