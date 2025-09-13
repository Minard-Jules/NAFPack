submodule(NAFPack_meshgrid:NAFPack_meshgrid_complex_3D) NAFPack_meshgrid_complex_3D_compute

    implicit none(type, external)

contains

    !=================================================================================
    ! Compute the meshgrid of three complex vectors in 32-bit (single precision)
    !=================================================================================

    module subroutine compute_meshgrid_cmplx_3D_sp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, &
        loop_method)
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        complex(sp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        integer(isp), intent(in) :: Nx, Ny, Nz
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx, Nz), Y(Ny, Nx, Nz), Z(Ny, Nx, Nz))

        if (loop_method%use_do_classic) then
            call compute_do_classic_cmplx_3D_sp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_cmplx_3D_sp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_cmplx_3D_sp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_cmplx_3D_sp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_cmplx_3D_sp

    pure subroutine compute_do_classic_cmplx_3D_sp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        complex(sp), dimension(:, :, :), intent(out) :: X, Y, Z
        integer(isp) :: i, j, k

        do k = 1, Nz
            do i = 1, Ny
                X(i, :, k) = x_vector(:)
            end do
        end do

        do k = 1, Nz
            do j = 1, Nx
                Y(:, j, k) = y_vector(:)
            end do
        end do

        do j = 1, Nx
            do i = 1, Ny
                Z(i, j, :) = z_vector(:)
            end do
        end do

    end subroutine compute_do_classic_cmplx_3D_sp

    pure subroutine compute_do_vectorized_cmplx_3D_sp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        complex(sp), dimension(:, :, :), intent(out) :: X, Y, Z

        X = spread(spread(x_vector, 1, Ny), 3, Nz)
        Y = spread(spread(y_vector, 2, Nx), 3, Nz)
        Z = spread(spread(z_vector, 1, Nx), 1, Ny)

    end subroutine compute_do_vectorized_cmplx_3D_sp

    pure subroutine compute_do_concurrent_cmplx_3D_sp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        complex(sp), dimension(:, :, :), intent(out) :: X, Y, Z
        integer(isp) :: i, j, k

        do concurrent(k=1:Nz, i=1:Ny)
            X(i, :, k) = x_vector(:)
        end do

        do concurrent(k=1:Nz, j=1:Nx)
            Y(:, j, k) = y_vector(:)
        end do

        do concurrent(j=1:Nx, i=1:Ny)
            Z(i, j, :) = z_vector(:)
        end do

    end subroutine compute_do_concurrent_cmplx_3D_sp

    subroutine compute_openmp_cmplx_3D_sp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, threads)
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz, threads
        complex(sp), dimension(:, :, :), intent(out) :: X, Y, Z
        integer(isp) :: i, j, k

        !$omp parallel default(none) private(i, j, k) &
        !$omp& shared(X, Y, Z, x_vector, y_vector, z_vector, Nx, Ny, Nz) &
        !$omp& num_threads(threads)

        !$omp do
        do k = 1, Nz
            do i = 1, Ny
                X(i, :, k) = x_vector(:)
            end do
        end do
        !$omp end do

        !$omp do
        do k = 1, Nz
            do j = 1, Nx
                Y(:, j, k) = y_vector(:)
            end do
        end do
        !$omp end do

        !$omp do
        do j = 1, Nx
            do i = 1, Ny
                Z(i, j, :) = z_vector(:)
            end do
        end do
        !$omp end do

        !$omp end parallel

    end subroutine compute_openmp_cmplx_3D_sp

    !=================================================================================
    ! Compute the meshgrid of three complex vectors in 64-bit (double precision)
    !=================================================================================

    module subroutine compute_meshgrid_cmplx_3D_dp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, &
        loop_method)
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        complex(dp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        integer(isp), intent(in) :: Nx, Ny, Nz
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx, Nz), Y(Ny, Nx, Nz), Z(Ny, Nx, Nz))

        if (loop_method%use_do_classic) then
            call compute_do_classic_cmplx_3D_dp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_cmplx_3D_dp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_cmplx_3D_dp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_cmplx_3D_dp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_cmplx_3D_dp

    pure subroutine compute_do_classic_cmplx_3D_dp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        complex(dp), dimension(:, :, :), intent(out) :: X, Y, Z
        integer(isp) :: i, j, k

        do k = 1, Nz
            do i = 1, Ny
                X(i, :, k) = x_vector(:)
            end do
        end do

        do k = 1, Nz
            do j = 1, Nx
                Y(:, j, k) = y_vector(:)
            end do
        end do

        do j = 1, Nx
            do i = 1, Ny
                Z(i, j, :) = z_vector(:)
            end do
        end do

    end subroutine compute_do_classic_cmplx_3D_dp

    pure subroutine compute_do_vectorized_cmplx_3D_dp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        complex(dp), dimension(:, :, :), intent(out) :: X, Y, Z

        X = spread(spread(x_vector, 1, Ny), 3, Nz)
        Y = spread(spread(y_vector, 2, Nx), 3, Nz)
        Z = spread(spread(z_vector, 1, Nx), 1, Ny)

    end subroutine compute_do_vectorized_cmplx_3D_dp

    pure subroutine compute_do_concurrent_cmplx_3D_dp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        complex(dp), dimension(:, :, :), intent(out) :: X, Y, Z
        integer(isp) :: i, j, k

        do concurrent(k=1:Nz, i=1:Ny)
            X(i, :, k) = x_vector(:)
        end do

        do concurrent(k=1:Nz, j=1:Nx)
            Y(:, j, k) = y_vector(:)
        end do

        do concurrent(j=1:Nx, i=1:Ny)
            Z(i, j, :) = z_vector(:)
        end do

    end subroutine compute_do_concurrent_cmplx_3D_dp

    subroutine compute_openmp_cmplx_3D_dp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, threads)
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz, threads
        complex(dp), dimension(:, :, :), intent(out) :: X, Y, Z
        integer(isp) :: i, j, k

        !$omp parallel default(none) private(i, j, k) &
        !$omp& shared(X, Y, Z, x_vector, y_vector, z_vector, Nx, Ny, Nz) &
        !$omp& num_threads(threads)

        !$omp do
        do k = 1, Nz
            do i = 1, Ny
                X(i, :, k) = x_vector(:)
            end do
        end do
        !$omp end do

        !$omp do
        do k = 1, Nz
            do j = 1, Nx
                Y(:, j, k) = y_vector(:)
            end do
        end do
        !$omp end do

        !$omp do
        do j = 1, Nx
            do i = 1, Ny
                Z(i, j, :) = z_vector(:)
            end do
        end do
        !$omp end do

        !$omp end parallel

    end subroutine compute_openmp_cmplx_3D_dp

    !=================================================================================
    ! Compute the meshgrid of three complex vectors in 128-bit (Quadruple precision)
    !=================================================================================

    module subroutine compute_meshgrid_cmplx_3D_qp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, &
        loop_method)
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        complex(qp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        integer(isp), intent(in) :: Nx, Ny, Nz
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx, Nz), Y(Ny, Nx, Nz), Z(Ny, Nx, Nz))

        if (loop_method%use_do_classic) then
            call compute_do_classic_cmplx_3D_qp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_cmplx_3D_qp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_cmplx_3D_qp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_cmplx_3D_qp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_cmplx_3D_qp

    pure subroutine compute_do_classic_cmplx_3D_qp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        complex(qp), dimension(:, :, :), intent(out) :: X, Y, Z
        integer(isp) :: i, j, k

        do k = 1, Nz
            do i = 1, Ny
                X(i, :, k) = x_vector(:)
            end do
        end do

        do k = 1, Nz
            do j = 1, Nx
                Y(:, j, k) = y_vector(:)
            end do
        end do

        do j = 1, Nx
            do i = 1, Ny
                Z(i, j, :) = z_vector(:)
            end do
        end do

    end subroutine compute_do_classic_cmplx_3D_qp

    pure subroutine compute_do_vectorized_cmplx_3D_qp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        complex(qp), dimension(:, :, :), intent(out) :: X, Y, Z

        X = spread(spread(x_vector, 1, Ny), 3, Nz)
        Y = spread(spread(y_vector, 2, Nx), 3, Nz)
        Z = spread(spread(z_vector, 1, Nx), 1, Ny)

    end subroutine compute_do_vectorized_cmplx_3D_qp

    pure subroutine compute_do_concurrent_cmplx_3D_qp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        complex(qp), dimension(:, :, :), intent(out) :: X, Y, Z
        integer(isp) :: i, j, k

        do concurrent(k=1:Nz, i=1:Ny)
            X(i, :, k) = x_vector(:)
        end do

        do concurrent(k=1:Nz, j=1:Nx)
            Y(:, j, k) = y_vector(:)
        end do

        do concurrent(j=1:Nx, i=1:Ny)
            Z(i, j, :) = z_vector(:)
        end do

    end subroutine compute_do_concurrent_cmplx_3D_qp

    subroutine compute_openmp_cmplx_3D_qp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, threads)
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz, threads
        complex(qp), dimension(:, :, :), intent(out) :: X, Y, Z
        integer(isp) :: i, j, k

        !$omp parallel default(none) private(i, j, k) &
        !$omp& shared(X, Y, Z, x_vector, y_vector, z_vector, Nx, Ny, Nz) &
        !$omp& num_threads(threads)

        !$omp do
        do k = 1, Nz
            do i = 1, Ny
                X(i, :, k) = x_vector(:)
            end do
        end do
        !$omp end do

        !$omp do
        do k = 1, Nz
            do j = 1, Nx
                Y(:, j, k) = y_vector(:)
            end do
        end do
        !$omp end do

        !$omp do
        do j = 1, Nx
            do i = 1, Ny
                Z(i, j, :) = z_vector(:)
            end do
        end do
        !$omp end do

        !$omp end parallel

    end subroutine compute_openmp_cmplx_3D_qp

end submodule NAFPack_meshgrid_complex_3D_compute
