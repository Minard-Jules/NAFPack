submodule(NAFPack_meshgrid:NAFPack_meshgrid_integer_3D) NAFPack_meshgrid_integer_3D_compute

    implicit none(type, external)

contains

    !=================================================================================
    ! Compute the meshgrid of three integer vectors in 8-bit
    !=================================================================================

    module subroutine compute_meshgrid_integer_3D_i8( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, &
        loop_method)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(i8), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        integer(isp), intent(in) :: Nx, Ny, Nz
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx, Nz), Y(Ny, Nx, Nz), Z(Ny, Nx, Nz))

        if (loop_method%use_do_classic) then
            call compute_do_classic_integer_3D_i8(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_integer_3D_i8(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_integer_3D_i8(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_integer_3D_i8(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_integer_3D_i8

    pure subroutine compute_do_classic_integer_3D_i8( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(i8), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_do_classic_integer_3D_i8

    pure subroutine compute_do_vectorized_integer_3D_i8( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(i8), dimension(:, :, :), intent(out) :: X, Y, Z

        X = spread(spread(x_vector, 1, Ny), 3, Nz)
        Y = spread(spread(y_vector, 2, Nx), 3, Nz)
        Z = spread(spread(z_vector, 1, Nx), 1, Ny)

    end subroutine compute_do_vectorized_integer_3D_i8

    pure subroutine compute_do_concurrent_integer_3D_i8( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(i8), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_do_concurrent_integer_3D_i8

    subroutine compute_openmp_integer_3D_i8( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, threads)
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz, threads
        integer(i8), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_openmp_integer_3D_i8

    !=================================================================================
    ! Compute the meshgrid of three integer vectors in 16-bit
    !=================================================================================

    module subroutine compute_meshgrid_integer_3D_i16( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, &
        loop_method)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(i16), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        integer(isp), intent(in) :: Nx, Ny, Nz
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx, Nz), Y(Ny, Nx, Nz), Z(Ny, Nx, Nz))

        if (loop_method%use_do_classic) then
            call compute_do_classic_integer_3D_i16(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_integer_3D_i16(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_integer_3D_i16(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_integer_3D_i16(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_integer_3D_i16

    pure subroutine compute_do_classic_integer_3D_i16( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(i16), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_do_classic_integer_3D_i16

    pure subroutine compute_do_vectorized_integer_3D_i16( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(i16), dimension(:, :, :), intent(out) :: X, Y, Z

        X = spread(spread(x_vector, 1, Ny), 3, Nz)
        Y = spread(spread(y_vector, 2, Nx), 3, Nz)
        Z = spread(spread(z_vector, 1, Nx), 1, Ny)

    end subroutine compute_do_vectorized_integer_3D_i16

    pure subroutine compute_do_concurrent_integer_3D_i16( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(i16), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_do_concurrent_integer_3D_i16

    subroutine compute_openmp_integer_3D_i16( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, threads)
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz, threads
        integer(i16), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_openmp_integer_3D_i16

    !=================================================================================
    ! Compute the meshgrid of three integer vectors in 32-bit (single precision)
    !=================================================================================

    module subroutine compute_meshgrid_integer_3D_isp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, &
        loop_method)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        integer(isp), intent(in) :: Nx, Ny, Nz
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx, Nz), Y(Ny, Nx, Nz), Z(Ny, Nx, Nz))

        if (loop_method%use_do_classic) then
            call compute_do_classic_integer_3D_isp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_integer_3D_isp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_integer_3D_isp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_integer_3D_isp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_integer_3D_isp

    pure subroutine compute_do_classic_integer_3D_isp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(isp), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_do_classic_integer_3D_isp

    pure subroutine compute_do_vectorized_integer_3D_isp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(isp), dimension(:, :, :), intent(out) :: X, Y, Z

        X = spread(spread(x_vector, 1, Ny), 3, Nz)
        Y = spread(spread(y_vector, 2, Nx), 3, Nz)
        Z = spread(spread(z_vector, 1, Nx), 1, Ny)

    end subroutine compute_do_vectorized_integer_3D_isp

    pure subroutine compute_do_concurrent_integer_3D_isp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(isp), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_do_concurrent_integer_3D_isp

    subroutine compute_openmp_integer_3D_isp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, threads)
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz, threads
        integer(isp), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_openmp_integer_3D_isp

    !=================================================================================
    ! Compute the meshgrid of three integer vectors in 64-bit (double precision)
    !=================================================================================

    module subroutine compute_meshgrid_integer_3D_idp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, &
        loop_method)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(idp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        integer(isp), intent(in) :: Nx, Ny, Nz
        type(LoopMethod), intent(in) :: loop_method

        allocate (X(Ny, Nx, Nz), Y(Ny, Nx, Nz), Z(Ny, Nx, Nz))

        if (loop_method%use_do_classic) then
            call compute_do_classic_integer_3D_idp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_vectorized) then
            call compute_do_vectorized_integer_3D_idp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%use_do_concurrent) then
            call compute_do_concurrent_integer_3D_idp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz)
        else if (loop_method%parallel%use_openmp) then
            call compute_openmp_integer_3D_idp(x_vector, y_vector, z_vector, X, Y, Z, Nx, Ny, Nz, &
                                            loop_method%parallel%num_threads)
        end if

    end subroutine compute_meshgrid_integer_3D_idp

    pure subroutine compute_do_classic_integer_3D_idp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(idp), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_do_classic_integer_3D_idp

    pure subroutine compute_do_vectorized_integer_3D_idp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(idp), dimension(:, :, :), intent(out) :: X, Y, Z

        X = spread(spread(x_vector, 1, Ny), 3, Nz)
        Y = spread(spread(y_vector, 2, Nx), 3, Nz)
        Z = spread(spread(z_vector, 1, Nx), 1, Ny)

    end subroutine compute_do_vectorized_integer_3D_idp

    pure subroutine compute_do_concurrent_integer_3D_idp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz
        integer(idp), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_do_concurrent_integer_3D_idp

    subroutine compute_openmp_integer_3D_idp( &
        x_vector, y_vector, z_vector, &
        X, Y, Z, &
        Nx, Ny, Nz, threads)
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        integer(isp), intent(in) :: Nx, Ny, Nz, threads
        integer(idp), dimension(:, :, :), intent(out) :: X, Y, Z
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

    end subroutine compute_openmp_integer_3D_idp

end submodule NAFPack_meshgrid_integer_3D_compute
