submodule(NAFPack_meshgrid) NAFPack_meshgrid_real_3D

    implicit none(type, external)

contains

    module subroutine meshgrid_sp_3D( &
        x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
        real(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        real(sp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: sX, sY, sZ, i, j, k

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        sX = size(x_vector)
        sY = size(y_vector)
        sZ = size(z_vector)
        if (use_ij_indexing) then

            allocate (X(sY, sX, sZ), Y(sY, sX, sZ), Z(sY, sX, sZ))

            if (loop_method_used%use_do_classic) then

                do k = 1, sZ
                    do i = 1, sY
                        X(i, :, k) = x_vector
                    end do
                end do

                do k = 1, sZ
                    do j = 1, sX
                        Y(:, j, k) = y_vector
                    end do
                end do

                do j = 1, sX
                    do i = 1, sY
                        Z(i, j, :) = z_vector
                    end do
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(spread(x_vector, 1, sY), 3, sZ)
                Y = spread(spread(y_vector, 2, sX), 3, sZ)
                Z = spread(spread(z_vector, 1, sX), 1, sY)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(k=1:sZ, i=1:sY)
                    X(i, :, k) = x_vector
                end do

                do concurrent(k=1:sZ, j=1:sX)
                    Y(:, j, k) = y_vector
                end do

                do concurrent(j=1:sX, i=1:sY)
                    Z(i, j, :) = z_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j, k) &
                !$omp& shared(X, Y, Z, x_vector, y_vector, z_vector, sX, sY, sZ) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do k = 1, sZ
                    do i = 1, sY
                        X(i, :, k) = x_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do k = 1, sZ
                    do j = 1, sX
                        Y(:, j, k) = y_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do j = 1, sX
                    do i = 1, sY
                        Z(i, j, :) = z_vector
                    end do
                end do
                !$omp end do

                !$omp end parallel
            end if

        else if (use_xy_indexing) then

            allocate (X(sX, sY, sZ), Y(sX, sY, sZ), Z(sX, sY, sZ))

            if (loop_method_used%use_do_classic) then

                do k = 1, sZ
                    do i = 1, sY
                        X(:, i, k) = x_vector
                    end do
                end do

                do k = 1, sZ
                    do j = 1, sX
                        Y(j, :, k) = y_vector
                    end do
                end do

                do i = 1, sY
                    do j = 1, sX
                        Z(j, i, :) = z_vector
                    end do
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(spread(x_vector, 2, sY), 3, sZ)
                Y = spread(spread(y_vector, 1, sX), 3, sZ)
                Z = spread(spread(z_vector, 1, sX), 2, sY)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(k=1:sZ, i=1:sY)
                    X(:, i, k) = x_vector
                end do

                do concurrent(k=1:sZ, j=1:sX)
                    Y(j, :, k) = y_vector
                end do

                do concurrent(j=1:sX, i=1:sY)
                    Z(j, i, :) = z_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j, k) &
                !$omp& shared(X, Y, Z, x_vector, y_vector, z_vector, sX, sY, sZ) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do k = 1, sZ
                    do i = 1, sY
                        X(:, i, k) = x_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do k = 1, sZ
                    do j = 1, sX
                        Y(j, :, k) = y_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do i = 1, sY
                    do j = 1, sX
                        Z(j, i, :) = z_vector
                    end do
                end do
                !$omp end do

                !$omp end parallel
            end if
        end if

    end subroutine meshgrid_sp_3D

    module subroutine meshgrid_dp_3D( &
        x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
        real(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        real(dp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: sX, sY, sZ, i, j, k

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        sX = size(x_vector)
        sY = size(y_vector)
        sZ = size(z_vector)
        if (use_ij_indexing) then

            allocate (X(sY, sX, sZ), Y(sY, sX, sZ), Z(sY, sX, sZ))

            if (loop_method_used%use_do_classic) then

                do k = 1, sZ
                    do i = 1, sY
                        X(i, :, k) = x_vector
                    end do
                end do

                do k = 1, sZ
                    do j = 1, sX
                        Y(:, j, k) = y_vector
                    end do
                end do

                do j = 1, sX
                    do i = 1, sY
                        Z(i, j, :) = z_vector
                    end do
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(spread(x_vector, 1, sY), 3, sZ)
                Y = spread(spread(y_vector, 2, sX), 3, sZ)
                Z = spread(spread(z_vector, 1, sX), 1, sY)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(k=1:sZ, i=1:sY)
                    X(i, :, k) = x_vector
                end do

                do concurrent(k=1:sZ, j=1:sX)
                    Y(:, j, k) = y_vector
                end do

                do concurrent(j=1:sX, i=1:sY)
                    Z(i, j, :) = z_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j, k) &
                !$omp& shared(X, Y, Z, x_vector, y_vector, z_vector, sX, sY, sZ) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do k = 1, sZ
                    do i = 1, sY
                        X(i, :, k) = x_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do k = 1, sZ
                    do j = 1, sX
                        Y(:, j, k) = y_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do j = 1, sX
                    do i = 1, sY
                        Z(i, j, :) = z_vector
                    end do
                end do
                !$omp end do

                !$omp end parallel
            end if

        else if (use_xy_indexing) then

            allocate (X(sX, sY, sZ), Y(sX, sY, sZ), Z(sX, sY, sZ))

            if (loop_method_used%use_do_classic) then

                do k = 1, sZ
                    do i = 1, sY
                        X(:, i, k) = x_vector
                    end do
                end do

                do k = 1, sZ
                    do j = 1, sX
                        Y(j, :, k) = y_vector
                    end do
                end do

                do i = 1, sY
                    do j = 1, sX
                        Z(j, i, :) = z_vector
                    end do
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(spread(x_vector, 2, sY), 3, sZ)
                Y = spread(spread(y_vector, 1, sX), 3, sZ)
                Z = spread(spread(z_vector, 1, sX), 2, sY)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(k=1:sZ, i=1:sY)
                    X(:, i, k) = x_vector
                end do

                do concurrent(k=1:sZ, j=1:sX)
                    Y(j, :, k) = y_vector
                end do

                do concurrent(j=1:sX, i=1:sY)
                    Z(j, i, :) = z_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j, k) &
                !$omp& shared(X, Y, Z, x_vector, y_vector, z_vector, sX, sY, sZ) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do k = 1, sZ
                    do i = 1, sY
                        X(:, i, k) = x_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do k = 1, sZ
                    do j = 1, sX
                        Y(j, :, k) = y_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do i = 1, sY
                    do j = 1, sX
                        Z(j, i, :) = z_vector
                    end do
                end do
                !$omp end do

                !$omp end parallel
            end if
        end if

    end subroutine meshgrid_dp_3D

    module subroutine meshgrid_qp_3D( &
        x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
        real(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        real(qp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: sX, sY, sZ, i, j, k

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        sX = size(x_vector)
        sY = size(y_vector)
        sZ = size(z_vector)
        if (use_ij_indexing) then

            allocate (X(sY, sX, sZ), Y(sY, sX, sZ), Z(sY, sX, sZ))

            if (loop_method_used%use_do_classic) then

                do k = 1, sZ
                    do i = 1, sY
                        X(i, :, k) = x_vector
                    end do
                end do

                do k = 1, sZ
                    do j = 1, sX
                        Y(:, j, k) = y_vector
                    end do
                end do

                do j = 1, sX
                    do i = 1, sY
                        Z(i, j, :) = z_vector
                    end do
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(spread(x_vector, 1, sY), 3, sZ)
                Y = spread(spread(y_vector, 2, sX), 3, sZ)
                Z = spread(spread(z_vector, 1, sX), 1, sY)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(k=1:sZ, i=1:sY)
                    X(i, :, k) = x_vector
                end do

                do concurrent(k=1:sZ, j=1:sX)
                    Y(:, j, k) = y_vector
                end do

                do concurrent(j=1:sX, i=1:sY)
                    Z(i, j, :) = z_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j, k) &
                !$omp& shared(X, Y, Z, x_vector, y_vector, z_vector, sX, sY, sZ) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do k = 1, sZ
                    do i = 1, sY
                        X(i, :, k) = x_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do k = 1, sZ
                    do j = 1, sX
                        Y(:, j, k) = y_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do j = 1, sX
                    do i = 1, sY
                        Z(i, j, :) = z_vector
                    end do
                end do
                !$omp end do

                !$omp end parallel
            end if

        else if (use_xy_indexing) then

            allocate (X(sX, sY, sZ), Y(sX, sY, sZ), Z(sX, sY, sZ))

            if (loop_method_used%use_do_classic) then

                do k = 1, sZ
                    do i = 1, sY
                        X(:, i, k) = x_vector
                    end do
                end do

                do k = 1, sZ
                    do j = 1, sX
                        Y(j, :, k) = y_vector
                    end do
                end do

                do i = 1, sY
                    do j = 1, sX
                        Z(j, i, :) = z_vector
                    end do
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(spread(x_vector, 2, sY), 3, sZ)
                Y = spread(spread(y_vector, 1, sX), 3, sZ)
                Z = spread(spread(z_vector, 1, sX), 2, sY)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(k=1:sZ, i=1:sY)
                    X(:, i, k) = x_vector
                end do

                do concurrent(k=1:sZ, j=1:sX)
                    Y(j, :, k) = y_vector
                end do

                do concurrent(j=1:sX, i=1:sY)
                    Z(j, i, :) = z_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j, k) &
                !$omp& shared(X, Y, Z, x_vector, y_vector, z_vector, sX, sY, sZ) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do k = 1, sZ
                    do i = 1, sY
                        X(:, i, k) = x_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do k = 1, sZ
                    do j = 1, sX
                        Y(j, :, k) = y_vector
                    end do
                end do
                !$omp end do

                !$omp do
                do i = 1, sY
                    do j = 1, sX
                        Z(j, i, :) = z_vector
                    end do
                end do
                !$omp end do

                !$omp end parallel
            end if
        end if

    end subroutine meshgrid_qp_3D

end submodule NAFPack_meshgrid_real_3D
