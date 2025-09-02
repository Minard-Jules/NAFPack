submodule(NAFPack_meshgrid) NAFPack_meshgrid_complex_2D

    implicit none(type, external)

contains

    module subroutine meshgrid_cmplx_sp_2D( &
        x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector
        complex(sp), dimension(:, :), allocatable, intent(out) :: X, Y
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: sX, sY, i, j

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        sX = size(x_vector)
        sY = size(y_vector)
        if (use_ij_indexing) then

            allocate (X(sY, sX), Y(sY, sX))

            if (loop_method_used%use_do_classic) then
                do i = 1, sY
                    X(i, :) = x_vector
                end do

                do j = 1, sX
                    Y(:, j) = y_vector
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(x_vector, 1, sY)
                Y = spread(y_vector, 2, sX)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(i=1:sY)
                    X(i, :) = x_vector
                end do

                do concurrent(j=1:sX)
                    Y(:, j) = y_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j) &
                !$omp& shared(X, Y, x_vector, y_vector, sX, sY) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do i = 1, sY
                    X(i, :) = x_vector
                end do
                !$omp end do

                !$omp do
                do j = 1, sX
                    Y(:, j) = y_vector
                end do
                !$omp end do

                !$omp end parallel
            end if
        else if (use_xy_indexing) then
            allocate (X(sX, sY), Y(sX, sY))

            if (loop_method_used%use_do_classic) then
                do i = 1, sY
                    X(:, i) = x_vector
                end do

                do j = 1, sX
                    Y(j, :) = y_vector
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(x_vector, 2, sY)
                Y = spread(y_vector, 1, sX)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(i=1:sY)
                    X(:, i) = x_vector
                end do

                do concurrent(j=1:sX)
                    Y(j, :) = y_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j) &
                !$omp& shared(X, Y, x_vector, y_vector, sX, sY) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do i = 1, sY
                    X(:, i) = x_vector
                end do
                !$omp end do

                !$omp do
                do j = 1, sX
                    Y(j, :) = y_vector
                end do
                !$omp end do

                !$omp end parallel
            end if
        end if

    end subroutine meshgrid_cmplx_sp_2D

    module subroutine meshgrid_cmplx_dp_2D( &
        x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector
        complex(dp), dimension(:, :), allocatable, intent(out) :: X, Y
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: sX, sY, i, j

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        sX = size(x_vector)
        sY = size(y_vector)
        if (use_ij_indexing) then

            allocate (X(sY, sX), Y(sY, sX))

            if (loop_method_used%use_do_classic) then
                do i = 1, sY
                    X(i, :) = x_vector
                end do

                do j = 1, sX
                    Y(:, j) = y_vector
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(x_vector, 1, sY)
                Y = spread(y_vector, 2, sX)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(i=1:sY)
                    X(i, :) = x_vector
                end do

                do concurrent(j=1:sX)
                    Y(:, j) = y_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j) &
                !$omp& shared(X, Y, x_vector, y_vector, sX, sY) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do i = 1, sY
                    X(i, :) = x_vector
                end do
                !$omp end do

                !$omp do
                do j = 1, sX
                    Y(:, j) = y_vector
                end do
                !$omp end do

                !$omp end parallel
            end if
        else if (use_xy_indexing) then
            allocate (X(sX, sY), Y(sX, sY))

            if (loop_method_used%use_do_classic) then
                do i = 1, sY
                    X(:, i) = x_vector
                end do

                do j = 1, sX
                    Y(j, :) = y_vector
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(x_vector, 2, sY)
                Y = spread(y_vector, 1, sX)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(i=1:sY)
                    X(:, i) = x_vector
                end do

                do concurrent(j=1:sX)
                    Y(j, :) = y_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j) &
                !$omp& shared(X, Y, x_vector, y_vector, sX, sY) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do i = 1, sY
                    X(:, i) = x_vector
                end do
                !$omp end do

                !$omp do
                do j = 1, sX
                    Y(j, :) = y_vector
                end do
                !$omp end do

                !$omp end parallel
            end if
        end if

    end subroutine meshgrid_cmplx_dp_2D

    module subroutine meshgrid_cmplx_qp_2D( &
        x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector
        complex(qp), dimension(:, :), allocatable, intent(out) :: X, Y
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        type(LoopMethod), optional, intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        logical :: use_ij_indexing, use_xy_indexing
        integer :: sX, sY, i, j

        call check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        sX = size(x_vector)
        sY = size(y_vector)
        if (use_ij_indexing) then

            allocate (X(sY, sX), Y(sY, sX))

            if (loop_method_used%use_do_classic) then
                do i = 1, sY
                    X(i, :) = x_vector
                end do

                do j = 1, sX
                    Y(:, j) = y_vector
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(x_vector, 1, sY)
                Y = spread(y_vector, 2, sX)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(i=1:sY)
                    X(i, :) = x_vector
                end do

                do concurrent(j=1:sX)
                    Y(:, j) = y_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j) &
                !$omp& shared(X, Y, x_vector, y_vector, sX, sY) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do i = 1, sY
                    X(i, :) = x_vector
                end do
                !$omp end do

                !$omp do
                do j = 1, sX
                    Y(:, j) = y_vector
                end do
                !$omp end do

                !$omp end parallel
            end if
        else if (use_xy_indexing) then
            allocate (X(sX, sY), Y(sX, sY))

            if (loop_method_used%use_do_classic) then
                do i = 1, sY
                    X(:, i) = x_vector
                end do

                do j = 1, sX
                    Y(j, :) = y_vector
                end do
            else if (loop_method_used%use_vectorized) then
                X = spread(x_vector, 2, sY)
                Y = spread(y_vector, 1, sX)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(i=1:sY)
                    X(:, i) = x_vector
                end do

                do concurrent(j=1:sX)
                    Y(j, :) = y_vector
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel default(none) private(i, j) &
                !$omp& shared(X, Y, x_vector, y_vector, sX, sY) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)

                !$omp do
                do i = 1, sY
                    X(:, i) = x_vector
                end do
                !$omp end do

                !$omp do
                do j = 1, sX
                    Y(j, :) = y_vector
                end do
                !$omp end do

                !$omp end parallel
            end if
        end if

    end subroutine meshgrid_cmplx_qp_2D

end submodule NAFPack_meshgrid_complex_2D
