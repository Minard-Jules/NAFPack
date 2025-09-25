submodule (test_meshgrid) test_meshgrid_complex_2D

    implicit none(type, external)

contains

    module subroutine test_meshgrid_cmplx_sp_2D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 5, ny = 4
        complex(sp), dimension(:), allocatable :: x_vector, y_vector
        complex(sp), dimension(:, :), allocatable :: X, Y
        type(LoopMethod) :: loop_method
        integer(isp) :: i

        allocate (x_vector(nx), y_vector(ny))

        x_vector = [(cmplx(i, -i, kind=sp), i=1, nx)]
        y_vector = [(cmplx(i, -i, kind=sp), i=1, ny)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        deallocate (x_vector, y_vector)

    end subroutine test_meshgrid_cmplx_sp_2D

    subroutine check_meshgrid_cmplx_sp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector
        complex(sp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer(isp) :: i, j

        call check(error, size(X, 1) == ny .and. size(X, 2) == nx)
        if (allocated(error)) return
        call check(error, size(Y, 1) == ny .and. size(Y, 2) == nx)
        if (allocated(error)) return

        do i = 1, ny
            call check(error, all(abs(X(i, :) - x_vector) < TOL_TEST_sp))
            if (allocated(error)) return
        end do

        do j = 1, nx
            call check(error, all(abs(Y(:, j) - y_vector) < TOL_TEST_sp))
            if (allocated(error)) return
        end do
    end subroutine check_meshgrid_cmplx_sp_2D_ij

    subroutine check_meshgrid_cmplx_sp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector
        complex(sp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer(isp) :: i, j

        call check(error, size(X, 1) == nx .and. size(X, 2) == ny)
        if (allocated(error)) return
        call check(error, size(Y, 1) == nx .and. size(Y, 2) == ny)
        if (allocated(error)) return

        do j = 1, ny
            call check(error, all(abs(X(:, j) - x_vector) < TOL_TEST_sp))
            if (allocated(error)) return
        end do
        do i = 1, nx
            call check(error, all(abs(Y(i, :) - y_vector) < TOL_TEST_sp))
            if (allocated(error)) return
        end do
    end subroutine check_meshgrid_cmplx_sp_2D_xy

    module subroutine test_meshgrid_cmplx_dp_2D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 5, ny = 4
        complex(dp), dimension(:), allocatable :: x_vector, y_vector
        complex(dp), dimension(:, :), allocatable :: X, Y
        type(LoopMethod) :: loop_method
        integer(isp) :: i

        allocate (x_vector(nx), y_vector(ny))

        x_vector = [(cmplx(i, -i, kind=dp), i=1, nx)]
        y_vector = [(cmplx(i, -i, kind=dp), i=1, ny)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        deallocate (x_vector, y_vector)

    end subroutine test_meshgrid_cmplx_dp_2D

    subroutine check_meshgrid_cmplx_dp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector
        complex(dp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer(isp) :: i, j

        call check(error, size(X, 1) == ny .and. size(X, 2) == nx)
        if (allocated(error)) return
        call check(error, size(Y, 1) == ny .and. size(Y, 2) == nx)
        if (allocated(error)) return

        do i = 1, ny
            call check(error, all(abs(X(i, :) - x_vector) < TOL_TEST_dp))
            if (allocated(error)) return
        end do

        do j = 1, nx
            call check(error, all(abs(Y(:, j) - y_vector) < TOL_TEST_dp))
            if (allocated(error)) return
        end do
    end subroutine check_meshgrid_cmplx_dp_2D_ij

    subroutine check_meshgrid_cmplx_dp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector
        complex(dp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer(isp) :: i, j

        call check(error, size(X, 1) == nx .and. size(X, 2) == ny)
        if (allocated(error)) return
        call check(error, size(Y, 1) == nx .and. size(Y, 2) == ny)
        if (allocated(error)) return

        do j = 1, ny
            call check(error, all(abs(X(:, j) - x_vector) < TOL_TEST_dp))
            if (allocated(error)) return
        end do
        do i = 1, nx
            call check(error, all(abs(Y(i, :) - y_vector) < TOL_TEST_dp))
            if (allocated(error)) return
        end do
    end subroutine check_meshgrid_cmplx_dp_2D_xy

    module subroutine test_meshgrid_cmplx_qp_2D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 5, ny = 4
        complex(qp), dimension(:), allocatable :: x_vector, y_vector
        complex(qp), dimension(:, :), allocatable :: X, Y
        type(LoopMethod) :: loop_method
        integer(isp) :: i

        allocate (x_vector(nx), y_vector(ny))

        x_vector = [(cmplx(i, -i, kind=qp), i=1, nx)]
        y_vector = [(cmplx(i, -i, kind=qp), i=1, ny)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        deallocate (x_vector, y_vector)

    end subroutine test_meshgrid_cmplx_qp_2D

    subroutine check_meshgrid_cmplx_qp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector
        complex(qp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer(isp) :: i, j

        call check(error, size(X, 1) == ny .and. size(X, 2) == nx)
        if (allocated(error)) return
        call check(error, size(Y, 1) == ny .and. size(Y, 2) == nx)
        if (allocated(error)) return

        do i = 1, ny
            call check(error, all(abs(X(i, :) - x_vector) < TOL_TEST_qp))
            if (allocated(error)) return
        end do

        do j = 1, nx
            call check(error, all(abs(Y(:, j) - y_vector) < TOL_TEST_qp))
            if (allocated(error)) return
        end do
    end subroutine check_meshgrid_cmplx_qp_2D_ij

    subroutine check_meshgrid_cmplx_qp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector
        complex(qp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer(isp) :: i, j

        call check(error, size(X, 1) == nx .and. size(X, 2) == ny)
        if (allocated(error)) return
        call check(error, size(Y, 1) == nx .and. size(Y, 2) == ny)
        if (allocated(error)) return

        do j = 1, ny
            call check(error, all(abs(X(:, j) - x_vector) < TOL_TEST_qp))
            if (allocated(error)) return
        end do
        do i = 1, nx
            call check(error, all(abs(Y(i, :) - y_vector) < TOL_TEST_qp))
            if (allocated(error)) return
        end do
    end subroutine check_meshgrid_cmplx_qp_2D_xy

end submodule test_meshgrid_complex_2D