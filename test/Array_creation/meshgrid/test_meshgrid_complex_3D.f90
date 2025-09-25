submodule (test_meshgrid) test_meshgrid_complex_3D 

    implicit none(type, external)

contains

    module subroutine test_meshgrid_cmplx_sp_3D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 3, ny = 2, nz = 4
        complex(sp), dimension(:), allocatable :: x_vector, y_vector, z_vector
        complex(sp), dimension(:, :, :), allocatable :: X, Y, Z
        type(LoopMethod) :: loop_method
        integer(isp) :: i

        allocate(x_vector(nx), y_vector(ny))

        x_vector = [(cmplx(i, -i, kind=sp), i=1, nx)]
        y_vector = [(cmplx(i, -i, kind=sp), i=1, ny)]
        z_vector = [(cmplx(i, -i, kind=sp), i=1, nz)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_sp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        deallocate(x_vector, y_vector, z_vector)

    end subroutine test_meshgrid_cmplx_sp_3D

    subroutine check_meshgrid_cmplx_sp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        type(error_type), allocatable, intent(out) :: error
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        complex(sp), dimension(:, :, :), intent(in) :: X, Y, Z
        integer, intent(in) :: nx, ny, nz
        integer(isp) :: i, j, k
        
        call check(error, size(X, 1) == ny .and. size(X, 2) == nx .and. size(X, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Y, 1) == ny .and. size(Y, 2) == nx .and. size(Y, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Z, 1) == ny .and. size(Z, 2) == nx .and. size(Z, 3) == nz)
        if (allocated(error)) return

        do k = 1, nz
            do i = 1, ny
                call check(error, all(abs(X(i, :, k) - x_vector) < TOL_TEST_sp))
                if (allocated(error)) return
            end do
        end do

        do k = 1, nz
            do j = 1, nx
                call check(error, all(abs(Y(:, j, k) - y_vector) < TOL_TEST_sp))
                if (allocated(error)) return
            end do
        end do

        do j = 1, nx
            do i = 1, ny
                call check(error, all(abs(Z(i, j, :) - z_vector) < TOL_TEST_sp))
                if (allocated(error)) return
            end do
        end do
    end subroutine check_meshgrid_cmplx_sp_3D_ij

    subroutine check_meshgrid_cmplx_sp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        type(error_type), allocatable, intent(out) :: error
        complex(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        complex(sp), dimension(:, :, :), intent(in) :: X, Y, Z
        integer, intent(in) :: nx, ny, nz
        integer(isp) :: i, j, k

        call check(error, size(X, 1) == nx .and. size(X, 2) == ny .and. size(X, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Y, 1) == nx .and. size(Y, 2) == ny .and. size(Y, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Z, 1) == nx .and. size(Z, 2) == ny .and. size(Z, 3) == nz)
        if (allocated(error)) return

        do k = 1, nz
            do j = 1, ny
                call check(error, all(abs(X(:, j, k) - x_vector) < TOL_TEST_sp))
                if (allocated(error)) return
            end do
        end do

        do k = 1, nz
            do i = 1, nx
                call check(error, all(abs(Y(i, :, k) - y_vector) < TOL_TEST_sp))
                if (allocated(error)) return
            end do
        end do

        do j = 1, ny
            do i = 1, nx
                call check(error, all(abs(Z(i, j, :) - z_vector) < TOL_TEST_sp))
                if (allocated(error)) return
            end do
        end do
    end subroutine check_meshgrid_cmplx_sp_3D_xy

    module subroutine test_meshgrid_cmplx_dp_3D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 3, ny = 2, nz = 4
        complex(dp), dimension(:), allocatable :: x_vector, y_vector, z_vector
        complex(dp), dimension(:, :, :), allocatable :: X, Y, Z
        type(LoopMethod) :: loop_method
        integer(isp) :: i

        allocate(x_vector(nx), y_vector(ny), z_vector(nz))

        x_vector = [(cmplx(i, -i, kind=dp), i=1, nx)]
        y_vector = [(cmplx(i, -i, kind=dp), i=1, ny)]
        z_vector = [(cmplx(i, -i, kind=dp), i=1, nz)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_dp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        deallocate(x_vector, y_vector, z_vector)

    end subroutine test_meshgrid_cmplx_dp_3D

    subroutine check_meshgrid_cmplx_dp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        type(error_type), allocatable, intent(out) :: error
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        complex(dp), dimension(:, :, :), intent(in) :: X, Y, Z
        integer, intent(in) :: nx, ny, nz
        integer(isp) :: i, j, k

        call check(error, size(X, 1) == ny .and. size(X, 2) == nx .and. size(X, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Y, 1) == ny .and. size(Y, 2) == nx .and. size(Y, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Z, 1) == ny .and. size(Z, 2) == nx .and. size(Z, 3) == nz)
        if (allocated(error)) return

        do k = 1, nz
            do i = 1, ny
                call check(error, all(abs(X(i, :, k) - x_vector) < TOL_TEST_dp))
                if (allocated(error)) return
            end do
        end do

        do k = 1, nz
            do j = 1, nx
                call check(error, all(abs(Y(:, j, k) - y_vector) < TOL_TEST_dp))
                if (allocated(error)) return
            end do
        end do

        do j = 1, nx
            do i = 1, ny
                call check(error, all(abs(Z(i, j, :) - z_vector) < TOL_TEST_dp))
                if (allocated(error)) return
            end do
        end do
    end subroutine check_meshgrid_cmplx_dp_3D_ij

    subroutine check_meshgrid_cmplx_dp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        type(error_type), allocatable, intent(out) :: error
        complex(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        complex(dp), dimension(:, :, :), intent(in) :: X, Y, Z
        integer, intent(in) :: nx, ny, nz
        integer(isp) :: i, j, k

        call check(error, size(X, 1) == nx .and. size(X, 2) == ny .and. size(X, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Y, 1) == nx .and. size(Y, 2) == ny .and. size(Y, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Z, 1) == nx .and. size(Z, 2) == ny .and. size(Z, 3) == nz)
        if (allocated(error)) return

        do k = 1, nz
            do j = 1, ny
                call check(error, all(abs(X(:, j, k) - x_vector) < TOL_TEST_dp))
                if (allocated(error)) return
            end do
        end do

        do k = 1, nz
            do i = 1, nx
                call check(error, all(abs(Y(i, :, k) - y_vector) < TOL_TEST_dp))
                if (allocated(error)) return
            end do
        end do

        do j = 1, ny
            do i = 1, nx
                call check(error, all(abs(Z(i, j, :) - z_vector) < TOL_TEST_dp))
                if (allocated(error)) return
            end do
        end do
    end subroutine check_meshgrid_cmplx_dp_3D_xy

    module subroutine test_meshgrid_cmplx_qp_3D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 3, ny = 2, nz = 4
        complex(qp), dimension(:), allocatable :: x_vector, y_vector, z_vector
        complex(qp), dimension(:, :, :), allocatable :: X, Y, Z
        type(LoopMethod) :: loop_method
        integer(isp) :: i

        allocate(x_vector(nx), y_vector(ny), z_vector(nz))

        x_vector = [(cmplx(i, -i, kind=qp), i=1, nx)]
        y_vector = [(cmplx(i, -i, kind=qp), i=1, ny)]
        z_vector = [(cmplx(i, -i, kind=qp), i=1, nz)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, z_vector, X, Y, Z, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_cmplx_qp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        deallocate(X, Y, Z)

        deallocate(x_vector, y_vector, z_vector)

    end subroutine test_meshgrid_cmplx_qp_3D

    subroutine check_meshgrid_cmplx_qp_3D_ij(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        type(error_type), allocatable, intent(out) :: error
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        complex(qp), dimension(:, :, :), intent(in) :: X, Y, Z
        integer, intent(in) :: nx, ny, nz
        integer(isp) :: i, j, k

        call check(error, size(X, 1) == ny .and. size(X, 2) == nx .and. size(X, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Y, 1) == ny .and. size(Y, 2) == nx .and. size(Y, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Z, 1) == ny .and. size(Z, 2) == nx .and. size(Z, 3) == nz)
        if (allocated(error)) return

        do k = 1, nz
            do i = 1, ny
                call check(error, all(abs(X(i, :, k) - x_vector) < TOL_TEST_qp))
                if (allocated(error)) return
            end do
        end do

        do k = 1, nz
            do j = 1, nx
                call check(error, all(abs(Y(:, j, k) - y_vector) < TOL_TEST_qp))
                if (allocated(error)) return
            end do
        end do

        do j = 1, nx
            do i = 1, ny
                call check(error, all(abs(Z(i, j, :) - z_vector) < TOL_TEST_qp))
                if (allocated(error)) return
            end do
        end do
    end subroutine check_meshgrid_cmplx_qp_3D_ij

    subroutine check_meshgrid_cmplx_qp_3D_xy(error, x_vector, y_vector, z_vector, X, Y, Z, nx, ny, nz)
        type(error_type), allocatable, intent(out) :: error
        complex(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
        complex(qp), dimension(:, :, :), intent(in) :: X, Y, Z
        integer, intent(in) :: nx, ny, nz
        integer(isp) :: i, j, k

        call check(error, size(X, 1) == nx .and. size(X, 2) == ny .and. size(X, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Y, 1) == nx .and. size(Y, 2) == ny .and. size(Y, 3) == nz)
        if (allocated(error)) return
        call check(error, size(Z, 1) == nx .and. size(Z, 2) == ny .and. size(Z, 3) == nz)
        if (allocated(error)) return

        do k = 1, nz
            do j = 1, ny
                call check(error, all(abs(X(:, j, k) - x_vector) < TOL_TEST_qp))
                if (allocated(error)) return
            end do
        end do

        do k = 1, nz
            do i = 1, nx
                call check(error, all(abs(Y(i, :, k) - y_vector) < TOL_TEST_qp))
                if (allocated(error)) return
            end do
        end do

        do j = 1, ny
            do i = 1, nx
                call check(error, all(abs(Z(i, j, :) - z_vector) < TOL_TEST_qp))
                if (allocated(error)) return
            end do
        end do
    end subroutine check_meshgrid_cmplx_qp_3D_xy

end submodule test_meshgrid_complex_3D