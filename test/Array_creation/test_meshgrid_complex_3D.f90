module test_meshgrid_complex_3D

    use NAFPack_kinds, only: dp, sp, qp
    use NAFPack_constant, only: TOL_TEST_sp, TOL_TEST_dp, TOL_TEST_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_loop_method, only: LoopMethod, init_loop_method
    use NAFPack_meshgrid, only: meshgrid, INDEXING_IJ, INDEXING_XY

    implicit none(type, external)

    private
    public :: collect_meshgrid

contains

    subroutine collect_meshgrid(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("meshgrid complex sp 3D", test_meshgrid_cmplx_sp_3D), &
                    new_unittest("meshgrid complex dp 3D", test_meshgrid_cmplx_dp_3D), &
                    new_unittest("meshgrid complex qp 3D", test_meshgrid_cmplx_qp_3D) &
                    ]

    end subroutine collect_meshgrid

    subroutine test_meshgrid_cmplx_sp_3D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 3, ny = 2, nz = 4
        complex(sp), dimension(:), allocatable :: x_vector, y_vector, z_vector
        complex(sp), dimension(:, :, :), allocatable :: X, Y, Z
        type(LoopMethod) :: loop_method
        integer :: i

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
        integer :: i, j, k
        
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
        integer :: i, j, k

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

    subroutine test_meshgrid_cmplx_dp_3D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 3, ny = 2, nz = 4
        complex(dp), dimension(:), allocatable :: x_vector, y_vector, z_vector
        complex(dp), dimension(:, :, :), allocatable :: X, Y, Z
        type(LoopMethod) :: loop_method
        integer :: i

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
        integer :: i, j, k

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
        integer :: i, j, k

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

    subroutine test_meshgrid_cmplx_qp_3D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 3, ny = 2, nz = 4
        complex(qp), dimension(:), allocatable :: x_vector, y_vector, z_vector
        complex(qp), dimension(:, :, :), allocatable :: X, Y, Z
        type(LoopMethod) :: loop_method
        integer :: i

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
        integer :: i, j, k

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
        integer :: i, j, k

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

end module test_meshgrid_complex_3D

program test
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_meshgrid_complex_3D

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("meshgrid complex 3D", collect_meshgrid) &
                 ]

    call init_color_output(.true.)

    do is = 1, size(testsuites)
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

end program test
