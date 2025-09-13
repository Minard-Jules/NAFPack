module test_meshgrid_integer_2D

    use NAFPack_kinds, only: i8, i16, isp, idp
    use NAFPack_constant, only: TOL_TEST_sp, TOL_TEST_dp, TOL_TEST_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_meshgrid, only: meshgrid, INDEXING_IJ, INDEXING_XY
    use NAFPack_loop_method, only: LoopMethod, init_loop_method

    implicit none(type, external)

    private
    public :: collect_meshgrid

contains

    subroutine collect_meshgrid(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("meshgrid integer i8 2D", test_meshgrid_i8_2D), &
                    new_unittest("meshgrid integer i16 2D", test_meshgrid_i16_2D), &
                    new_unittest("meshgrid integer isp 2D", test_meshgrid_isp_2D), &
                    new_unittest("meshgrid integer idp 2D", test_meshgrid_idp_2D) &
                    ]

    end subroutine collect_meshgrid

    subroutine test_meshgrid_i8_2D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 5, ny = 4
        integer(i8), dimension(:), allocatable :: x_vector, y_vector
        integer(i8), dimension(:, :), allocatable :: X, Y
        type(LoopMethod) :: loop_method
        integer :: i

        allocate (x_vector(nx), y_vector(ny))

        x_vector = [(int(i, i8), i=1, nx)]
        y_vector = [(int(i, i8), i=1, ny)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_i8_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_i8_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_i8_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_i8_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_i8_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_i8_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_i8_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_i8_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        deallocate (x_vector, y_vector)

    end subroutine test_meshgrid_i8_2D

    subroutine check_meshgrid_i8_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector
        integer(i8), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer :: i, j

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
    end subroutine check_meshgrid_i8_2D_ij

    subroutine check_meshgrid_i8_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        integer(i8), dimension(:), intent(in) :: x_vector, y_vector
        integer(i8), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer :: i, j

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
    end subroutine check_meshgrid_i8_2D_xy

    subroutine test_meshgrid_i16_2D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 5, ny = 4
        integer(i16), dimension(:), allocatable :: x_vector, y_vector
        integer(i16), dimension(:, :), allocatable :: X, Y
        type(LoopMethod) :: loop_method
        integer :: i

        allocate (x_vector(nx), y_vector(ny))

        x_vector = [(int(i, i16), i=1, nx)]
        y_vector = [(int(i, i16), i=1, ny)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_i16_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_i16_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_i16_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_i16_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_i16_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_i16_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_i16_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_i16_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        deallocate (x_vector, y_vector)

    end subroutine test_meshgrid_i16_2D

    subroutine check_meshgrid_i16_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector
        integer(i16), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer :: i, j

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
    end subroutine check_meshgrid_i16_2D_ij

    subroutine check_meshgrid_i16_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        integer(i16), dimension(:), intent(in) :: x_vector, y_vector
        integer(i16), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer :: i, j

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
    end subroutine check_meshgrid_i16_2D_xy

    subroutine test_meshgrid_isp_2D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 5, ny = 4
        integer(isp), dimension(:), allocatable :: x_vector, y_vector
        integer(isp), dimension(:, :), allocatable :: X, Y
        type(LoopMethod) :: loop_method
        integer :: i

        allocate (x_vector(nx), y_vector(ny))

        x_vector = [(int(i, isp), i=1, nx)]
        y_vector = [(int(i, isp), i=1, ny)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_isp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_isp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_isp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_isp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_isp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_isp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_isp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_isp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        deallocate (x_vector, y_vector)

    end subroutine test_meshgrid_isp_2D

    subroutine check_meshgrid_isp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer :: i, j

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
    end subroutine check_meshgrid_isp_2D_ij

    subroutine check_meshgrid_isp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        integer(isp), dimension(:), intent(in) :: x_vector, y_vector
        integer(isp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer :: i, j

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
    end subroutine check_meshgrid_isp_2D_xy

    subroutine test_meshgrid_idp_2D(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 5, ny = 4
        integer(idp), dimension(:), allocatable :: x_vector, y_vector
        integer(idp), dimension(:, :), allocatable :: X, Y
        type(LoopMethod) :: loop_method
        integer :: i

        allocate (x_vector(nx), y_vector(ny))

        x_vector = [(int(i, idp), i=1, nx)]
        y_vector = [(int(i, idp), i=1, ny)]

        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_idp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_idp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_idp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_IJ, loop_method=loop_method)
        call check_meshgrid_idp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)


        loop_method = init_loop_method(use_do_classic=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_idp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_vectorized=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_idp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_idp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        call meshgrid(x_vector, y_vector, X, Y, indexing=INDEXING_XY, loop_method=loop_method)
        call check_meshgrid_idp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        deallocate (X, Y)

        deallocate (x_vector, y_vector)

    end subroutine test_meshgrid_idp_2D

    subroutine check_meshgrid_idp_2D_ij(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector
        integer(idp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer :: i, j

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
    end subroutine check_meshgrid_idp_2D_ij

    subroutine check_meshgrid_idp_2D_xy(error, x_vector, y_vector, X, Y, nx, ny)
        type(error_type), allocatable, intent(out) :: error
        integer(idp), dimension(:), intent(in) :: x_vector, y_vector
        integer(idp), dimension(:, :), intent(in) :: X, Y
        integer, intent(in) :: nx, ny
        integer :: i, j

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
    end subroutine check_meshgrid_idp_2D_xy

end module test_meshgrid_integer_2D

program test
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_meshgrid_integer_2D

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("meshgrid integer 2D", collect_meshgrid) &
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
