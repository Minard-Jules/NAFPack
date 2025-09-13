module test_Fourier_Transform_fft_cmplx

    use NAFPack_kinds, only: dp, sp, qp
    use NAFPack_constant, only: TOL_TEST_sp, TOL_TEST_dp, TOL_TEST_qp
    use NAFPack_constant, only: pi_sp, im_sp, pi_dp, im_dp, pi_qp, im_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_loop_method, only: LoopMethod, init_loop_method
    use NAFPack_Fourier_Transform, only: Fourier_Transform

    implicit none(type, external)

    private
    public :: collect_Fourier_Transform_fft

contains

    subroutine collect_Fourier_Transform_fft(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("fft complex sp 1D", test_cmplx_fft_sp)&
                    ]

    end subroutine collect_Fourier_Transform_fft

    subroutine test_cmplx_fft_sp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 8
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method

        real(sp), dimension(Nx) :: xlist
        complex(sp), dimension(Nx) :: signal
        complex(sp), dimension(Nx) :: fs_exact, fs_fft
        integer :: i

        xlist = [(real((i - 1), kind=sp), i=1, Nx)]
        xlist = xlist * 2 * pi_sp / Nx

        signal = (0._sp, 0._sp)
        signal(1) = (1._sp, 0._sp)

        fs_exact = (1._sp, 0._sp)

        loop_method = init_loop_method(use_do_classic=.true.)
        call FT%init_fft_plan(Nx)
        fs_fft = FT%fft(signal, loop_method)
        call check(error, maxval(abs(fs_exact - fs_fft)) < TOL_TEST_sp)
        if (allocated(error)) return
        call FT%destroy_fft_plan()

        loop_method = init_loop_method(use_do_concurrent=.true.)
        call FT%init_fft_plan(Nx)
        fs_fft = FT%fft(signal, loop_method)
        call check(error, maxval(abs(fs_exact - fs_fft)) < TOL_TEST_sp)
        if (allocated(error)) return
        call FT%destroy_fft_plan()

        loop_method = init_loop_method(use_vectorized=.true.)
        call FT%init_fft_plan(Nx)
        fs_fft = FT%fft(signal, loop_method)
        call check(error, maxval(abs(fs_exact - fs_fft)) < TOL_TEST_sp)
        if (allocated(error)) return
        call FT%destroy_fft_plan()

        loop_method = init_loop_method(use_openmp=.true.)
        call FT%init_fft_plan(Nx)
        fs_fft = FT%fft(signal, loop_method)
        call check(error, maxval(abs(fs_exact - fs_fft)) < TOL_TEST_sp)
        if (allocated(error)) return
        call FT%destroy_fft_plan()

    end subroutine test_cmplx_fft_sp

end module test_Fourier_Transform_fft_cmplx

program test
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_Fourier_Transform_fft_cmplx

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("Fourier Transform fft 1D complex", collect_Fourier_Transform_fft) &
                 ]

    call init_color_output(.true.)

    do is = 1, size(testsuites)
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(signal) failed!"
        error stop
    end if

end program test
