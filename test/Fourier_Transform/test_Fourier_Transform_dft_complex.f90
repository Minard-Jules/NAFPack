module test_Fourier_Transform_dft_cmplx

    use NAFPack_kinds, only: dp, sp, qp
    use NAFPack_constant, only: TOL_TEST_sp, TOL_TEST_dp, TOL_TEST_qp
    use NAFPack_constant, only: pi_sp, im_sp, pi_dp, im_dp, pi_qp, im_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_loop_method, only: LoopMethod, init_loop_method
    use NAFPack_Fourier_Transform, only: Fourier_Transform

    implicit none(type, external)

    private
    public :: collect_Fourier_Transform_dft

contains

    subroutine collect_Fourier_Transform_dft(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("dft complex sp 1D", test_cmplx_dft_sp), &
                    new_unittest("idft complex sp 1D", test_cmplx_idft_sp), &
                    new_unittest("dft complex dp 1D", test_cmplx_dft_dp), &
                    new_unittest("idft complex dp 1D", test_cmplx_idft_dp), &
                    new_unittest("dft complex qp 1D", test_cmplx_dft_qp), &
                    new_unittest("idft complex qp 1D", test_cmplx_idft_qp) &
                    ]

    end subroutine collect_Fourier_Transform_dft

    subroutine test_cmplx_dft_sp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method

        real(sp), dimension(Nx) :: xlist
        complex(sp), dimension(Nx) :: signal
        complex(sp), dimension(Nx) :: fs_exact, fs_DFT
        integer :: i

        xlist = [(real((i - 1), kind=sp), i=1, Nx)]
        xlist = xlist * 2 * pi_sp / Nx

        signal = (0._sp, 0._sp)
        signal(1) = (1._sp, 0._sp)

        fs_exact = (1._sp, 0._sp)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_sp)
        if (allocated(error)) return

    end subroutine test_cmplx_dft_sp

    subroutine test_cmplx_idft_sp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method

        real(sp), dimension(Nx) :: xlist
        complex(sp), dimension(Nx) :: fs
        complex(sp), dimension(Nx) :: signal_exact, signal_IDFT
        integer :: i

        xlist = [(real((i - 1), kind=sp), i=1, Nx)]
        xlist = xlist * 2 * pi_sp / Nx

        fs = (1._sp, 0._sp)

        signal_exact = (0._sp, 0._sp)
        signal_exact(1) = (1._sp, 0._sp)

        loop_method = init_loop_method(use_do_classic=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_sp)
        if (allocated(error)) return

    end subroutine test_cmplx_idft_sp

    subroutine test_cmplx_dft_dp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method

        real(dp), dimension(Nx) :: xlist
        complex(dp), dimension(Nx) :: signal
        complex(dp), dimension(Nx) :: fs_exact, fs_DFT
        integer :: i

        xlist = [(real((i - 1), kind=sp), i=1, Nx)]
        xlist = xlist * 2 * pi_dp / Nx

        signal = (0._dp, 0._dp)
        signal(1) = (1._dp, 0._dp)

        fs_exact = (1._dp, 0._dp)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_dp)
        if (allocated(error)) return

    end subroutine test_cmplx_dft_dp

    subroutine test_cmplx_idft_dp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method

        real(dp), dimension(Nx) :: xlist
        complex(dp), dimension(Nx) :: fs
        complex(dp), dimension(Nx) :: signal_exact, signal_IDFT
        integer :: i

        xlist = [(real((i - 1), kind=dp), i=1, Nx)]
        xlist = xlist * 2 * pi_dp / Nx

        fs = (1._dp, 0._dp)

        signal_exact = (0._dp, 0._dp)
        signal_exact(1) = (1._dp, 0._dp)

        loop_method = init_loop_method(use_do_classic=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_dp)
        if (allocated(error)) return

    end subroutine test_cmplx_idft_dp

    subroutine test_cmplx_dft_qp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method

        real(qp), dimension(Nx) :: xlist
        complex(qp), dimension(Nx) :: signal
        complex(qp), dimension(Nx) :: fs_exact, fs_DFT
        integer :: i

        xlist = [(real((i - 1), kind=sp), i=1, Nx)]
        xlist = xlist * 2 * pi_qp / Nx

        signal = (0._qp, 0._qp)
        signal(1) = (1._qp, 0._qp)

        fs_exact = (1._qp, 0._qp)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true.)
        fs_DFT = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT)) < TOL_TEST_qp)
        if (allocated(error)) return

    end subroutine test_cmplx_dft_qp

    subroutine test_cmplx_idft_qp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method

        real(qp), dimension(Nx) :: xlist
        complex(qp), dimension(Nx) :: fs
        complex(qp), dimension(Nx) :: signal_exact, signal_IDFT
        integer :: i

        xlist = [(real((i - 1), kind=qp), i=1, Nx)]
        xlist = xlist * 2 * pi_qp / Nx

        fs = (1._qp, 0._qp)

        signal_exact = (0._qp, 0._qp)
        signal_exact(1) = (1._qp, 0._qp)

        loop_method = init_loop_method(use_do_classic=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true.)
        signal_IDFT = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT)) < TOL_TEST_qp)
        if (allocated(error)) return

    end subroutine test_cmplx_idft_qp

end module test_Fourier_Transform_dft_cmplx

program test
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_Fourier_Transform_dft_cmplx

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("Fourier Transform DFT 1D complex", collect_Fourier_Transform_dft) &
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
