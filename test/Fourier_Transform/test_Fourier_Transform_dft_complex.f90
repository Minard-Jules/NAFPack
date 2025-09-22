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

        integer, parameter :: Nx = 5, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(sp), dimension(Nx) :: xlist
        complex(sp), dimension(number_tests, Nx) :: signal
        complex(sp), dimension(Nx) :: signal_used
        complex(sp), dimension(number_tests, Nx) :: fs_exact
        complex(sp), dimension(Nx) :: fs_DFT
        integer :: i, freq

        xlist = [(real((i - 1), kind=sp), i=1, Nx)]
        xlist = xlist * 2 * pi_sp / Nx

        ! test 1: Dirac delta function
        signal(1, :) = (0._sp, 0._sp)
        signal(1, 1) = (1._sp, 0._sp)

        fs_exact(1, :) = (1._sp, 0._sp)

        !test 2: sinusoidal signal of frequency 2
        freq = 2
        signal(2, :) = exp(im_sp * freq * xlist)

        fs_exact(2, :) = (0._sp, 0._sp)
        fs_exact(2, freq + 1) = cmplx(real(Nx, kind=sp), 0._sp, kind=sp)

        ! test 3: constant signal (DC component)
        signal(3, :) = (1._sp, 0._sp)

        fs_exact(3, :) = (0._sp, 0._sp)
        fs_exact(3, 1) = cmplx(real(Nx, kind=sp), 0._sp, kind=sp)

        do i = 1, number_tests
            signal_used = signal(i, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_sp)
            if (allocated(error)) return
        end do

    end subroutine test_cmplx_dft_sp

    subroutine test_cmplx_idft_sp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(sp), dimension(Nx) :: xlist
        complex(sp), dimension(number_tests, Nx) :: fs
        complex(sp), dimension(Nx) :: fs_used
        complex(sp), dimension(number_tests, Nx) :: signal_exact
        complex(sp), dimension(Nx) :: signal_IDFT
        integer :: i, freq

        xlist = [(real((i - 1), kind=sp), i=1, Nx)]
        xlist = xlist * 2 * pi_sp / Nx

        ! test 1: Dirac delta function
        fs = (1._sp, 0._sp)

        signal_exact = (0._sp, 0._sp)
        signal_exact(1, 1) = (1._sp, 0._sp)

        !test 2: sinusoidal signal of frequency 2
        freq = 2
        fs(2, :) = (0._sp, 0._sp)
        fs(2, freq + 1) = cmplx(real(Nx, kind=sp), 0._sp, kind=sp)

        signal_exact(2, :) = exp(im_sp * freq * xlist)

        ! test 3: constant signal (DC component)
        fs(3, :) = (0._sp, 0._sp)
        fs(3, 1) = cmplx(real(Nx, kind=sp), 0._sp, kind=sp)

        signal_exact(3, :) = (1._sp, 0._sp)

        do i = 1, number_tests
            fs_used = fs(i, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_sp)
            if (allocated(error)) return
        end do

    end subroutine test_cmplx_idft_sp

    subroutine test_cmplx_dft_dp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(dp), dimension(Nx) :: xlist
        complex(dp), dimension(number_tests, Nx) :: signal
        complex(dp), dimension(Nx) :: signal_used
        complex(dp), dimension(number_tests, Nx) :: fs_exact
        complex(dp), dimension(Nx) :: fs_DFT
        integer :: i, freq

        xlist = [(real((i - 1), kind=dp), i=1, Nx)]
        xlist = xlist * 2 * pi_dp / Nx

        signal(1, :) = (0._dp, 0._dp)
        signal(1, 1) = (1._dp, 0._dp)
        fs_exact(1, :) = (1._dp, 0._dp)

        freq = 2
        signal(2, :) = exp(im_dp * freq * xlist)
        fs_exact(2, :) = (0._dp, 0._dp)
        fs_exact(2, freq + 1) = cmplx(real(Nx, kind=dp), 0._dp, kind=dp)

        signal(3, :) = (1._dp, 0._dp)
        fs_exact(3, :) = (0._dp, 0._dp)
        fs_exact(3, 1) = cmplx(real(Nx, kind=dp), 0._dp, kind=dp)

        do i = 1, number_tests
            signal_used = signal(i, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_dp)
            if (allocated(error)) return
        end do
    end subroutine test_cmplx_dft_dp

    subroutine test_cmplx_idft_dp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(dp), dimension(Nx) :: xlist
        complex(dp), dimension(number_tests, Nx) :: fs
        complex(dp), dimension(Nx) :: fs_used
        complex(dp), dimension(number_tests, Nx) :: signal_exact
        complex(dp), dimension(Nx) :: signal_IDFT
        integer :: i, freq

        xlist = [(real((i - 1), kind=dp), i=1, Nx)]
        xlist = xlist * 2 * pi_dp / Nx

        fs = (1._dp, 0._dp)
        signal_exact = (0._dp, 0._dp)
        signal_exact(1, 1) = (1._dp, 0._dp)

        freq = 2
        fs(2, :) = (0._dp, 0._dp)
        fs(2, freq + 1) = cmplx(real(Nx, kind=dp), 0._dp, kind=dp)
        signal_exact(2, :) = exp(im_dp * freq * xlist)

        fs(3, :) = (0._dp, 0._dp)
        fs(3, 1) = cmplx(real(Nx, kind=dp), 0._dp, kind=dp)
        signal_exact(3, :) = (1._dp, 0._dp)

        do i = 1, number_tests
            fs_used = fs(i, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_dp)
            if (allocated(error)) return
        end do
    end subroutine test_cmplx_idft_dp

    subroutine test_cmplx_dft_qp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(qp), dimension(Nx) :: xlist
        complex(qp), dimension(number_tests, Nx) :: signal
        complex(qp), dimension(Nx) :: signal_used
        complex(qp), dimension(number_tests, Nx) :: fs_exact
        complex(qp), dimension(Nx) :: fs_DFT
        integer :: i, freq

        xlist = [(real((i - 1), kind=qp), i=1, Nx)]
        xlist = xlist * 2 * pi_qp / Nx

        signal(1, :) = (0._qp, 0._qp)
        signal(1, 1) = (1._qp, 0._qp)
        fs_exact(1, :) = (1._qp, 0._qp)

        freq = 2
        signal(2, :) = exp(im_qp * freq * xlist)
        fs_exact(2, :) = (0._qp, 0._qp)
        fs_exact(2, freq + 1) = cmplx(real(Nx, kind=qp), 0._qp, kind=qp)

        signal(3, :) = (1._qp, 0._qp)
        fs_exact(3, :) = (0._qp, 0._qp)
        fs_exact(3, 1) = cmplx(real(Nx, kind=qp), 0._qp, kind=qp)

        do i = 1, number_tests
            signal_used = signal(i, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true.)
            fs_DFT = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :) - fs_DFT)) < TOL_TEST_qp)
            if (allocated(error)) return
        end do
    end subroutine test_cmplx_dft_qp

    subroutine test_cmplx_idft_qp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(qp), dimension(Nx) :: xlist
        complex(qp), dimension(number_tests, Nx) :: fs
        complex(qp), dimension(Nx) :: fs_used
        complex(qp), dimension(number_tests, Nx) :: signal_exact
        complex(qp), dimension(Nx) :: signal_IDFT
        integer :: i, freq

        xlist = [(real((i - 1), kind=qp), i=1, Nx)]
        xlist = xlist * 2 * pi_qp / Nx

        fs = (1._qp, 0._qp)
        signal_exact = (0._qp, 0._qp)
        signal_exact(1, 1) = (1._qp, 0._qp)

        freq = 2
        fs(2, :) = (0._qp, 0._qp)
        fs(2, freq + 1) = cmplx(real(Nx, kind=qp), 0._qp, kind=qp)
        signal_exact(2, :) = exp(im_qp * freq * xlist)

        fs(3, :) = (0._qp, 0._qp)
        fs(3, 1) = cmplx(real(Nx, kind=qp), 0._qp, kind=qp)
        signal_exact(3, :) = (1._qp, 0._qp)

        do i = 1, number_tests
            fs_used = fs(i, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true.)
            signal_IDFT = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_qp)
            if (allocated(error)) return
        end do
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
