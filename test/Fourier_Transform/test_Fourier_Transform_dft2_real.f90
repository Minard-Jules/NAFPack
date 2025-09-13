module test_Fourier_Transform_dft2_real

    use NAFPack_kinds, only: dp, sp, qp
    use NAFPack_constant, only: TOL_TEST_sp, TOL_TEST_dp, TOL_TEST_qp
    use NAFPack_constant, only: pi_sp, im_sp, pi_dp, im_dp, pi_qp, im_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_loop_method, only: LoopMethod, init_loop_method
    use NAFPack_Fourier_Transform, only: Fourier_Transform

    implicit none(type, external)

    private
    public :: collect_Fourier_Transform_dft2

contains

    subroutine collect_Fourier_Transform_dft2(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("dft real sp 2D", test_real_dft2_sp), &
                    new_unittest("idft real sp 2D", test_real_idft2_sp), &
                    new_unittest("dft real dp 2D", test_real_dft2_dp), &
                    new_unittest("idft real dp 2D", test_real_idft2_dp), &
                    new_unittest("dft real qp 2D", test_real_dft2_qp), &
                    new_unittest("idft real qp 2D", test_real_idft2_qp) &
                    ]

    end subroutine collect_Fourier_Transform_dft2

    subroutine test_real_dft2_sp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(sp), dimension(Nx, Ny) :: signal
        complex(sp), dimension(Nx, Ny) :: fs_exact, fs_DFT2

        signal = 0._sp
        signal(1, 1) = 1._sp

        fs_exact = (1._sp, 0._sp)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_sp)
        if (allocated(error)) return

    end subroutine test_real_dft2_sp

    subroutine test_real_idft2_sp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(sp), dimension(Nx, Ny) :: fs
        complex(sp), dimension(Nx, Ny) :: signal_exact, signal_IDFT2

        fs = 1._sp

        signal_exact = (0._sp, 0._sp)
        signal_exact(1, 1) = (1._sp, 0._sp)

        loop_method = init_loop_method(use_do_classic=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_sp)
        if (allocated(error)) return

    end subroutine test_real_idft2_sp

    subroutine test_real_dft2_dp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(dp), dimension(Nx, Ny) :: signal
        complex(dp), dimension(Nx, Ny) :: fs_exact, fs_DFT2

        signal = 0._dp
        signal(1, 1) = 1._dp

        fs_exact = (1._dp, 0._dp)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_dp)
        if (allocated(error)) return

    end subroutine test_real_dft2_dp

    subroutine test_real_idft2_dp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(dp), dimension(Nx, Ny) :: fs
        complex(dp), dimension(Nx, Ny) :: signal_exact, signal_IDFT2

        fs = 1._dp

        signal_exact = (0._dp, 0._dp)
        signal_exact(1, 1) = (1._dp, 0._dp)

        loop_method = init_loop_method(use_do_classic=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_dp)
        if (allocated(error)) return

    end subroutine test_real_idft2_dp

    subroutine test_real_dft2_qp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(sp), dimension(Nx, Ny) :: signal
        complex(sp), dimension(Nx, Ny) :: fs_exact, fs_DFT2

        signal = 0._qp
        signal(1, 1) = 1._qp

        fs_exact = (1._qp, 0._qp)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        fs_DFT2 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_DFT2)) < TOL_TEST_qp)
        if (allocated(error)) return

    end subroutine test_real_dft2_qp

    subroutine test_real_idft2_qp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(qp), dimension(Nx, Ny) :: fs
        complex(qp), dimension(Nx, Ny) :: signal_exact, signal_IDFT2

        fs = 1._qp

        signal_exact = (0._qp, 0._qp)
        signal_exact(1, 1) = (1._qp, 0._qp)

        loop_method = init_loop_method(use_do_classic=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_IDFT2 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_IDFT2)) < TOL_TEST_qp)
        if (allocated(error)) return

    end subroutine test_real_idft2_qp

end module test_Fourier_Transform_dft2_real

program test
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_Fourier_Transform_dft2_real

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("Fourier Transform DFT 2D real", collect_Fourier_Transform_dft2) &
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