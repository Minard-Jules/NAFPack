module test_Fourier_Transform_dft3_complex

    use NAFPack_kinds, only: dp, sp, qp
    use NAFPack_constant, only: TOL_TEST_sp, TOL_TEST_dp, TOL_TEST_qp
    use NAFPack_constant, only: pi_sp, im_sp, pi_dp, im_dp, pi_qp, im_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_loop_method, only: LoopMethod, init_loop_method
    use NAFPack_Fourier_Transform, only: Fourier_Transform

    implicit none(type, external)

    private
    public :: collect_Fourier_Transform_dft3

contains

    subroutine collect_Fourier_Transform_dft3(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("dft complex sp 3D", test_cmplx_dft3_sp), &
                    new_unittest("idft complex sp 3D", test_cmplx_idft3_sp), &
                    new_unittest("dft complex dp 3D", test_cmplx_dft3_dp), &
                    new_unittest("idft complex dp 3D", test_cmplx_idft3_dp), &
                    new_unittest("dft complex qp 3D", test_cmplx_dft3_qp), &
                    new_unittest("idft complex qp 3D", test_cmplx_idft3_qp) &
                    ]

    end subroutine collect_Fourier_Transform_dft3

    subroutine test_cmplx_dft3_sp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10, Nz = 15
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        complex(sp), dimension(Nx, Ny, Nz) :: signal
        complex(sp), dimension(Nx, Ny, Nz) :: fs_exact, fs_dft3

        signal = (0._sp, 0._sp)
        signal(1, 1, 1) = (1._sp, 0._sp)

        fs_exact = (1._sp, 0._sp)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_sp)
        if (allocated(error)) return

    end subroutine test_cmplx_dft3_sp

    subroutine test_cmplx_idft3_sp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10, Nz = 15
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        complex(sp), dimension(Nx, Ny, Nz) :: fs
        complex(sp), dimension(Nx, Ny, Nz) :: signal_exact, signal_Idft3

        fs = (1._sp, 0._sp)

        signal_exact = (0._sp, 0._sp)
        signal_exact(1, 1, 1) = (1._sp, 0._sp)

        loop_method = init_loop_method(use_do_classic=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_sp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_sp)
        if (allocated(error)) return

    end subroutine test_cmplx_idft3_sp

    subroutine test_cmplx_dft3_dp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10, Nz = 15
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        complex(dp), dimension(Nx, Ny, Nz) :: signal
        complex(dp), dimension(Nx, Ny, Nz) :: fs_exact, fs_dft3

        signal = (0._dp, 0._dp)
        signal(1, 1, 1) = (1._dp, 0._dp)

        fs_exact = (1._dp, 0._dp)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_dp)
        if (allocated(error)) return

    end subroutine test_cmplx_dft3_dp

    subroutine test_cmplx_idft3_dp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10, Nz = 15
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        complex(dp), dimension(Nx, Ny, Nz) :: fs
        complex(dp), dimension(Nx, Ny, Nz) :: signal_exact, signal_Idft3

        fs = (1._dp, 0._dp)

        signal_exact = (0._dp, 0._dp)
        signal_exact(1, 1, 1) = (1._dp, 0._dp)

        loop_method = init_loop_method(use_do_classic=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_dp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_dp)
        if (allocated(error)) return

    end subroutine test_cmplx_idft3_dp

    subroutine test_cmplx_dft3_qp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10, Nz = 15
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        complex(sp), dimension(Nx, Ny, Nz) :: signal
        complex(sp), dimension(Nx, Ny, Nz) :: fs_exact, fs_dft3

        signal = (0._qp, 0._qp)
        signal(1, 1, 1) = (1._qp, 0._qp)

        fs_exact = (1._qp, 0._qp)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true., num_threads=4)
        fs_dft3 = FT%dft(signal, loop_method=loop_method)
        call check(error, maxval(abs(fs_exact - fs_dft3)) < TOL_TEST_qp)
        if (allocated(error)) return

    end subroutine test_cmplx_dft3_qp

    subroutine test_cmplx_idft3_qp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Nx = 5, Ny = 10, Nz = 15
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        complex(qp), dimension(Nx, Ny, Nz) :: fs
        complex(qp), dimension(Nx, Ny, Nz) :: signal_exact, signal_Idft3

        fs = (1._qp, 0._qp)

        signal_exact = (0._qp, 0._qp)
        signal_exact(1, 1, 1) = (1._qp, 0._qp)

        loop_method = init_loop_method(use_do_classic=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_qp)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        signal_Idft3 = FT%idft(fs, loop_method=loop_method)
        call check(error, maxval(abs(signal_exact - signal_Idft3)) < TOL_TEST_qp)
        if (allocated(error)) return

    end subroutine test_cmplx_idft3_qp

end module test_Fourier_Transform_dft3_complex

program test
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_Fourier_Transform_dft3_complex

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("Fourier Transform DFT 3D complex", collect_Fourier_Transform_dft3) &
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
