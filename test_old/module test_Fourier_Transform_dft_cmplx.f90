module test_dft_cmplx

    use testdrive, only: new_unittest, unittest_type, error_type, check
    use DFT_naive, only: dft
    use kinds, only: dp

    implicit none(type, external)

    private
    public :: collect_dft

contains

    subroutine collect_dft(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("dft complex(dp)  1D", test_cmplx_dft_dp), &
                    new_unittest("idft complex(dp) 1D", test_cmplx_idft_dp)
                    ]

    end subroutine collect_dft

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
            fs_DFT = dft(signal_used)
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
            signal_IDFT = FT%idft(fs_used)
            call check(error, maxval(abs(signal_exact(i, :) - signal_IDFT)) < TOL_TEST_dp)
            if (allocated(error)) return
        end do
    end subroutine test_cmplx_idft_dp

end module test_dft_cmplx

program test
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_dft_cmplx

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

    deallocate(testsuites)

end program test