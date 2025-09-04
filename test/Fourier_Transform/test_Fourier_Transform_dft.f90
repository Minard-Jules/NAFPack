module test_Fourier_Transform_dft

    use NAFPack_kinds, only: dp, sp, qp
    use NAFPack_constant, only: TOL_TEST, pi_sp, im_sp, pi_dp, im_dp, pi_qp, im_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_loop_method, only: LoopMethod, init_loop_method
    use NAFPack_Fourier_Transform, only: Fourier_Transform

    implicit none(type, external)

    private
    public :: collect_Fourier_Transform_dft

contains

    subroutine collect_Fourier_Transform_dft(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [&
            new_unittest("test_dft_dp", test_cmplx_dft_dp) &
        ]

    end subroutine collect_Fourier_Transform_dft

    subroutine test_cmplx_dft_dp(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: Mx = 10, My = 20
        integer, parameter :: Mx2 = 5, My2 = 6
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method

        real(dp), dimension(Mx) :: xlist
        real(dp), dimension(My) :: ylist
        complex(dp), dimension(Mx) :: S
        complex(dp), dimension(Mx) :: fs_exact_DFT, fs_DFT, diff_DFT
        integer :: i

        do i = 1, Mx
            xlist(i) = i - 1
        end do
        xlist = xlist * 2 * pi_dp / Mx

        do i = 1, My
            ylist(i) = i - 1
        end do
        ylist = ylist * 2 * pi_dp / My

        !==================================================================
        !test DFT
        S = exp(-im_dp * xlist)

        fs_exact_DFT = (0.d0, 0.d0)
        fs_exact_DFT(Mx) = (10.d0, 0.d0)

        loop_method = init_loop_method(use_do_classic=.true.)
        fs_DFT = FT%dft(S, loop_method=loop_method)

        diff_DFT = fs_exact_DFT - fs_DFT

        call check(error, maxval(abs(diff_DFT)) < TOL_TEST)
        if (allocated(error)) return

        loop_method = init_loop_method(use_vectorized=.true.)
        fs_DFT = FT%dft(S, loop_method=loop_method)

        diff_DFT = fs_exact_DFT - fs_DFT

        call check(error, maxval(abs(diff_DFT)) < TOL_TEST)
        if (allocated(error)) return

        loop_method = init_loop_method(use_do_concurrent=.true.)
        fs_DFT = FT%dft(S, loop_method=loop_method)

        diff_DFT = fs_exact_DFT - fs_DFT

        call check(error, maxval(abs(diff_DFT)) < TOL_TEST)
        if (allocated(error)) return

        loop_method = init_loop_method(use_openmp=.true.)
        fs_DFT = FT%dft(S, loop_method=loop_method)

        diff_DFT = fs_exact_DFT - fs_DFT

        call check(error, maxval(abs(diff_DFT)) < TOL_TEST)
        if (allocated(error)) return
    end subroutine test_cmplx_dft_dp
    
end module test_Fourier_Transform_dft

program test
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_Fourier_Transform_dft

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("Fourier Transform DFT", collect_Fourier_Transform_dft) &
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
