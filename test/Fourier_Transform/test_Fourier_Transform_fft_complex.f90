module test_Fourier_Transform_fft_cmplx

    use NAFPack_kinds, only: dp, sp, qp
    use NAFPack_constant, only: TOL_TEST_sp, TOL_TEST_dp, TOL_TEST_qp
    use NAFPack_constant, only: pi_sp, im_sp, pi_dp, im_dp, pi_qp, im_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_loop_method_type, only: LoopMethod
    use NAFPack_loop_method, only: init_loop_method
    use NAFPack_Fourier_Transform, only: Fourier_Transform, FFTAlgorithm, &
                                         ALG_RADIX2_DIT, ALG_RADIX2_DIF, &
                                         ALG_RADIX3_DIT, ALG_RADIX3_DIF, &
                                         ALG_RADIX4_DIT, ALG_RADIX4_DIF, &
                                         ALG_RADIX5_DIT, ALG_RADIX5_DIF, &
                                         ALG_MIXED_DIT, ALG_MIXED_DIF, &
                                         ALG_SPLIT_DIT, ALG_SPLIT_DIF
    use NAFPack_implementation_type, only: ImplementationType, ITERATIVE, recursive

    implicit none(type, external)

    private
    public :: collect_Fourier_Transform_fft

contains

    subroutine collect_Fourier_Transform_fft(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("fft iterative radix-2     DIT complex(sp) 1D", &
                                 test_fft_radix2_iterative_DIT_cmplx_sp), &
                    new_unittest("fft iterative radix-2     DIF complex(sp) 1D", &
                                 test_fft_radix2_iterative_DIF_cmplx_sp), &
                    new_unittest("fft recursive radix-2     DIT complex(sp) 1D", &
                                 test_fft_radix2_recursive_DIT_cmplx_sp), &
                    new_unittest("fft recursive radix-2     DIF complex(sp) 1D", &
                                 test_fft_radix2_recursive_DIF_cmplx_sp), &
                    new_unittest("fft iterative radix-3     DIT complex(sp) 1D", &
                                 test_fft_radix3_iterative_DIT_cmplx_sp), &
                    new_unittest("fft iterative radix-3     DIF complex(sp) 1D", &
                                 test_fft_radix3_iterative_DIF_cmplx_sp), &
                    new_unittest("fft recursive radix-3     DIT complex(sp) 1D", &
                                 test_fft_radix3_recursive_DIT_cmplx_sp), &
                    new_unittest("fft recursive radix-3     DIF complex(sp) 1D", &
                                 test_fft_radix3_recursive_DIF_cmplx_sp), &
                    new_unittest("fft iterative radix-4     DIT complex(sp) 1D", &
                                 test_fft_radix4_iterative_DIT_cmplx_sp), &
                    new_unittest("fft iterative radix-4     DIF complex(sp) 1D", &
                                 test_fft_radix4_iterative_DIF_cmplx_sp), &
                    new_unittest("fft recursive radix-4     DIT complex(sp) 1D", &
                                 test_fft_radix4_recursive_DIT_cmplx_sp), &
                    new_unittest("fft recursive radix-4     DIF complex(sp) 1D", &
                                 test_fft_radix4_recursive_DIF_cmplx_sp), &
                    new_unittest("fft iterative radix-5     DIT complex(sp) 1D", &
                                 test_fft_radix5_iterative_DIT_cmplx_sp), &
                    new_unittest("fft iterative radix-5     DIF complex(sp) 1D", &
                                 test_fft_radix5_iterative_DIF_cmplx_sp), &
                    new_unittest("fft recursive radix-5     DIT complex(sp) 1D", &
                                 test_fft_radix5_recursive_DIT_cmplx_sp), &
                    new_unittest("fft recursive radix-5     DIF complex(sp) 1D", &
                                 test_fft_radix5_recursive_DIF_cmplx_sp), &
                    new_unittest("fft iterative mixed radix DIT complex(sp) 1D", &
                                 test_fft_mixed_radix_iterative_DIT_cmplx_sp), &
                    new_unittest("fft iterative mixed radix DIF complex(sp) 1D", &
                                 test_fft_mixed_radix_iterative_DIF_cmplx_sp), &
                    new_unittest("fft recursive mixed radix DIT complex(sp) 1D", &
                                 test_fft_mixed_radix_recursive_DIT_cmplx_sp), &
                    new_unittest("fft recursive mixed radix DIF complex(sp) 1D", &
                                 test_fft_mixed_radix_recursive_DIF_cmplx_sp), &
                    new_unittest("fft iterative split radix DIT complex(sp) 1D", &
                                 test_fft_split_radix_iterative_DIT_cmplx_sp), &
                    new_unittest("fft iterative split radix DIF complex(sp) 1D", &
                                 test_fft_split_radix_iterative_DIF_cmplx_sp), &
                    new_unittest("fft recursive split radix DIT complex(sp) 1D", &
                                 test_fft_split_radix_recursive_DIT_cmplx_sp), &
                    new_unittest("fft recursive split radix DIF complex(sp) 1D", &
                                 test_fft_split_radix_recursive_DIF_cmplx_sp) &
                    ]

    end subroutine collect_Fourier_Transform_fft

    subroutine test_fft_radix2_iterative_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(5) :: loop_method
        integer, parameter :: Nx = 16
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)
        loop_method(2) = init_loop_method(use_openmp_simd=.true.)
        loop_method(3) = init_loop_method(use_array_syntax=.true.)
        loop_method(4) = init_loop_method(use_do_concurrent=.true.)
        loop_method(5) = init_loop_method(use_openmp=.true.)

        do i = 1, 5
            call test_fft_cmplx_sp(error, ALG_RADIX2_DIT, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_radix2_iterative_DIT_cmplx_sp

    subroutine test_fft_radix2_iterative_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(5) :: loop_method
        integer, parameter :: Nx = 16
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)
        loop_method(2) = init_loop_method(use_openmp_simd=.true.)
        loop_method(3) = init_loop_method(use_array_syntax=.true.)
        loop_method(4) = init_loop_method(use_do_concurrent=.true.)
        loop_method(5) = init_loop_method(use_openmp=.true.)

        do i = 1, 5
            call test_fft_cmplx_sp(error, ALG_RADIX2_DIF, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_radix2_iterative_DIF_cmplx_sp

    subroutine test_fft_radix2_recursive_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 16

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_RADIX2_DIT, Nx, recursive, loop_method)
    end subroutine test_fft_radix2_recursive_DIT_cmplx_sp

    subroutine test_fft_radix2_recursive_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 16

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_RADIX2_DIF, Nx, recursive, loop_method)
    end subroutine test_fft_radix2_recursive_DIF_cmplx_sp

    subroutine test_fft_radix3_iterative_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 27
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_RADIX3_DIT, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_radix3_iterative_DIT_cmplx_sp

    subroutine test_fft_radix3_iterative_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 27
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_RADIX3_DIF, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_radix3_iterative_DIF_cmplx_sp

    subroutine test_fft_radix3_recursive_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 27

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_RADIX3_DIT, Nx, recursive, loop_method)
    end subroutine test_fft_radix3_recursive_DIT_cmplx_sp

    subroutine test_fft_radix3_recursive_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 27

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_RADIX3_DIF, Nx, recursive, loop_method)
    end subroutine test_fft_radix3_recursive_DIF_cmplx_sp

    subroutine test_fft_radix4_iterative_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 16
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_RADIX4_DIT, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_radix4_iterative_DIT_cmplx_sp

    subroutine test_fft_radix4_iterative_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 16
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_RADIX4_DIF, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_radix4_iterative_DIF_cmplx_sp

    subroutine test_fft_radix4_recursive_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 16

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_RADIX4_DIT, Nx, recursive, loop_method)
    end subroutine test_fft_radix4_recursive_DIT_cmplx_sp

    subroutine test_fft_radix4_recursive_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 16

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_RADIX4_DIF, Nx, recursive, loop_method)
    end subroutine test_fft_radix4_recursive_DIF_cmplx_sp

    subroutine test_fft_radix5_iterative_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 25
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_RADIX5_DIT, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_radix5_iterative_DIT_cmplx_sp

    subroutine test_fft_radix5_iterative_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 25
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_RADIX5_DIF, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_radix5_iterative_DIF_cmplx_sp

    subroutine test_fft_radix5_recursive_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 25

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_RADIX5_DIT, Nx, recursive, loop_method)
    end subroutine test_fft_radix5_recursive_DIT_cmplx_sp

    subroutine test_fft_radix5_recursive_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 25

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_RADIX5_DIF, Nx, recursive, loop_method)
    end subroutine test_fft_radix5_recursive_DIF_cmplx_sp

    subroutine test_fft_mixed_radix_iterative_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 105
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_MIXED_DIT, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_mixed_radix_iterative_DIT_cmplx_sp

    subroutine test_fft_mixed_radix_iterative_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 105
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_MIXED_DIF, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_mixed_radix_iterative_DIF_cmplx_sp

    subroutine test_fft_mixed_radix_recursive_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 105

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_MIXED_DIT, Nx, recursive, loop_method)
    end subroutine test_fft_mixed_radix_recursive_DIT_cmplx_sp

    subroutine test_fft_mixed_radix_recursive_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 105

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_MIXED_DIF, Nx, recursive, loop_method)
    end subroutine test_fft_mixed_radix_recursive_DIF_cmplx_sp

    subroutine test_fft_split_radix_iterative_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 16
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_SPLIT_DIT, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_split_radix_iterative_DIT_cmplx_sp

    subroutine test_fft_split_radix_iterative_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod), dimension(1) :: loop_method
        integer, parameter :: Nx = 16
        integer :: i

        loop_method(1) = init_loop_method(use_do_classic=.true.)

        do i = 1, 1
            call test_fft_cmplx_sp(error, ALG_SPLIT_DIF, Nx, ITERATIVE, loop_method(i))
        end do
    end subroutine test_fft_split_radix_iterative_DIF_cmplx_sp

    subroutine test_fft_split_radix_recursive_DIT_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 16

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_SPLIT_DIT, Nx, recursive, loop_method)
    end subroutine test_fft_split_radix_recursive_DIT_cmplx_sp

    subroutine test_fft_split_radix_recursive_DIF_cmplx_sp(error)
        type(error_type), allocatable, intent(out) :: error
        type(LoopMethod) :: loop_method
        integer, parameter :: Nx = 16

        loop_method = init_loop_method()
        call test_fft_cmplx_sp(error, ALG_SPLIT_DIF, Nx, recursive, loop_method)
    end subroutine test_fft_split_radix_recursive_DIF_cmplx_sp

    subroutine test_fft_cmplx_sp(error, algorithm, Nx, implementation_type, loop_method)
        type(error_type), allocatable, intent(out) :: error
        type(FFTAlgorithm), intent(in) :: algorithm
        integer, intent(in) :: Nx
        type(ImplementationType), intent(in) :: implementation_type
        type(LoopMethod), intent(in) :: loop_method

        integer, parameter :: number_tests = 3
        type(Fourier_Transform) :: FT

        real(sp), dimension(Nx) :: xlist
        real(sp), dimension(number_tests, Nx) :: signal_re, signal_im
        real(sp), dimension(Nx) :: signal_used_re, signal_used_im
        real(sp), dimension(number_tests, Nx) :: fs_exact_re, fs_exact_im
        real(sp), dimension(Nx) :: fs_fft_re, fs_fft_im
        integer :: i, freq, j

        xlist = [(real((i - 1), kind=sp), i=1, Nx)]
        xlist = xlist * 2 * pi_sp / Nx

        ! test 1: Dirac delta function
        signal_re(1, :) = 0._sp
        signal_im(1, :) = 0._sp
        signal_re(1, 1) = 1._sp

        fs_exact_re(1, :) = 1._sp
        fs_exact_im(1, :) = 0._sp

        !test 2: sinusoidal signal of frequency 2
        freq = 2
        signal_re(2, :) = cos(freq * xlist)
        signal_im(2, :) = sin(freq * xlist)

        fs_exact_re(2, :) = 0._sp
        fs_exact_im(2, :) = 0._sp
        fs_exact_re(2, freq + 1) = real(Nx, kind=sp)
        fs_exact_im(2, freq + 1) = 0._sp

        ! test 3: constant signal (DC component)
        signal_re(3, :) = 1._sp
        signal_im(3, :) = 0._sp

        fs_exact_re(3, :) = 0._sp
        fs_exact_im(3, :) = 0._sp
        fs_exact_re(3, 1) = real(Nx, kind=sp)
        fs_exact_im(3, 1) = 0._sp

        ! call FT%init_fft_plan(Nx, algorithm=algorithm)
        call FT%create_fft_plan_sp(Nx, algorithm=algorithm)
        do i = 1, number_tests
            signal_used_re = signal_re(i, :)
            signal_used_im = signal_im(i, :)
            ! fs_fft = FT%fft(signal_used, loop_method, implementation_type)
            call FT%execute_fft_plan_sp(signal_used_re, signal_used_im, &
                                        fs_fft_re, fs_fft_im, &
                                        loop_method=loop_method, implementation_type=implementation_type)
            ! print*, "  compute re     ", "  exact re       ", '|', "   compute im       ", "exact im"
            ! do j = 1, Nx
            !     print*, fs_fft_re(j), fs_exact_re(i, j), '|', fs_fft_im(j), fs_exact_im(i, j)
            ! end do
            ! print*, ""
            call check(error, maxval(abs(fs_exact_re(i, :) - fs_fft_re)) < TOL_TEST_sp)
            call check(error, maxval(abs(fs_exact_im(i, :) - fs_fft_im)) < TOL_TEST_sp)
            if (allocated(error)) return
        end do
        call FT%destroy_fft_plan_sp()

    end subroutine test_fft_cmplx_sp

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

    deallocate (testsuites)

end program test
