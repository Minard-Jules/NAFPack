module NAFPack_Fourier_Transform

    use NAFPack_kinds, only: dp, sp, qp, i8, i16, isp, idp
    use NAFPack_constant, only: pi_sp, pi_dp, pi_qp, im_sp, im_dp, im_qp
    use NAFPack_loop_method_type, only: LoopMethod, default_loop_method
    use NAFPack_loop_method, only: check_loop_method
    use NAFPack_math_utils, only: sieve_of_eratosthenes, is_power_of_two, power_of_p_exponent
    use NAFPack_implementation_type, only: ImplementationType, recursive, ITERATIVE, DEFAULT_IMPLEMENTATION_TYPE
    use NAFPack_plan_type, only: DecimationMethod, DIT, DIF
    use NAFPack_plan_type, only: FFTAlgorithm, ALG_NONE, ALG_AUTO, ALG_RADIX2_DIT, ALG_RADIX2_DIF, &
                                 ALG_RADIX3_DIT, ALG_RADIX3_DIF, &
                                 ALG_RADIX4_DIT, ALG_RADIX4_DIF, &
                                 ALG_RADIX5_DIT, ALG_RADIX5_DIF, &
                                 ALG_MIXED_DIT, ALG_MIXED_DIF, ALG_SPLIT_DIT, ALG_SPLIT_DIF
    use NAFPack_plan_type, only: FFTPlan, StageView

    implicit none(type, external)

    private
    public :: Fourier_Transform
    public :: dft, idft
    public :: dft2, idft2
    public :: dft3, idft3

    public :: DecimationMethod, DIT, DIF
    public :: FFTAlgorithm, ALG_NONE, ALG_AUTO, &
              ALG_RADIX2_DIT, ALG_RADIX2_DIF, &
              ALG_RADIX3_DIT, ALG_RADIX3_DIF, &
              ALG_RADIX4_DIT, ALG_RADIX4_DIF, &
              ALG_RADIX5_DIT, ALG_RADIX5_DIF, &
              ALG_MIXED_DIT, ALG_MIXED_DIF, ALG_SPLIT_DIT, ALG_SPLIT_DIF

    type :: Fourier_Transform
        type(FFTPlan) :: plan
    contains
        procedure :: create_fft_plan_sp
        procedure :: execute_fft_plan_sp
        procedure :: destroy_fft_plan_sp

        procedure, nopass, private :: dft_cmplx_sp, dft_cmplx_dp, dft_cmplx_qp
        generic :: dft => dft_cmplx_sp, dft_cmplx_dp, dft_cmplx_qp
        procedure, nopass, private :: dft_real_sp, dft_real_dp, dft_real_qp
        generic :: dft => dft_real_sp, dft_real_dp, dft_real_qp
        procedure, nopass, private :: dft2_cmplx_sp, dft2_cmplx_dp, dft2_cmplx_qp
        generic :: dft => dft2_cmplx_sp, dft2_cmplx_dp, dft2_cmplx_qp
        procedure, nopass, private :: dft2_real_sp, dft2_real_dp, dft2_real_qp
        generic :: dft => dft2_real_sp, dft2_real_dp, dft2_real_qp
        procedure, nopass, private :: dft3_cmplx_sp, dft3_cmplx_dp, dft3_cmplx_qp
        generic :: dft => dft3_cmplx_sp, dft3_cmplx_dp, dft3_cmplx_qp
        procedure, nopass, private :: dft3_real_sp, dft3_real_dp, dft3_real_qp
        generic :: dft => dft3_real_sp, dft3_real_dp, dft3_real_qp

        procedure, nopass, private :: idft_cmplx_sp, idft_cmplx_dp, idft_cmplx_qp
        generic :: idft => idft_cmplx_sp, idft_cmplx_dp, idft_cmplx_qp
        procedure, nopass, private :: idft_real_sp, idft_real_dp, idft_real_qp
        generic :: idft => idft_real_sp, idft_real_dp, idft_real_qp
        procedure, nopass, private :: idft2_cmplx_sp, idft2_cmplx_dp, idft2_cmplx_qp
        generic :: idft => idft2_cmplx_sp, idft2_cmplx_dp, idft2_cmplx_qp
        procedure, nopass, private :: idft2_real_sp, idft2_real_dp, idft2_real_qp
        generic :: idft => idft2_real_sp, idft2_real_dp, idft2_real_qp
        procedure, nopass, private :: idft3_cmplx_sp, idft3_cmplx_dp, idft3_cmplx_qp
        generic :: idft => idft3_cmplx_sp, idft3_cmplx_dp, idft3_cmplx_qp
        procedure, nopass, private :: idft3_real_sp, idft3_real_dp, idft3_real_qp
        generic :: idft => idft3_real_sp, idft3_real_dp, idft3_real_qp
    end type Fourier_Transform

    interface dft
        module procedure dft_cmplx_sp, dft_cmplx_dp, dft_cmplx_qp
        module procedure dft_real_sp, dft_real_dp, dft_real_qp
    end interface dft

    interface idft
        module procedure idft_cmplx_sp, idft_cmplx_dp, idft_cmplx_qp
        module procedure idft_real_sp, idft_real_dp, idft_real_qp
    end interface idft

    interface dft2
        module procedure dft2_cmplx_sp, dft2_cmplx_dp, dft2_cmplx_qp
        module procedure dft2_real_sp, dft2_real_dp, dft2_real_qp
    end interface dft2

    interface idft2
        module procedure idft2_cmplx_sp, idft2_cmplx_dp, idft2_cmplx_qp
        module procedure idft2_real_sp, idft2_real_dp, idft2_real_qp
    end interface idft2

    interface dft3
        module procedure dft3_cmplx_sp, dft3_cmplx_dp, dft3_cmplx_qp
        module procedure dft3_real_sp, dft3_real_dp, dft3_real_qp
    end interface dft3

    interface idft3
        module procedure idft3_cmplx_sp, idft3_cmplx_dp, idft3_cmplx_qp
        module procedure idft3_real_sp, idft3_real_dp, idft3_real_qp
    end interface idft3

    !====================================================================================
    ! FFT and IFFT for 1D, 2D, and 3D signals
    !====================================================================================

    interface
        pure module subroutine create_fft_plan_sp(this, N, algorithm, decimation_method)
            class(Fourier_Transform), intent(inout) :: this
            integer(isp), intent(in) :: N
            type(FFTAlgorithm), optional, intent(in) :: algorithm
            type(DecimationMethod), optional, intent(in) :: decimation_method
        end subroutine create_fft_plan_sp
    end interface

    interface
        module subroutine execute_fft_plan_sp(this, signal_re, signal_im, result_re, result_im, loop_method, implementation_type)
            class(Fourier_Transform), intent(in) :: this
            real(sp), dimension(:), intent(in) :: signal_re, signal_im
            real(sp), dimension(:), intent(out) :: result_re, result_im
            type(LoopMethod), optional, intent(in) :: loop_method
            type(ImplementationType), optional, intent(in) :: implementation_type
        end subroutine execute_fft_plan_sp
    end interface

    interface
        pure module subroutine destroy_fft_plan_sp(this)
            class(Fourier_Transform), intent(inout) :: this
        end subroutine destroy_fft_plan_sp
    end interface

    !====================================================================================
    ! DFT and IDFT for 1D, 2D, and 3D signals
    !====================================================================================

    interface
        module function dft_cmplx_sp(signal, loop_method) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:), allocatable :: result
        end function dft_cmplx_sp

        module function dft_cmplx_dp(signal, loop_method) result(result)
            complex(dp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:), allocatable :: result
        end function dft_cmplx_dp

        module function dft_cmplx_qp(signal, loop_method) result(result)
            complex(qp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:), allocatable :: result
        end function dft_cmplx_qp
    end interface
    interface
        module function idft_cmplx_sp(f_signal, loop_method) result(result)
            complex(sp), dimension(:), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:), allocatable :: result
        end function idft_cmplx_sp

        module function idft_cmplx_dp(f_signal, loop_method) result(result)
            complex(dp), dimension(:), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:), allocatable :: result
        end function idft_cmplx_dp

        module function idft_cmplx_qp(f_signal, loop_method) result(result)
            complex(qp), dimension(:), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:), allocatable :: result
        end function idft_cmplx_qp
    end interface

    interface
        module function dft_real_sp(signal, loop_method) result(result)
            real(sp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:), allocatable :: result
        end function dft_real_sp

        module function dft_real_dp(signal, loop_method) result(result)
            real(dp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:), allocatable :: result
        end function dft_real_dp

        module function dft_real_qp(signal, loop_method) result(result)
            real(qp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:), allocatable :: result
        end function dft_real_qp
    end interface
    interface
        module function idft_real_sp(f_signal, loop_method) result(result)
            real(sp), dimension(:), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:), allocatable :: result
        end function idft_real_sp

        module function idft_real_dp(f_signal, loop_method) result(result)
            real(dp), dimension(:), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:), allocatable :: result
        end function idft_real_dp

        module function idft_real_qp(f_signal, loop_method) result(result)
            real(qp), dimension(:), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:), allocatable :: result
        end function idft_real_qp
    end interface

    interface
        module function dft2_cmplx_sp(signal, loop_method) result(result)
            complex(sp), dimension(:, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:, :), allocatable :: result
        end function dft2_cmplx_sp

        module function dft2_cmplx_dp(signal, loop_method) result(result)
            complex(dp), dimension(:, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:, :), allocatable :: result
        end function dft2_cmplx_dp

        module function dft2_cmplx_qp(signal, loop_method) result(result)
            complex(qp), dimension(:, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:, :), allocatable :: result
        end function dft2_cmplx_qp
    end interface
    interface
        module function idft2_cmplx_sp(f_signal, loop_method) result(result)
            complex(sp), dimension(:, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:, :), allocatable :: result
        end function idft2_cmplx_sp

        module function idft2_cmplx_dp(f_signal, loop_method) result(result)
            complex(dp), dimension(:, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:, :), allocatable :: result
        end function idft2_cmplx_dp

        module function idft2_cmplx_qp(f_signal, loop_method) result(result)
            complex(qp), dimension(:, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:, :), allocatable :: result
        end function idft2_cmplx_qp
    end interface

    interface
        module function dft2_real_sp(signal, loop_method) result(result)
            real(sp), dimension(:, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:, :), allocatable :: result
        end function dft2_real_sp

        module function dft2_real_dp(signal, loop_method) result(result)
            real(dp), dimension(:, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:, :), allocatable :: result
        end function dft2_real_dp

        module function dft2_real_qp(signal, loop_method) result(result)
            real(qp), dimension(:, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:, :), allocatable :: result
        end function dft2_real_qp
    end interface
    interface
        module function idft2_real_sp(f_signal, loop_method) result(result)
            real(sp), dimension(:, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:, :), allocatable :: result
        end function idft2_real_sp

        module function idft2_real_dp(f_signal, loop_method) result(result)
            real(dp), dimension(:, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:, :), allocatable :: result
        end function idft2_real_dp

        module function idft2_real_qp(f_signal, loop_method) result(result)
            real(qp), dimension(:, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:, :), allocatable :: result
        end function idft2_real_qp
    end interface

    interface
        module function dft3_cmplx_sp(signal, loop_method) result(result)
            complex(sp), dimension(:, :, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:, :, :), allocatable :: result
        end function dft3_cmplx_sp

        module function dft3_cmplx_dp(signal, loop_method) result(result)
            complex(dp), dimension(:, :, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:, :, :), allocatable :: result
        end function dft3_cmplx_dp

        module function dft3_cmplx_qp(signal, loop_method) result(result)
            complex(qp), dimension(:, :, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:, :, :), allocatable :: result
        end function dft3_cmplx_qp
    end interface
    interface
        module function idft3_cmplx_sp(f_signal, loop_method) result(result)
            complex(sp), dimension(:, :, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:, :, :), allocatable :: result
        end function idft3_cmplx_sp

        module function idft3_cmplx_dp(f_signal, loop_method) result(result)
            complex(dp), dimension(:, :, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:, :, :), allocatable :: result
        end function idft3_cmplx_dp

        module function idft3_cmplx_qp(f_signal, loop_method) result(result)
            complex(qp), dimension(:, :, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:, :, :), allocatable :: result
        end function idft3_cmplx_qp
    end interface

    interface
        module function dft3_real_sp(signal, loop_method) result(result)
            real(sp), dimension(:, :, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:, :, :), allocatable :: result
        end function dft3_real_sp

        module function dft3_real_dp(signal, loop_method) result(result)
            real(dp), dimension(:, :, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:, :, :), allocatable :: result
        end function dft3_real_dp

        module function dft3_real_qp(signal, loop_method) result(result)
            real(qp), dimension(:, :, :), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:, :, :), allocatable :: result
        end function dft3_real_qp
    end interface
    interface
        module function idft3_real_sp(f_signal, loop_method) result(result)
            real(sp), dimension(:, :, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:, :, :), allocatable :: result
        end function idft3_real_sp

        module function idft3_real_dp(f_signal, loop_method) result(result)
            real(dp), dimension(:, :, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(:, :, :), allocatable :: result
        end function idft3_real_dp

        module function idft3_real_qp(f_signal, loop_method) result(result)
            real(qp), dimension(:, :, :), intent(in) :: f_signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(:, :, :), allocatable :: result
        end function idft3_real_qp
    end interface

end module NAFPack_Fourier_Transform
