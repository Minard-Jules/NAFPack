module NAFPack_Fourier_Transform

    use NAFPack_kinds, only: dp, sp, qp, i8, i16, isp, idp
    use NAFPack_constant, only: pi_sp, pi_dp, pi_qp, im_sp, im_dp, im_qp
    use NAFPack_loop_method, only: LoopMethod, default_loop_method, check_loop_method
    use NAFPack_math_utils, only: sieve_of_eratosthenes, is_power_of_two, power_of_p_exponent

    implicit none(type, external)

    private
    public :: Fourier_Transform
    public :: dft, idft
    public :: dft2, idft2

    type :: Twiddles_sp
        integer(isp) :: block_size = 0
        integer(isp) :: radix = 0
        complex(sp), dimension(:), allocatable :: twiddles_factor
    end type Twiddles_sp

    type :: FFTPlan
        integer(isp) :: N = 0
        integer(isp), dimension(:), allocatable :: radix_plan
        type(Twiddles_sp), dimension(:), allocatable :: twiddles
        logical :: is_initialized = .false.
        logical :: is_pure_radix2 = .false.
    end type FFTPlan

    type :: Fourier_Transform
        type(FFTPlan) :: fft_plan
    contains
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

        procedure, private :: init_fft_plan_sp
        generic :: init_fft_plan => init_fft_plan_sp

        procedure, private :: fft_cmplx_sp
        generic :: fft => fft_cmplx_sp

        procedure, private :: destroy_fft_plan_sp
        generic :: destroy_fft_plan => destroy_fft_plan_sp
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
    ! FFT plan creation and initialization
    !====================================================================================

    interface
        module subroutine init_fft_plan_sp(this, N)
            class(Fourier_Transform), intent(inout) :: this
            integer(isp), intent(in) :: N
        end subroutine init_fft_plan_sp
    end interface

    !====================================================================================
    ! compute FFT
    !====================================================================================

    interface
        module function fft_cmplx_sp(this, signal, loop_method) result(result)
            class(Fourier_Transform), intent(inout) :: this
            complex(sp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(:), allocatable :: result
        end function fft_cmplx_sp
    end interface

    !====================================================================================
    ! destr

    interface
        module subroutine destroy_fft_plan_sp(this)
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

    !====================================================================================
    ! FFT and IFFT for 1D, 2D, and 3D signals
    !====================================================================================

end module NAFPack_Fourier_Transform
