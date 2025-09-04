module NAFPack_Fourier_Transform

    use NAFPack_kinds, only: dp, sp, qp, i8, i16, isp, idp
    use NAFPack_constant, only: pi_sp, pi_dp, pi_qp, im_sp, im_dp, im_qp
    use NAFPack_loop_method, only: LoopMethod, default_loop_method, check_loop_method

    implicit none(type, external)

    private
    public :: Fourier_Transform
    public :: dft

    type :: Fourier_Transform
    contains
        procedure, nopass, private :: dft_cmplx_sp, dft_cmplx_dp, dft_cmplx_qp
        procedure, nopass, private :: dft_real_sp, dft_real_dp, dft_real_qp
        generic :: dft => &
            dft_cmplx_sp, dft_cmplx_dp, dft_cmplx_qp, &
            dft_real_sp, dft_real_dp, dft_real_qp
    end type Fourier_Transform

    interface dft
        module procedure dft_cmplx_sp, dft_cmplx_dp, dft_cmplx_qp
        module procedure dft_real_sp, dft_real_dp, dft_real_qp
    end interface dft

    interface
        module function dft_cmplx_sp(signal, loop_method) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(size(signal)) :: result
        end function dft_cmplx_sp

        module function dft_cmplx_dp(signal, loop_method) result(result)
            complex(dp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(size(signal)) :: result
        end function dft_cmplx_dp

        module function dft_cmplx_qp(signal, loop_method) result(result)
            complex(qp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(size(signal)) :: result
        end function dft_cmplx_qp
    end interface

    interface
        module function dft_real_sp(signal, loop_method) result(result)
            real(sp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(sp), dimension(size(signal)) :: result
        end function dft_real_sp

        module function dft_real_dp(signal, loop_method) result(result)
            real(dp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(dp), dimension(size(signal)) :: result
        end function dft_real_dp

        module function dft_real_qp(signal, loop_method) result(result)
            real(qp), dimension(:), intent(in) :: signal
            type(LoopMethod), optional, intent(in) :: loop_method
            complex(qp), dimension(size(signal)) :: result
        end function dft_real_qp
    end interface
end module NAFPack_Fourier_Transform
