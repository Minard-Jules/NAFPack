submodule(NAFPack_Fourier_Transform) NAFPack_Fourier_Transform_dft

    implicit none(type, external)

    interface
        module function compute_dft_cmplx_sp(signal, N, loop_method) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            integer(isp), intent(in) :: N
            type(LoopMethod), intent(in) :: loop_method
            complex(sp), dimension(N) :: result
        end function compute_dft_cmplx_sp
    end interface

    interface
        module function compute_dft_cmplx_dp(signal, N, loop_method) result(result)
            complex(dp), dimension(:), intent(in) :: signal
            integer(isp), intent(in) :: N
            type(LoopMethod), intent(in) :: loop_method
            complex(dp), dimension(N) :: result
        end function compute_dft_cmplx_dp
    end interface

    interface
        module function compute_dft_cmplx_qp(signal, N, loop_method) result(result)
            complex(qp), dimension(:), intent(in) :: signal
            integer(isp), intent(in) :: N
            type(LoopMethod), intent(in) :: loop_method
            complex(qp), dimension(N) :: result
        end function compute_dft_cmplx_qp
    end interface 

contains

    module function dft_cmplx_sp(signal, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:), allocatable :: result
        type(LoopMethod) :: loop_method_used
        integer(isp) :: N

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        N = size(signal)
        allocate (result(N))

        result = compute_dft_cmplx_sp(signal, N, loop_method_used)

    end function dft_cmplx_sp

    module function idft_cmplx_sp(f_signal, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:), allocatable :: result
        integer(isp) :: N

        N = size(f_signal)
        allocate (result(N))
        result = conjg(dft_cmplx_sp(conjg(f_signal), loop_method)) / N
    end function idft_cmplx_sp

    module function dft_cmplx_dp(signal, loop_method) result(result)
        complex(dp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:), allocatable :: result
        type(LoopMethod) :: loop_method_used
        integer(isp) :: N

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        N = size(signal)
        allocate (result(N))

        result = compute_dft_cmplx_dp(signal, N, loop_method_used)

    end function dft_cmplx_dp

    module function idft_cmplx_dp(f_signal, loop_method) result(result)
        complex(dp), dimension(:), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:), allocatable :: result
        integer(isp) :: N

        N = size(f_signal)
        allocate (result(N))
        result = conjg(dft_cmplx_dp(conjg(f_signal), loop_method)) / N
    end function idft_cmplx_dp

    module function dft_cmplx_qp(signal, loop_method) result(result)
        complex(qp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:), allocatable :: result
        type(LoopMethod) :: loop_method_used
        integer(isp) :: N

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        N = size(signal)
        allocate (result(N))

       result = compute_dft_cmplx_qp(signal, N, loop_method_used)

    end function dft_cmplx_qp

    module function idft_cmplx_qp(f_signal, loop_method) result(result)
        complex(qp), dimension(:), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:), allocatable :: result
        integer(isp) :: N

        N = size(f_signal)
        allocate (result(N))
        result = conjg(dft_cmplx_qp(conjg(f_signal), loop_method)) / N
    end function idft_cmplx_qp

    module function dft_real_sp(signal, loop_method) result(result)
        real(sp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:), allocatable :: result
        complex(sp), dimension(size(signal)) :: signal_cmplx
        integer(isp) :: N

        N = size(signal)
        allocate (result(N))

        signal_cmplx = cmplx(signal, 0.0_dp, kind=sp)

        if (present(loop_method)) then
            result = dft_cmplx_sp(signal_cmplx, loop_method)
        else
            result = dft_cmplx_sp(signal_cmplx)
        end if
    end function dft_real_sp

    module function idft_real_sp(f_signal, loop_method) result(result)
        real(sp), dimension(:), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:), allocatable :: result
        complex(sp), dimension(size(f_signal)) :: f_signal_cmplx
        integer(isp) :: N

        N = size(f_signal)
        allocate (result(N))

        f_signal_cmplx = cmplx(f_signal, 0.0_sp, kind=sp)
        result = conjg(dft_cmplx_sp(conjg(f_signal_cmplx), loop_method)) / N
    end function idft_real_sp

    module function dft_real_dp(signal, loop_method) result(result)
        real(dp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:), allocatable :: result
        complex(dp), dimension(size(signal)) :: signal_cmplx
        integer(isp) :: N

        N = size(signal)
        allocate (result(N))

        signal_cmplx = cmplx(signal, 0.0_dp, kind=dp)

        if (present(loop_method)) then
            result = dft_cmplx_dp(signal_cmplx, loop_method)
        else
            result = dft_cmplx_dp(signal_cmplx)
        end if
    end function dft_real_dp

    module function idft_real_dp(f_signal, loop_method) result(result)
        real(dp), dimension(:), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:), allocatable :: result
        complex(dp), dimension(size(f_signal)) :: f_signal_cmplx
        integer(isp) :: N

        N = size(f_signal)
        allocate (result(N))

        f_signal_cmplx = cmplx(f_signal, 0.0_dp, kind=dp)
        result = conjg(dft_cmplx_dp(conjg(f_signal_cmplx), loop_method)) / N
    end function idft_real_dp

    module function dft_real_qp(signal, loop_method) result(result)
        real(qp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:), allocatable :: result
        complex(qp), dimension(size(signal)) :: signal_cmplx
        integer(isp) :: N

        N = size(signal)
        allocate (result(N))

        signal_cmplx = cmplx(signal, 0.0_qp, kind=qp)

        if (present(loop_method)) then
            result = dft_cmplx_qp(signal_cmplx, loop_method)
        else
            result = dft_cmplx_qp(signal_cmplx)
        end if
    end function dft_real_qp

    module function idft_real_qp(f_signal, loop_method) result(result)
        real(qp), dimension(:), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:), allocatable :: result
        complex(qp), dimension(size(f_signal)) :: f_signal_cmplx
        integer(isp) :: N

        N = size(f_signal)
        allocate (result(N))

        f_signal_cmplx = cmplx(f_signal, 0.0_qp, kind=qp)
        result = conjg(dft_cmplx_qp(conjg(f_signal_cmplx), loop_method)) / N
    end function idft_real_qp

end submodule NAFPack_Fourier_Transform_dft
