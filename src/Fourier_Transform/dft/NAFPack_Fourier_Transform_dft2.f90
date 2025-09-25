submodule(NAFPack_Fourier_Transform) NAFPack_Fourier_Transform_dft2

    implicit none(type, external)

contains

    module function dft2_cmplx_sp(signal, loop_method) result(result)
        complex(sp), dimension(:, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:, :), allocatable :: result
        complex(sp), dimension(:, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, Nx, Ny

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        allocate (result(Nx, Ny))
        allocate (tmp(Nx, Ny))

        do i = 1, Ny
            tmp(:, i) = dft_cmplx_sp(signal(:, i), loop_method)
        end do

        do i = 1, Nx
            result(i, :) = dft_cmplx_sp(tmp(i, :), loop_method)
        end do

    end function dft2_cmplx_sp

    module function idft2_cmplx_sp(f_signal, loop_method) result(result)
        complex(sp), dimension(:, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:, :), allocatable :: result
        complex(sp), dimension(:, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, Nx, Ny

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        allocate (result(Nx, Ny))
        allocate (tmp(Nx, Ny))

        do i = 1, Ny
            tmp(:, i) = idft_cmplx_sp(f_signal(:, i), loop_method)
        end do

        do i = 1, Nx
            result(i, :) = idft_cmplx_sp(tmp(i, :), loop_method)
        end do

    end function idft2_cmplx_sp

    module function dft2_cmplx_dp(signal, loop_method) result(result)
        complex(dp), dimension(:, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:, :), allocatable :: result
        complex(dp), dimension(:, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, Nx, Ny

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        allocate (result(Nx, Ny))
        allocate (tmp(Nx, Ny))

        do i = 1, Ny
            tmp(:, i) = dft_cmplx_dp(signal(:, i), loop_method)
        end do

        do i = 1, Nx
            result(i, :) = dft_cmplx_dp(tmp(i, :), loop_method)
        end do

    end function dft2_cmplx_dp

    module function idft2_cmplx_dp(f_signal, loop_method) result(result)
        complex(dp), dimension(:, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:, :), allocatable :: result
        complex(dp), dimension(:, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, Nx, Ny

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        allocate (result(Nx, Ny))
        allocate (tmp(Nx, Ny))

        do i = 1, Ny
            tmp(:, i) = idft_cmplx_dp(f_signal(:, i), loop_method)
        end do

        do i = 1, Nx
            result(i, :) = idft_cmplx_dp(tmp(i, :), loop_method)
        end do

    end function idft2_cmplx_dp

    module function dft2_cmplx_qp(signal, loop_method) result(result)
        complex(qp), dimension(:, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:, :), allocatable :: result
        complex(qp), dimension(:, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, Nx, Ny

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        allocate (result(Nx, Ny))
        allocate (tmp(Nx, Ny))

        do i = 1, Ny
            tmp(:, i) = dft_cmplx_qp(signal(:, i), loop_method)
        end do

        do i = 1, Nx
            result(i, :) = dft_cmplx_qp(tmp(i, :), loop_method)
        end do

    end function dft2_cmplx_qp

    module function idft2_cmplx_qp(f_signal, loop_method) result(result)
        complex(qp), dimension(:, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:, :), allocatable :: result
        complex(qp), dimension(:, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, Nx, Ny

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        allocate (result(Nx, Ny))
        allocate (tmp(Nx, Ny))

        do i = 1, Ny
            tmp(:, i) = idft_cmplx_qp(f_signal(:, i), loop_method)
        end do

        do i = 1, Nx
            result(i, :) = idft_cmplx_qp(tmp(i, :), loop_method)
        end do

    end function idft2_cmplx_qp

    module function dft2_real_sp(signal, loop_method) result(result)
        real(sp), dimension(:, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:, :), allocatable :: result
        complex(sp), dimension(size(signal, 1), &
                               size(signal, 2)) :: signal_cmplx
        integer(isp) :: Nx, Ny

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        allocate (result(Nx, Ny))

        signal_cmplx = cmplx(signal, 0.0_dp, kind=sp)

        if (present(loop_method)) then
            result = dft2_cmplx_sp(signal_cmplx, loop_method)
        else
            result = dft2_cmplx_sp(signal_cmplx)
        end if
    end function dft2_real_sp

    module function idft2_real_sp(f_signal, loop_method) result(result)
        real(sp), dimension(:, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:, :), allocatable :: result
        complex(sp), dimension(size(f_signal, 1), &
                               size(f_signal, 2)) :: f_signal_cmplx
        integer(isp) :: Nx, Ny

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        allocate (result(Nx, Ny))

        f_signal_cmplx = cmplx(f_signal, 0.0_dp, kind=sp)

        if (present(loop_method)) then
            result = idft2_cmplx_sp(f_signal_cmplx, loop_method)
        else
            result = idft2_cmplx_sp(f_signal_cmplx)
        end if
    end function idft2_real_sp

    module function dft2_real_dp(signal, loop_method) result(result)
        real(dp), dimension(:, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:, :), allocatable :: result
        complex(dp), dimension(size(signal, 1), &
                               size(signal, 2)) :: signal_cmplx
        integer(isp) :: Nx, Ny

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        allocate (result(Nx, Ny))

        signal_cmplx = cmplx(signal, 0.0_dp, kind=dp)

        if (present(loop_method)) then
            result = dft2_cmplx_dp(signal_cmplx, loop_method)
        else
            result = dft2_cmplx_dp(signal_cmplx)
        end if
    end function dft2_real_dp

    module function idft2_real_dp(f_signal, loop_method) result(result)
        real(dp), dimension(:, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:, :), allocatable :: result
        complex(dp), dimension(size(f_signal, 1), &
                               size(f_signal, 2)) :: f_signal_cmplx
        integer(isp) :: Nx, Ny

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        allocate (result(Nx, Ny))

        f_signal_cmplx = cmplx(f_signal, 0.0_dp, kind=dp)

        if (present(loop_method)) then
            result = idft2_cmplx_dp(f_signal_cmplx, loop_method)
        else
            result = idft2_cmplx_dp(f_signal_cmplx)
        end if
    end function idft2_real_dp

    module function dft2_real_qp(signal, loop_method) result(result)
        real(qp), dimension(:, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:, :), allocatable :: result
        complex(qp), dimension(size(signal, 1), &
                               size(signal, 2)) :: signal_cmplx
        integer(isp) :: Nx, Ny

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        allocate (result(Nx, Ny))

        signal_cmplx = cmplx(signal, 0.0_dp, kind=qp)

        if (present(loop_method)) then
            result = dft2_cmplx_qp(signal_cmplx, loop_method)
        else
            result = dft2_cmplx_qp(signal_cmplx)
        end if
    end function dft2_real_qp

    module function idft2_real_qp(f_signal, loop_method) result(result)
        real(qp), dimension(:, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:, :), allocatable :: result
        complex(qp), dimension(size(f_signal, 1), &
                               size(f_signal, 2)) :: f_signal_cmplx
        integer(isp) :: Nx, Ny

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        allocate (result(Nx, Ny))

        f_signal_cmplx = cmplx(f_signal, 0.0_dp, kind=qp)

        if (present(loop_method)) then
            result = idft2_cmplx_qp(f_signal_cmplx, loop_method)
        else
            result = idft2_cmplx_qp(f_signal_cmplx)
        end if
    end function idft2_real_qp

end submodule NAFPack_Fourier_Transform_dft2
