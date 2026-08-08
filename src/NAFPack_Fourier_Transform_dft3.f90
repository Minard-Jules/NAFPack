submodule(NAFPack_Fourier_Transform) NAFPack_Fourier_Transform_dft3

    implicit none(type, external)

contains

    module function dft3_cmplx_sp(signal, loop_method) result(result)
        complex(sp), dimension(:, :, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:, :, :), allocatable :: result
        complex(sp), dimension(:, :, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, j, Nx, Ny, Nz

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        Nz = size(signal, 3)
        allocate (result(Nx, Ny, Nz))
        allocate (tmp(Nx, Ny, Nz))

        do i = 1, Nx
            do j = 1, Ny
                tmp(i, j, :) = dft_cmplx_sp(signal(i, j, :), loop_method_used)
            end do
        end do

        do i = 1, Nx
            do j = 1, Nz
                result(i, :, j) = dft_cmplx_sp(tmp(i, :, j), loop_method_used)
            end do
        end do

        tmp = result
        do j = 1, Ny
            do i = 1, Nz
                result(:, j, i) = dft_cmplx_sp(tmp(:, j, i), loop_method_used)
            end do
        end do

    end function dft3_cmplx_sp

    module function idft3_cmplx_sp(f_signal, loop_method) result(result)
        complex(sp), dimension(:, :, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:, :, :), allocatable :: result
        complex(sp), dimension(:, :, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, j, Nx, Ny, Nz

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        Nz = size(f_signal, 3)
        allocate (result(Nx, Ny, Nz))
        allocate (tmp(Nx, Ny, Nz))

        do i = 1, Nx
            do j = 1, Ny
                tmp(i, j, :) = idft_cmplx_sp(f_signal(i, j, :), loop_method_used)
            end do
        end do

        do i = 1, Nx
            do j = 1, Nz
                result(i, :, j) = idft_cmplx_sp(tmp(i, :, j), loop_method_used)
            end do
        end do

        tmp = result
        do j = 1, Ny
            do i = 1, Nz
                result(:, j, i) = idft_cmplx_sp(tmp(:, j, i), loop_method_used)
            end do
        end do

    end function idft3_cmplx_sp

    module function dft3_cmplx_dp(signal, loop_method) result(result)
        complex(dp), dimension(:, :, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:, :, :), allocatable :: result
        complex(dp), dimension(:, :, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, j, Nx, Ny, Nz

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        Nz = size(signal, 3)
        allocate (result(Nx, Ny, Nz))
        allocate (tmp(Nx, Ny, Nz))

        do i = 1, Nx
            do j = 1, Ny
                tmp(i, j, :) = dft_cmplx_dp(signal(i, j, :), loop_method_used)
            end do
        end do

        do i = 1, Nx
            do j = 1, Nz
                result(i, :, j) = dft_cmplx_dp(tmp(i, :, j), loop_method_used)
            end do
        end do

        tmp = result
        do j = 1, Ny
            do i = 1, Nz
                result(:, j, i) = dft_cmplx_dp(tmp(:, j, i), loop_method_used)
            end do
        end do

    end function dft3_cmplx_dp

    module function idft3_cmplx_dp(f_signal, loop_method) result(result)
        complex(dp), dimension(:, :, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:, :, :), allocatable :: result
        complex(dp), dimension(:, :, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, j, Nx, Ny, Nz

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        Nz = size(f_signal, 3)
        allocate (result(Nx, Ny, Nz))
        allocate (tmp(Nx, Ny, Nz))

        do i = 1, Nx
            do j = 1, Ny
                tmp(i, j, :) = idft_cmplx_dp(f_signal(i, j, :), loop_method_used)
            end do
        end do

        do i = 1, Nx
            do j = 1, Nz
                result(i, :, j) = idft_cmplx_dp(tmp(i, :, j), loop_method_used)
            end do
        end do

        tmp = result
        do j = 1, Ny
            do i = 1, Nz
                result(:, j, i) = idft_cmplx_dp(tmp(:, j, i), loop_method_used)
            end do
        end do

    end function idft3_cmplx_dp

    module function dft3_cmplx_qp(signal, loop_method) result(result)
        complex(qp), dimension(:, :, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:, :, :), allocatable :: result
        complex(qp), dimension(:, :, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, j, Nx, Ny, Nz

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        Nz = size(signal, 3)
        allocate (result(Nx, Ny, Nz))
        allocate (tmp(Nx, Ny, Nz))

        do i = 1, Nx
            do j = 1, Ny
                tmp(i, j, :) = dft_cmplx_qp(signal(i, j, :), loop_method_used)
            end do
        end do

        do i = 1, Nx
            do j = 1, Nz
                result(i, :, j) = dft_cmplx_qp(tmp(i, :, j), loop_method_used)
            end do
        end do

        tmp = result
        do j = 1, Ny
            do i = 1, Nz
                result(:, j, i) = dft_cmplx_qp(tmp(:, j, i), loop_method_used)
            end do
        end do

    end function dft3_cmplx_qp

    module function idft3_cmplx_qp(f_signal, loop_method) result(result)
        complex(qp), dimension(:, :, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:, :, :), allocatable :: result
        complex(qp), dimension(:, :, :), allocatable :: tmp
        type(LoopMethod) :: loop_method_used
        integer(isp) :: i, j, Nx, Ny, Nz

        if (present(loop_method)) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        Nz = size(f_signal, 3)
        allocate (result(Nx, Ny, Nz))
        allocate (tmp(Nx, Ny, Nz))

        do i = 1, Nx
            do j = 1, Ny
                tmp(i, j, :) = idft_cmplx_qp(f_signal(i, j, :), loop_method_used)
            end do
        end do

        do i = 1, Nx
            do j = 1, Nz
                result(i, :, j) = idft_cmplx_qp(tmp(i, :, j), loop_method_used)
            end do
        end do

        tmp = result
        do j = 1, Ny
            do i = 1, Nz
                result(:, j, i) = idft_cmplx_qp(tmp(:, j, i), loop_method_used)
            end do
        end do

    end function idft3_cmplx_qp

    module function dft3_real_sp(signal, loop_method) result(result)
        real(sp), dimension(:, :, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:, :, :), allocatable :: result
        complex(sp), dimension(size(signal, 1), &
                               size(signal, 2), &
                               size(signal, 3)) :: signal_cmplx
        integer(isp) :: Nx, Ny, Nz

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        Nz = size(signal, 3)
        allocate (result(Nx, Ny, Nz))

        signal_cmplx = cmplx(signal, 0.0_dp, kind=sp)

        if (present(loop_method)) then
            result = dft3_cmplx_sp(signal_cmplx, loop_method)
        else
            result = dft3_cmplx_sp(signal_cmplx)
        end if
    end function dft3_real_sp

    module function idft3_real_sp(f_signal, loop_method) result(result)
        real(sp), dimension(:, :, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:, :, :), allocatable :: result
        complex(sp), dimension(size(f_signal, 1), &
                               size(f_signal, 2), &
                               size(f_signal, 3)) :: f_signal_cmplx
        integer(isp) :: Nx, Ny, Nz

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        Nz = size(f_signal, 3)
        allocate (result(Nx, Ny, Nz))

        f_signal_cmplx = cmplx(f_signal, 0.0_dp, kind=sp)

        if (present(loop_method)) then
            result = idft3_cmplx_sp(f_signal_cmplx, loop_method)
        else
            result = idft3_cmplx_sp(f_signal_cmplx)
        end if
    end function idft3_real_sp

    module function dft3_real_dp(signal, loop_method) result(result)
        real(dp), dimension(:, :, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:, :, :), allocatable :: result
        complex(dp), dimension(size(signal, 1), &
                               size(signal, 2), &
                               size(signal, 3)) :: signal_cmplx
        integer(isp) :: Nx, Ny, Nz

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        Nz = size(signal, 3)
        allocate (result(Nx, Ny, Nz))

        signal_cmplx = cmplx(signal, 0.0_dp, kind=dp)

        if (present(loop_method)) then
            result = dft3_cmplx_dp(signal_cmplx, loop_method)
        else
            result = dft3_cmplx_dp(signal_cmplx)
        end if
    end function dft3_real_dp

    module function idft3_real_dp(f_signal, loop_method) result(result)
        real(dp), dimension(:, :, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(:, :, :), allocatable :: result
        complex(dp), dimension(size(f_signal, 1), &
                               size(f_signal, 2), &
                               size(f_signal, 3)) :: f_signal_cmplx
        integer(isp) :: Nx, Ny, Nz

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        Nz = size(f_signal, 3)
        allocate (result(Nx, Ny, Nz))

        f_signal_cmplx = cmplx(f_signal, 0.0_dp, kind=dp)

        if (present(loop_method)) then
            result = idft3_cmplx_dp(f_signal_cmplx, loop_method)
        else
            result = idft3_cmplx_dp(f_signal_cmplx)
        end if
    end function idft3_real_dp

    module function dft3_real_qp(signal, loop_method) result(result)
        real(qp), dimension(:, :, :), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:, :, :), allocatable :: result
        complex(qp), dimension(size(signal, 1), &
                               size(signal, 2), &
                               size(signal, 3)) :: signal_cmplx
        integer(isp) :: Nx, Ny, Nz

        Nx = size(signal, 1)
        Ny = size(signal, 2)
        Nz = size(signal, 3)
        allocate (result(Nx, Ny, Nz))

        signal_cmplx = cmplx(signal, 0.0_dp, kind=qp)

        if (present(loop_method)) then
            result = dft3_cmplx_qp(signal_cmplx, loop_method)
        else
            result = dft3_cmplx_qp(signal_cmplx)
        end if
    end function dft3_real_qp

    module function idft3_real_qp(f_signal, loop_method) result(result)
        real(qp), dimension(:, :, :), intent(in) :: f_signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(:, :, :), allocatable :: result
        complex(qp), dimension(size(f_signal, 1), &
                               size(f_signal, 2), &
                               size(f_signal, 3)) :: f_signal_cmplx
        integer(isp) :: Nx, Ny, Nz

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)
        Nz = size(f_signal, 3)
        allocate (result(Nx, Ny, Nz))

        f_signal_cmplx = cmplx(f_signal, 0.0_dp, kind=qp)

        if (present(loop_method)) then
            result = idft3_cmplx_qp(f_signal_cmplx, loop_method)
        else
            result = idft3_cmplx_qp(f_signal_cmplx)
        end if
    end function idft3_real_qp

end submodule NAFPack_Fourier_Transform_dft3
