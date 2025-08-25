!> Module for Fourier Transform
!> \[ F(\omega) = \int_{-\infty}^{\infty} f(t) e^{-i \omega t} dt \]
!> This module provides an interface for performing Fourier Transforms (FFT or DFT, IFFT) on 1D, 2D, and 3D signals.
!> It supports both forward and inverse transforms.
!> It allows users to choose between different methods for the Fourier Transform, such as NAFPack and FFTW.
module NAFPack_fft

    use FFTW3

    use NAFPack_constant
    implicit none(type, external)

    private
    public :: FFT_1D, FFT_2D, FFT_3D
    public :: IFFT_1D, IFFT_2D, IFFT_3D

contains

    !> Perform a 1D Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "NAFPack_DFT": Direct Discrete Fourier Transform
    !> - "NAFPack_FFT_1D": Fast Fourier Transform using NAFPack
    !> - "FFTW_FFT_1D": Fast Fourier Transform using FFTW
    !> - "FFTW_FFT_1D" + threads: Fast Fourier Transform using FFTW with multithreading
    function FFT_1D(signal, method, threads) result(result)

        complex(dp), dimension(:), intent(inout) :: signal
        character(*), intent(in) :: method
        integer, optional, intent(in) :: threads
        complex(dp), dimension(size(signal)) :: result

        if (method == "NAFPack_DFT") then
            result = NAFPack_DFT_1D(signal)
        else if (method == "NAFPack_FFT_1D") then
            result = NAFPack_FFT_1D(signal)
        else if (method == "FFTW_FFT_1D" .and. .not. present(threads)) then
            result = FFTW_FFT_1D(signal)
        else if (method == "FFTW_FFT_1D" .and. present(threads)) then
            result = FFTW_FFT_1D_threads(signal, threads)
        else
            stop "ERROR : Wrong method for FFT_1D"
        end if

    end function FFT_1D

    !> Perform a 1D inverse Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a  inverse fast Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "NAFPack_IFFT_1D": Fast Fourier Transform using NAFPack
    !> - "FFTW_IFFT_1D": Fast Fourier Transform using FFTW
    !> - "FFTW_IFFT_1D" + threads: Fast Fourier Transform using FFTW with multithreading
    function IFFT_1D(signal, method, threads) result(result)

        complex(dp), dimension(:), intent(inout) :: signal
        character(*), intent(in) :: method
        integer, optional, intent(in) :: threads
        complex(dp), dimension(size(signal)) :: result

        if (method == "NAFPack_IFFT_1D") then
            result = NAFPack_IFFT_1D(signal)
        else if (method == "FFTW_IFFT_1D" .and. .not. present(threads)) then
            result = FFTW_IFFT_1D(signal)
        else if (method == "FFTW_IFFT_1D" .and. present(threads)) then
            result = FFTW_IFFT_1D_threads(signal, threads)
        else
            stop "ERROR : Wrong method for IFFT_1D"
        end if

    end function IFFT_1D

    !> Perform a 2D Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "NAFPack_FFT_2D": Fast Fourier Transform using NAFPack
    !> - "FFTW_FFT_2D": Fast Fourier Transform using FFTW
    !> - "FFTW_FFT_2D" + threads: Fast Fourier Transform using FFTW with multithreading
    function FFT_2D(signal, method, threads) result(result)

        complex(dp), dimension(:, :), intent(inout) :: signal
        character(*), intent(in) :: method
        integer, optional, intent(in) :: threads
        complex(dp), dimension(size(signal, 1), size(signal, 2)) :: result

        if (method == "NAFPack_FFT_2D") then
            result = NAFPack_FFT_2D(signal)
        else if (method == "FFTW_FFT_2D" .and. .not. present(threads)) then
            result = FFTW_FFT_2D(signal)
        else if (method == "FFTW_FFT_2D" .and. present(threads)) then
            result = FFTW_FFT_2D_threads(signal, threads)
        else
            stop "ERROR : Wrong method for FFT_2D"
        end if

    end function FFT_2D

    !> Perform a 2D inverse Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a  inverse fast Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "NAFPack_IFFT_2D": Fast Fourier Transform using NAFPack
    !> - "FFTW_IFFT_2D": Fast Fourier Transform using FFTW
    !> - "FFTW_IFFT_2D" + threads: Fast Fourier Transform using FFTW with multithreading
    function IFFT_2D(signal, method, threads) result(result)

        complex(dp), dimension(:, :), intent(inout) :: signal
        character(*), intent(in) :: method
        integer, optional, intent(in) :: threads
        complex(dp), dimension(size(signal, 1), size(signal, 2)) :: result

        if (method == "NAFPack_IFFT_2D") then
            result = NAFPack_IFFT_2D(signal)
        else if (method == "FFTW_IFFT_2D" .and. .not. present(threads)) then
            result = FFTW_IFFT_2D(signal)
        else if (method == "FFTW_IFFT_2D" .and. present(threads)) then
            result = FFTW_IFFT_2D_threads(signal, threads)
        else
            stop "ERROR : Wrong method for IFFT_1D"
        end if

    end function IFFT_2D

    !> Perform a 3D Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "FFTW_FFT_3D": Fast Fourier Transform using FFTW
    !> - "FFTW_FFT_3D" + threads: Fast Fourier Transform using FFTW with multithreading
    function FFT_3D(signal, method, threads) result(result)

        complex(dp), dimension(:, :, :), intent(inout) :: signal
        character(*), intent(in) :: method
        integer, optional, intent(in) :: threads
        complex(dp), dimension(size(signal, 1), size(signal, 2), size(signal, 3)) :: result

        if (method == "FFTW_FFT_3D" .and. .not. present(threads)) then
            result = FFTW_FFT_3D(signal)
        else if (method == "FFTW_FFT_3D" .and. present(threads)) then
            result = FFTW_FFT_3D_threads(signal, threads)
        else
            stop "ERROR : Wrong method for FFT_2D"
        end if

    end function FFT_3D

    !> Perform a 3D inverse Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a  inverse fast Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "FFTW_IFFT_3D": Fast Fourier Transform using FFTW
    !> - "FFTW_IFFT_3D" + threads: Fast Fourier Transform using FFTW with multithreading
    function IFFT_3D(signal, method, threads) result(result)

        complex(dp), dimension(:, :, :), intent(inout) :: signal
        character(*), intent(in) :: method
        integer, optional, intent(in) :: threads
        complex(dp), dimension(size(signal, 1), size(signal, 2), size(signal, 3)) :: result

        if (method == "FFTW_IFFT_3D" .and. .not. present(threads)) then
            result = FFTW_IFFT_3D(signal)
        else if (method == "IFFTW_IFFT_3D" .and. present(threads)) then
            result = FFTW_IFFT_3D_threads(signal, threads)
        else
            stop "ERROR : Wrong method for IFFT_1D"
        end if

    end function IFFT_3D

!################### FFTW ##########################################

    !> Perform a 1D Fast Fourier Transform on a signal using FFTW with multithreading
    function FFTW_FFT_1D_threads(signal, threads) result(result)

        complex(c_double_complex), dimension(:), intent(inout) :: signal
        integer, intent(in) :: threads
        complex(c_double_complex), dimension(size(signal)) :: result
        integer :: error_init_thread
        type(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        if (error_init_thread == 0) stop "ERROR : Thread FFTW initialization problem"

        call fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_1d(size(signal), signal, result, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        call fftw_destroy_plan(plan)

        call fftw_cleanup_threads()

    end function FFTW_FFT_1D_threads

    !> Perform a 1D Fast Fourier Transform on a signal using FFTW
    function FFTW_FFT_1D(signal) result(result)

        complex(c_double_complex), dimension(:), intent(inout) :: signal
        complex(c_double_complex), dimension(size(signal)) :: result
        type(c_ptr) :: plan

        plan = fftw_plan_dft_1d(size(signal), signal, result, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        call fftw_destroy_plan(plan)

    end function FFTW_FFT_1D

    !> Perform a 1D inverse Fast Fourier Transform on a signal using FFTW
    function FFTW_IFFT_1D(signal) result(result)

        complex(c_double_complex), dimension(:), intent(inout) :: signal
        complex(c_double_complex), dimension(size(signal)) :: result
        type(c_ptr) :: plan

        plan = fftw_plan_dft_1d(size(signal), signal, result, FFTW_BACKWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        result = result / size(signal)
        call fftw_destroy_plan(plan)

    end function FFTW_IFFT_1D

    !> Perform a 1D inverse Fast Fourier Transform on a signal using FFTW with multithreading
    function FFTW_IFFT_1D_threads(signal, threads) result(result)

        complex(c_double_complex), dimension(:), intent(inout) :: signal
        integer, intent(in) :: threads
        complex(c_double_complex), dimension(size(signal)) :: result
        integer :: error_init_thread
        type(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        if (error_init_thread == 0) stop "ERROR : Thread FFTW initialization problem"

        call fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_1d(size(signal), signal, result, FFTW_BACKWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        result = result / size(signal)
        call fftw_destroy_plan(plan)

        call fftw_cleanup_threads()

    end function FFTW_IFFT_1D_threads

    !> Perform a 2D Fast Fourier Transform on a signal using FFTW
    function FFTW_FFT_2D(signal) result(result)

        complex(c_double_complex), dimension(:, :), intent(inout) :: signal
        complex(c_double_complex), dimension(size(signal, 1), size(signal, 2)) :: result
        type(c_ptr) :: plan

        plan = fftw_plan_dft_2d(size(signal, 2), size(signal, 1), signal, result, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        call fftw_destroy_plan(plan)

    end function FFTW_FFT_2D

    !> Perform a 2D Fast Fourier Transform on a signal using FFTW with multithreading
    function FFTW_FFT_2D_threads(signal, threads) result(result)

        complex(c_double_complex), dimension(:, :), intent(inout) :: signal
        integer, intent(in) :: threads
        complex(c_double_complex), dimension(size(signal, 1), size(signal, 2)) :: result
        integer :: error_init_thread
        type(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        if (error_init_thread == 0) stop "ERROR : Thread FFTW initialization problem"

        call fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_2d(size(signal, 2), size(signal, 1), signal, result, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        call fftw_destroy_plan(plan)

        call fftw_cleanup_threads()

    end function FFTW_FFT_2D_threads

    !> Perform a 2D inverse Fast Fourier Transform on a signal using FFTW
    function FFTW_IFFT_2D(signal) result(result)

        complex(c_double_complex), dimension(:, :), intent(inout) :: signal
        complex(c_double_complex), dimension(size(signal, 1), size(signal, 2)) :: result
        type(c_ptr) :: plan

        plan = fftw_plan_dft_2d(size(signal, 2), size(signal, 1), signal, result, FFTW_BACKWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        result = result / (size(signal, 1) * size(signal, 2))
        call fftw_destroy_plan(plan)

    end function FFTW_IFFT_2D

    !> Perform a 2D inverse Fast Fourier Transform on a signal using FFTW with multithreading
    function FFTW_IFFT_2D_threads(signal, threads) result(result)

        complex(c_double_complex), dimension(:, :), intent(inout) :: signal
        integer, intent(in) :: threads
        complex(c_double_complex), dimension(size(signal, 1), size(signal, 2)) :: result
        integer :: error_init_thread
        type(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        if (error_init_thread == 0) stop "ERROR : Thread FFTW initialization problem"

        call fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_2d(size(signal, 2), size(signal, 1), signal, result, FFTW_BACKWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        result = result / (size(signal, 1) * size(signal, 2))
        call fftw_destroy_plan(plan)

        call fftw_cleanup_threads()

    end function FFTW_IFFT_2D_threads

    !> Perform a 3D Fast Fourier Transform on a signal using FFTW
    function FFTW_FFT_3D(signal) result(result)

        complex(c_double_complex), dimension(:, :, :), intent(inout) :: signal
        complex(c_double_complex), dimension(size(signal, 1), size(signal, 2), size(signal, 3)) :: result
        type(c_ptr) :: plan

        plan = fftw_plan_dft_3d(size(signal, 3), size(signal, 2), size(signal, 1), signal, result, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        call fftw_destroy_plan(plan)

    end function FFTW_FFT_3D

    !> Perform a 3D Fast Fourier Transform on a signal using FFTW with multithreading
    function FFTW_FFT_3D_threads(signal, threads) result(result)

        complex(c_double_complex), dimension(:, :, :), intent(inout) :: signal
        integer, intent(in) :: threads
        complex(c_double_complex), dimension(size(signal, 1), size(signal, 2), size(signal, 3)) :: result
        integer :: error_init_thread
        type(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        if (error_init_thread == 0) stop "ERROR : Thread FFTW initialization problem"

        call fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_3d(size(signal, 3), size(signal, 2), size(signal, 1), signal, result, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        call fftw_destroy_plan(plan)

        call fftw_cleanup_threads()

    end function FFTW_FFT_3D_threads

    !> Perform a 3D inverse Fast Fourier Transform on a signal using FFTW
    function FFTW_IFFT_3D(signal) result(result)

        complex(c_double_complex), dimension(:, :, :), intent(inout) :: signal
        complex(c_double_complex), dimension(size(signal, 1), size(signal, 2), size(signal, 3)) :: result
        type(c_ptr) :: plan

        plan = fftw_plan_dft_3d(size(signal, 3), size(signal, 2), size(signal, 1), signal, result, FFTW_BACKWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        result = result / (size(signal, 1) * size(signal, 2) * size(signal, 3))
        call fftw_destroy_plan(plan)

    end function FFTW_IFFT_3D

    !> Perform a 3D inverse Fast Fourier Transform on a signal using FFTW with multithreading
    function FFTW_IFFT_3D_threads(signal, threads) result(result)

        complex(c_double_complex), dimension(:, :, :), intent(inout) :: signal
        integer, intent(in) :: threads
        complex(c_double_complex), dimension(size(signal, 1), size(signal, 2), size(signal, 3)) :: result
        integer :: error_init_thread
        type(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        if (error_init_thread == 0) stop "ERROR : Thread FFTW initialization problem"

        call fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_3d(size(signal, 3), size(signal, 2), size(signal, 1), signal, result, FFTW_BACKWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, signal, result)
        result = result / (size(signal, 1) * size(signal, 2) * size(signal, 3))
        call fftw_destroy_plan(plan)

        call fftw_cleanup_threads()

    end function FFTW_IFFT_3D_threads

!################### NAFPack ##########################################

    !> Perform a 1D Discrete Fourier Transform on a signal
    function NAFPack_DFT_1D(signal) result(result)

        complex(dp), dimension(:), intent(in) :: signal
        complex(dp), dimension(size(signal)) :: result
        complex(dp) :: S
        integer :: N, i, k, j

        N = size(signal)

        if (N == 1) then
            result = signal
        else
            do i = 1, N

                k = i - 1
                S = (0.d0, 0.d0)

                do j = 1, N
                    S = S + signal(j) * exp((-2 * pi * im * k * (j - 1)) / N)
                end do

                result(i) = S

            end do
        end if

    end function NAFPack_DFT_1D

    !> Compute the complex exponential factors for the FFT
    function fun_omega(N) result(result)

        integer, intent(in) :: N
        complex(dp), dimension(N/2) :: result
        integer :: i

        do i = 1, N / 2
            result(i) = exp(-2 * Im * pi * (i - 1) / N)
        end do

    end function fun_omega

    !> Perform a 1D Fast Fourier Transform (Cooley-Tukey) on a signal
    recursive function NAFPack_FFT_1D(signal) result(result)

        complex(dp), dimension(:), intent(in) :: signal
        complex(dp), dimension(size(signal)) :: result
        complex(dp), dimension(size(signal)/2) :: f_pair, f_impair, omega
        integer :: N

        N = size(signal)

        if (mod(N, 2) == 0) then
            f_pair = NAFPack_FFT_1D(signal(1 :: 2))
            f_impair = NAFPack_FFT_1D(signal(2 :: 2))

            omega = fun_omega(N)

            result(1:N / 2) = f_pair + f_impair * omega
            result(N / 2 + 1:) = f_pair - f_impair * omega
        else
            result = NAFPack_DFT_1D(signal)
        end if
    end function NAFPack_FFT_1D

    !> Perform a 1D inverse Fast Fourier Transform on a signal
    function NAFPack_IFFT_1D(f_signal) result(result)

        complex(dp), dimension(:), intent(in) :: f_signal
        complex(dp), dimension(size(f_signal)) :: result
        complex(dp), dimension(size(f_signal)) :: f_conjugate
        integer :: N

        N = size(f_signal)

        f_conjugate = conjg(f_signal)

        result = NAFPack_FFT_1D(f_conjugate)
        result = conjg(result)
        result = result / N

    end function NAFPack_IFFT_1D

    !> Perform a 2D Fast Fourier Transform on a signal
    function NAFPack_FFT_2D(signal) result(result)

        complex(dp), dimension(:, :), intent(in) :: signal
        complex(dp), dimension(size(signal, 1), size(signal, 2)) :: result
        integer :: Nx, Ny, i

        Nx = size(signal, 1)
        Ny = size(signal, 2)

        do i = 1, Nx
            result(i, :) = NAFPack_FFT_1D(signal(i, :))
        end do

        do i = 1, Ny
            result(:, i) = NAFPack_FFT_1D(result(:, i))
        end do

    end function NAFPack_FFT_2D

    !> Perform a 2D inverse Fast Fourier Transform on a signal
    function NAFPack_IFFT_2D(f_signal) result(result)

        complex(dp), dimension(:, :), intent(in) :: f_signal
        complex(dp), dimension(size(f_signal, 1), size(f_signal, 2)) :: result
        integer :: Nx, Ny, i

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)

        do i = 1, Nx
            result(i, :) = NAFPack_IFFT_1D(f_signal(i, :))
        end do

        do i = 1, Ny
            result(:, i) = NAFPack_IFFT_1D(result(:, i))
        end do

    end function NAFPack_IFFT_2D

end module NAFPack_fft
