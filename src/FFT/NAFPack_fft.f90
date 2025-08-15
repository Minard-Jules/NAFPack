!> Module for Fourier Transform
!> \[ F(\omega) = \int_{-\infty}^{\infty} f(t) e^{-i \omega t} dt \]
!> This module provides an interface for performing Fourier Transforms (FFT or DFT, IFFT) on 1D, 2D, and 3D signals.
!> It supports both forward and inverse transforms.
!> It allows users to choose between different methods for the Fourier Transform, such as NAFPack and FFTW.
MODULE NAFPack_fft

    USE FFTW3

    USE NAFPack_constant
    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE
    PUBLIC :: FFT_1D, FFT_2D, FFT_3D
    PUBLIC :: IFFT_1D, IFFT_2D, IFFT_3D

CONTAINS

    !> Perform a 1D Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "NAFPack_DFT": Direct Discrete Fourier Transform
    !> - "NAFPack_FFT_1D": Fast Fourier Transform using NAFPack
    !> - "FFTW_FFT_1D": Fast Fourier Transform using FFTW
    !> - "FFTW_FFT_1D" + threads: Fast Fourier Transform using FFTW with multithreading
    FUNCTION FFT_1D(signal, method, threads) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:), INTENT(INOUT) :: signal
        CHARACTER(*), INTENT(IN) :: method
        INTEGER, OPTIONAL, INTENT(IN) :: threads
        COMPLEX(dp), DIMENSION(size(signal)) :: RESULT

        IF (method == "NAFPack_DFT") THEN
            RESULT = NAFPack_DFT_1D(signal)
        ELSE IF (method == "NAFPack_FFT_1D") THEN
            RESULT = NAFPack_FFT_1D(signal)
        ELSE IF (method == "FFTW_FFT_1D" .AND. .NOT. present(threads)) THEN
            RESULT = FFTW_FFT_1D(signal)
        ELSE IF (method == "FFTW_FFT_1D" .AND. present(threads)) THEN
            RESULT = FFTW_FFT_1D_threads(signal, threads)
        ELSE
            STOP "ERROR : Wrong method for FFT_1D"
        END IF

    END FUNCTION FFT_1D

    !> Perform a 1D inverse Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a  inverse fast Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "NAFPack_IFFT_1D": Fast Fourier Transform using NAFPack
    !> - "FFTW_IFFT_1D": Fast Fourier Transform using FFTW
    !> - "FFTW_IFFT_1D" + threads: Fast Fourier Transform using FFTW with multithreading
    FUNCTION IFFT_1D(signal, method, threads) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:), INTENT(INOUT) :: signal
        CHARACTER(*), INTENT(IN) :: method
        INTEGER, OPTIONAL, INTENT(IN) :: threads
        COMPLEX(dp), DIMENSION(size(signal)) :: RESULT

        IF (method == "NAFPack_IFFT_1D") THEN
            RESULT = NAFPack_IFFT_1D(signal)
        ELSE IF (method == "FFTW_IFFT_1D" .AND. .NOT. present(threads)) THEN
            RESULT = FFTW_IFFT_1D(signal)
        ELSE IF (method == "FFTW_IFFT_1D" .AND. present(threads)) THEN
            RESULT = FFTW_IFFT_1D_threads(signal, threads)
        ELSE
            STOP "ERROR : Wrong method for IFFT_1D"
        END IF

    END FUNCTION IFFT_1D

    !> Perform a 2D Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "NAFPack_FFT_2D": Fast Fourier Transform using NAFPack
    !> - "FFTW_FFT_2D": Fast Fourier Transform using FFTW
    !> - "FFTW_FFT_2D" + threads: Fast Fourier Transform using FFTW with multithreading
    FUNCTION FFT_2D(signal, method, threads) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:, :), INTENT(INOUT) :: signal
        CHARACTER(*), INTENT(IN) :: method
        INTEGER, OPTIONAL, INTENT(IN) :: threads
        COMPLEX(dp), DIMENSION(size(signal, 1), size(signal, 2)) :: RESULT

        IF (method == "NAFPack_FFT_2D") THEN
            RESULT = NAFPack_FFT_2D(signal)
        ELSE IF (method == "FFTW_FFT_2D" .AND. .NOT. present(threads)) THEN
            RESULT = FFTW_FFT_2D(signal)
        ELSE IF (method == "FFTW_FFT_2D" .AND. present(threads)) THEN
            RESULT = FFTW_FFT_2D_threads(signal, threads)
        ELSE
            STOP "ERROR : Wrong method for FFT_2D"
        END IF

    END FUNCTION FFT_2D

    !> Perform a 2D inverse Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a  inverse fast Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "NAFPack_IFFT_2D": Fast Fourier Transform using NAFPack
    !> - "FFTW_IFFT_2D": Fast Fourier Transform using FFTW
    !> - "FFTW_IFFT_2D" + threads: Fast Fourier Transform using FFTW with multithreading
    FUNCTION IFFT_2D(signal, method, threads) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:, :), INTENT(INOUT) :: signal
        CHARACTER(*), INTENT(IN) :: method
        INTEGER, OPTIONAL, INTENT(IN) :: threads
        COMPLEX(dp), DIMENSION(size(signal, 1), size(signal, 2)) :: RESULT

        IF (method == "NAFPack_IFFT_2D") THEN
            RESULT = NAFPack_IFFT_2D(signal)
        ELSE IF (method == "FFTW_IFFT_2D" .AND. .NOT. present(threads)) THEN
            RESULT = FFTW_IFFT_2D(signal)
        ELSE IF (method == "FFTW_IFFT_2D" .AND. present(threads)) THEN
            RESULT = FFTW_IFFT_2D_threads(signal, threads)
        ELSE
            STOP "ERROR : Wrong method for IFFT_1D"
        END IF

    END FUNCTION IFFT_2D

    !> Perform a 3D Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "FFTW_FFT_3D": Fast Fourier Transform using FFTW
    !> - "FFTW_FFT_3D" + threads: Fast Fourier Transform using FFTW with multithreading
    FUNCTION FFT_3D(signal, method, threads) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:, :, :), INTENT(INOUT) :: signal
        CHARACTER(*), INTENT(IN) :: method
        INTEGER, OPTIONAL, INTENT(IN) :: threads
        COMPLEX(dp), DIMENSION(size(signal, 1), size(signal, 2), size(signal, 3)) :: RESULT

        IF (method == "FFTW_FFT_3D" .AND. .NOT. present(threads)) THEN
            RESULT = FFTW_FFT_3D(signal)
        ELSE IF (method == "FFTW_FFT_3D" .AND. present(threads)) THEN
            RESULT = FFTW_FFT_3D_threads(signal, threads)
        ELSE
            STOP "ERROR : Wrong method for FFT_2D"
        END IF

    END FUNCTION FFT_3D

    !> Perform a 3D inverse Fast Fourier Transform on a signal
    !>
    !> This function takes a signal and performs a  inverse fast Fourier Transform using the specified method.
    !> The available methods are:
    !>
    !> - "FFTW_IFFT_3D": Fast Fourier Transform using FFTW
    !> - "FFTW_IFFT_3D" + threads: Fast Fourier Transform using FFTW with multithreading
    FUNCTION IFFT_3D(signal, method, threads) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:, :, :), INTENT(INOUT) :: signal
        CHARACTER(*), INTENT(IN) :: method
        INTEGER, OPTIONAL, INTENT(IN) :: threads
        COMPLEX(dp), DIMENSION(size(signal, 1), size(signal, 2), size(signal, 3)) :: RESULT

        IF (method == "FFTW_IFFT_3D" .AND. .NOT. present(threads)) THEN
            RESULT = FFTW_IFFT_3D(signal)
        ELSE IF (method == "IFFTW_IFFT_3D" .AND. present(threads)) THEN
            RESULT = FFTW_IFFT_3D_threads(signal, threads)
        ELSE
            STOP "ERROR : Wrong method for IFFT_1D"
        END IF

    END FUNCTION IFFT_3D

!################### FFTW ##########################################

    !> Perform a 1D Fast Fourier Transform on a signal using FFTW with multithreading
    FUNCTION FFTW_FFT_1D_threads(signal, threads) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:), INTENT(INOUT) :: signal
        INTEGER, INTENT(IN) :: threads
        COMPLEX(c_double_complex), DIMENSION(size(signal)) :: RESULT
        INTEGER :: error_init_thread
        TYPE(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        IF (error_init_thread == 0) STOP "ERROR : Thread FFTW initialization problem"

        CALL fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_1d(size(signal), signal, RESULT, FFTW_FORWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        CALL fftw_destroy_plan(plan)

        CALL fftw_cleanup_threads()

    END FUNCTION FFTW_FFT_1D_threads

    !> Perform a 1D Fast Fourier Transform on a signal using FFTW
    FUNCTION FFTW_FFT_1D(signal) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:), INTENT(INOUT) :: signal
        COMPLEX(c_double_complex), DIMENSION(size(signal)) :: RESULT
        TYPE(c_ptr) :: plan

        plan = fftw_plan_dft_1d(size(signal), signal, RESULT, FFTW_FORWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        CALL fftw_destroy_plan(plan)

    END FUNCTION FFTW_FFT_1D

    !> Perform a 1D inverse Fast Fourier Transform on a signal using FFTW
    FUNCTION FFTW_IFFT_1D(signal) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:), INTENT(INOUT) :: signal
        COMPLEX(c_double_complex), DIMENSION(size(signal)) :: RESULT
        TYPE(c_ptr) :: plan

        plan = fftw_plan_dft_1d(size(signal), signal, RESULT, FFTW_BACKWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        RESULT = RESULT / size(signal)
        CALL fftw_destroy_plan(plan)

    END FUNCTION FFTW_IFFT_1D

    !> Perform a 1D inverse Fast Fourier Transform on a signal using FFTW with multithreading
    FUNCTION FFTW_IFFT_1D_threads(signal, threads) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:), INTENT(INOUT) :: signal
        INTEGER, INTENT(IN) :: threads
        COMPLEX(c_double_complex), DIMENSION(size(signal)) :: RESULT
        INTEGER :: error_init_thread
        TYPE(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        IF (error_init_thread == 0) STOP "ERROR : Thread FFTW initialization problem"

        CALL fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_1d(size(signal), signal, RESULT, FFTW_BACKWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        RESULT = RESULT / size(signal)
        CALL fftw_destroy_plan(plan)

        CALL fftw_cleanup_threads()

    END FUNCTION FFTW_IFFT_1D_threads

    !> Perform a 2D Fast Fourier Transform on a signal using FFTW
    FUNCTION FFTW_FFT_2D(signal) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:, :), INTENT(INOUT) :: signal
        COMPLEX(c_double_complex), DIMENSION(size(signal, 1), size(signal, 2)) :: RESULT
        TYPE(c_ptr) :: plan

        plan = fftw_plan_dft_2d(size(signal, 2), size(signal, 1), signal, RESULT, FFTW_FORWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        CALL fftw_destroy_plan(plan)

    END FUNCTION FFTW_FFT_2D

    !> Perform a 2D Fast Fourier Transform on a signal using FFTW with multithreading
    FUNCTION FFTW_FFT_2D_threads(signal, threads) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:, :), INTENT(INOUT) :: signal
        INTEGER, INTENT(IN) :: threads
        COMPLEX(c_double_complex), DIMENSION(size(signal, 1), size(signal, 2)) :: RESULT
        INTEGER :: error_init_thread
        TYPE(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        IF (error_init_thread == 0) STOP "ERROR : Thread FFTW initialization problem"

        CALL fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_2d(size(signal, 2), size(signal, 1), signal, RESULT, FFTW_FORWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        CALL fftw_destroy_plan(plan)

        CALL fftw_cleanup_threads()

    END FUNCTION FFTW_FFT_2D_threads

    !> Perform a 2D inverse Fast Fourier Transform on a signal using FFTW
    FUNCTION FFTW_IFFT_2D(signal) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:, :), INTENT(INOUT) :: signal
        COMPLEX(c_double_complex), DIMENSION(size(signal, 1), size(signal, 2)) :: RESULT
        TYPE(c_ptr) :: plan

        plan = fftw_plan_dft_2d(size(signal, 2), size(signal, 1), signal, RESULT, FFTW_BACKWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        RESULT = RESULT / (size(signal, 1) * size(signal, 2))
        CALL fftw_destroy_plan(plan)

    END FUNCTION FFTW_IFFT_2D

    !> Perform a 2D inverse Fast Fourier Transform on a signal using FFTW with multithreading
    FUNCTION FFTW_IFFT_2D_threads(signal, threads) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:, :), INTENT(INOUT) :: signal
        INTEGER, INTENT(IN) :: threads
        COMPLEX(c_double_complex), DIMENSION(size(signal, 1), size(signal, 2)) :: RESULT
        INTEGER :: error_init_thread
        TYPE(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        IF (error_init_thread == 0) STOP "ERROR : Thread FFTW initialization problem"

        CALL fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_2d(size(signal, 2), size(signal, 1), signal, RESULT, FFTW_BACKWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        RESULT = RESULT / (size(signal, 1) * size(signal, 2))
        CALL fftw_destroy_plan(plan)

        CALL fftw_cleanup_threads()

    END FUNCTION FFTW_IFFT_2D_threads

    !> Perform a 3D Fast Fourier Transform on a signal using FFTW
    FUNCTION FFTW_FFT_3D(signal) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:, :, :), INTENT(INOUT) :: signal
        COMPLEX(c_double_complex), DIMENSION(size(signal, 1), size(signal, 2), size(signal, 3)) :: RESULT
        TYPE(c_ptr) :: plan

        plan = fftw_plan_dft_3d(size(signal, 3), size(signal, 2), size(signal, 1), signal, RESULT, FFTW_FORWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        CALL fftw_destroy_plan(plan)

    END FUNCTION FFTW_FFT_3D

    !> Perform a 3D Fast Fourier Transform on a signal using FFTW with multithreading
    FUNCTION FFTW_FFT_3D_threads(signal, threads) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:, :, :), INTENT(INOUT) :: signal
        INTEGER, INTENT(IN) :: threads
        COMPLEX(c_double_complex), DIMENSION(size(signal, 1), size(signal, 2), size(signal, 3)) :: RESULT
        INTEGER :: error_init_thread
        TYPE(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        IF (error_init_thread == 0) STOP "ERROR : Thread FFTW initialization problem"

        CALL fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_3d(size(signal, 3), size(signal, 2), size(signal, 1), signal, RESULT, FFTW_FORWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        CALL fftw_destroy_plan(plan)

        CALL fftw_cleanup_threads()

    END FUNCTION FFTW_FFT_3D_threads

    !> Perform a 3D inverse Fast Fourier Transform on a signal using FFTW
    FUNCTION FFTW_IFFT_3D(signal) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:, :, :), INTENT(INOUT) :: signal
        COMPLEX(c_double_complex), DIMENSION(size(signal, 1), size(signal, 2), size(signal, 3)) :: RESULT
        TYPE(c_ptr) :: plan

        plan = fftw_plan_dft_3d(size(signal, 3), size(signal, 2), size(signal, 1), signal, RESULT, FFTW_BACKWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        RESULT = RESULT / (size(signal, 1) * size(signal, 2) * size(signal, 3))
        CALL fftw_destroy_plan(plan)

    END FUNCTION FFTW_IFFT_3D

    !> Perform a 3D inverse Fast Fourier Transform on a signal using FFTW with multithreading
    FUNCTION FFTW_IFFT_3D_threads(signal, threads) RESULT(RESULT)

        COMPLEX(c_double_complex), DIMENSION(:, :, :), INTENT(INOUT) :: signal
        INTEGER, INTENT(IN) :: threads
        COMPLEX(c_double_complex), DIMENSION(size(signal, 1), size(signal, 2), size(signal, 3)) :: RESULT
        INTEGER :: error_init_thread
        TYPE(c_ptr) :: plan

        error_init_thread = fftw_init_threads()
        IF (error_init_thread == 0) STOP "ERROR : Thread FFTW initialization problem"

        CALL fftw_plan_with_nthreads(threads)

        plan = fftw_plan_dft_3d(size(signal, 3), size(signal, 2), size(signal, 1), signal, RESULT, FFTW_BACKWARD, FFTW_ESTIMATE)
        CALL fftw_execute_dft(plan, signal, RESULT)
        RESULT = RESULT / (size(signal, 1) * size(signal, 2) * size(signal, 3))
        CALL fftw_destroy_plan(plan)

        CALL fftw_cleanup_threads()

    END FUNCTION FFTW_IFFT_3D_threads

!################### NAFPack ##########################################

    !> Perform a 1D Discrete Fourier Transform on a signal
    FUNCTION NAFPack_DFT_1D(signal) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:), INTENT(IN) :: signal
        COMPLEX(dp), DIMENSION(size(signal)) :: RESULT
        COMPLEX(dp) :: S
        INTEGER :: N, i, k, j

        N = size(signal)

        IF (N == 1) THEN
            RESULT = signal
        ELSE
            DO i = 1, N

                k = i - 1
                S = (0.d0, 0.d0)

                DO j = 1, N
                    S = S + signal(j) * exp((-2 * pi * im * k * (j - 1)) / N)
                END DO

                RESULT(i) = S

            END DO
        END IF

    END FUNCTION NAFPack_DFT_1D

    !> Compute the complex exponential factors for the FFT
    FUNCTION fun_omega(N) RESULT(RESULT)

        INTEGER, INTENT(IN) :: N
        COMPLEX(dp), DIMENSION(N/2) :: RESULT
        INTEGER :: i

        DO i = 1, N / 2
            RESULT(i) = exp(-2 * Im * pi * (i - 1) / N)
        END DO

    END FUNCTION fun_omega

    !> Perform a 1D Fast Fourier Transform (Cooley-Tukey) on a signal
    RECURSIVE FUNCTION NAFPack_FFT_1D(signal) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:), INTENT(IN) :: signal
        COMPLEX(dp), DIMENSION(size(signal)) :: RESULT
        COMPLEX(dp), DIMENSION(size(signal)/2) :: f_pair, f_impair, omega
        INTEGER :: N

        N = size(signal)

        IF (mod(N, 2) == 0) THEN
            f_pair = NAFPack_FFT_1D(signal(1 :: 2))
            f_impair = NAFPack_FFT_1D(signal(2 :: 2))

            omega = fun_omega(N)

            RESULT(1:N / 2) = f_pair + f_impair * omega
            RESULT(N / 2 + 1:) = f_pair - f_impair * omega
        ELSE
            RESULT = NAFPack_DFT_1D(signal)
        END IF
    END FUNCTION NAFPack_FFT_1D

    !> Perform a 1D inverse Fast Fourier Transform on a signal
    FUNCTION NAFPack_IFFT_1D(f_signal) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:), INTENT(IN) :: f_signal
        COMPLEX(dp), DIMENSION(size(f_signal)) :: RESULT
        COMPLEX(dp), DIMENSION(size(f_signal)) :: f_conjugate
        INTEGER :: N

        N = size(f_signal)

        f_conjugate = conjg(f_signal)

        RESULT = NAFPack_FFT_1D(f_conjugate)
        RESULT = conjg(RESULT)
        RESULT = RESULT / N

    END FUNCTION NAFPack_IFFT_1D

    !> Perform a 2D Fast Fourier Transform on a signal
    FUNCTION NAFPack_FFT_2D(signal) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:, :), INTENT(IN) :: signal
        COMPLEX(dp), DIMENSION(size(signal, 1), size(signal, 2)) :: RESULT
        INTEGER :: Nx, Ny, i

        Nx = size(signal, 1)
        Ny = size(signal, 2)

        DO i = 1, Nx
            RESULT(i, :) = NAFPack_FFT_1D(signal(i, :))
        END DO

        DO i = 1, Ny
            RESULT(:, i) = NAFPack_FFT_1D(RESULT(:, i))
        END DO

    END FUNCTION NAFPack_FFT_2D

    !> Perform a 2D inverse Fast Fourier Transform on a signal
    FUNCTION NAFPack_IFFT_2D(f_signal) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:, :), INTENT(IN) :: f_signal
        COMPLEX(dp), DIMENSION(size(f_signal, 1), size(f_signal, 2)) :: RESULT
        INTEGER :: Nx, Ny, i

        Nx = size(f_signal, 1)
        Ny = size(f_signal, 2)

        DO i = 1, Nx
            RESULT(i, :) = NAFPack_IFFT_1D(f_signal(i, :))
        END DO

        DO i = 1, Ny
            RESULT(:, i) = NAFPack_IFFT_1D(RESULT(:, i))
        END DO

    END FUNCTION NAFPack_IFFT_2D

END MODULE NAFPack_fft
