MODULE NAFPack_fft

    USE NAFPack_constant
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: FFT, FFT2
    PUBLIC :: IFFT, IFFT2
    PUBLIC :: DFT

    CONTAINS
    
    !Discrete Fourier transform 1D
    FUNCTION DFT(signal) RESULT(result)

        COMPLEX(dp), DIMENSION(:), INTENT(IN) :: signal
        COMPLEX(dp), DIMENSION(SIZE(signal)) :: result
        COMPLEX(dp) :: S
        INTEGER :: N, i, k, j

        N = SIZE(signal)

        IF (N == 1) THEN
            result = signal
        ELSE
            DO i = 1, N

                k=i-1
                S=(0.d0, 0.d0)

                DO j = 1, N
                    S = S + signal(j) * EXP((-2 * pi * im * k * (j - 1)) / N)
                END DO

                result(i) = S

            END DO
        END IF
    
    END FUNCTION DFT

    FUNCTION fun_omega(N) RESULT(result)

        INTEGER, INTENT(IN) :: N
        COMPLEX(dp), DIMENSION(N / 2) :: result
        INTEGER :: i

        DO i = 1, N / 2
            result(i) = EXP(-2 * Im * pi * (i - 1) / N)
        END DO

    END FUNCTION fun_omega

    !Fast Fourier transform 1D Cooley-Tukey
    RECURSIVE FUNCTION FFT(signal) RESULT(result)

        COMPLEX(dp), DIMENSION(:), INTENT(IN) :: signal
        COMPLEX(dp), DIMENSION(SIZE(signal)) :: result
        COMPLEX(dp), DIMENSION(SIZE(signal)/2) :: f_pair, f_impair,omega
        INTEGER :: N

        N = SIZE(signal)

        IF(MOD(N, 2) == 0)THEN
            f_pair = FFT(signal(1::2))
            f_impair = FFT(signal(2::2))

            omega = fun_omega(N)

            result(1:N/2) = f_pair + f_impair * omega
            result(N/2+1:) = f_pair - f_impair * omega
        ELSE
            result = DFT(signal)
        END IF
    END FUNCTION FFT

    !The inverse of fft
    FUNCTION IFFT(f_signal) RESULT(result)

        COMPLEX(dp), DIMENSION(:), INTENT(IN) :: f_signal
        COMPLEX(dp), DIMENSION(SIZE(f_signal)) :: result
        COMPLEX(dp), DIMENSION(SIZE(f_signal)) :: f_conjugate
        INTEGER :: N

        N = SIZE(f_signal)

        f_conjugate = CONJG(f_signal)

        result = FFT(f_conjugate)
        result = CONJG(result)
        result = result / N

    END FUNCTION IFFT

    !Fast Fourier transform 2D Cooley-Tukey
    FUNCTION FFT2(signal) RESULT(result)

        COMPLEX(dp), DIMENSION(:, :), INTENT(IN) :: signal
        COMPLEX(dp), DIMENSION(SIZE(signal, 1), SIZE(signal, 2)) :: result
        INTEGER :: Nx, Ny, i

        Nx = SIZE(signal,1)
        Ny = SIZE(signal,2)

        DO i = 1, Nx
            result(i, :) = FFT(signal(i, :))
        END DO

        DO i = 1, Ny
            result(:, i) = FFT(result(:, i))
        END DO

    END FUNCTION FFT2


    !The inverse of fft2
    FUNCTION IFFT2(f_signal) RESULT(result)

        COMPLEX(dp), DIMENSION(:, :), INTENT(IN) :: f_signal
        COMPLEX(dp), DIMENSION(SIZE(f_signal, 1), SIZE(f_signal, 2)) :: result
        INTEGER :: Nx, Ny, i

        Nx = SIZE(f_signal, 1)
        Ny = SIZE(f_signal, 2)

        DO i = 1, Nx
            result(i, :) = IFFT(f_signal(i, :))
        END DO

        DO i = 1, Ny
            result(:, i) = IFFT(result(:, i))
        END DO

    END FUNCTION IFFT2

END MODULE NAFPack_fft