MODULE test_NAFPack_fft

    USE NAFPack_fft
    USE NAFPack_meshgrid
    USE NAFPack_constant

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: test_FFT

    CONTAINS 

    SUBROUTINE test_FFT(stat)

        LOGICAL, INTENT(INOUT) :: stat

        INTEGER, PARAMETER :: Mx = 10, My = 20
        INTEGER, PARAMETER :: Mx2 = 5, My2 = 6
        REAL(dp), DIMENSION(Mx) :: xlist
        REAL(dp), DIMENSION(My) :: ylist
        COMPLEX(dp), DIMENSION(Mx) :: S,fS
        COMPLEX(dp), DIMENSION(Mx) :: fs_exact_DFT, fs_DFT, diff_DFT
        COMPLEX(dp), DIMENSION(Mx) :: fs_exact_FFT, fs_FFT, diff_FFT
        COMPLEX(dp), DIMENSION(Mx) :: s_exact_IFFT, s_IFFT, diff_IFFT
        COMPLEX(dp), DIMENSION(My, Mx) :: fs_exact_FFT2, fs_FFT2, diff_FFT2
        COMPLEX(dp), DIMENSION(My2, Mx2) :: s_exact_IFFT2, s_IFFT2, diff_IFFT2
        REAL(dp), DIMENSION(My, Mx) :: X, Y
        COMPLEX(dp), DIMENSION(My, Mx) :: S2D
        COMPLEX(dp), DIMENSION(My2, Mx2) :: fS2S
        INTEGER :: i


        DO i = 1, Mx
            xlist(i) = i - 1
        END DO
        xlist = xlist * 2 * pi / Mx

        DO i = 1, My
            ylist(i) = i - 1
        END DO
        ylist = ylist * 2 * pi / My


        !==================================================================
        !test DFT
        S = EXP(- im * xlist)

        fs_exact_DFT = (0.d0, 0.d0)
        fs_exact_DFT(Mx) = (10.d0, 0.d0)

        fs_DFT = FFT_1D(S, "NAFPack_DFT")

        diff_DFT = fs_exact_DFT - fs_DFT
        IF(MAXVAL(ABS(diff_DFT)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"DFT", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"DFT", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test FFT
        S = 1 + (COS(xlist))**2 + 3 * SIN(4 * xlist)

        fs_exact_FFT = (0.d0, 0.d0)
        fs_exact_FFT(1) = (15.d0, 0.d0)
        fs_exact_FFT(3) = (2.5d0, 0.d0)
        fs_exact_FFT(5) = (0.d0, -15.d0)
        fs_exact_FFT(7) = (0.d0, 15.d0)
        fs_exact_FFT(9) = (2.5d0, 0.d0)

        fs_FFT = FFT_1D(S, "NAFPack_FFT_1D")

        diff_FFT = fs_exact_FFT - fs_FFT
        IF(MAXVAL(ABS(diff_FFT)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"FFT", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"FFT", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test FFTW
        S = 1 + (COS(xlist))**2 + 3 * SIN(4 * xlist)

        fs_exact_FFT = (0.d0, 0.d0)
        fs_exact_FFT(1) = (15.d0, 0.d0)
        fs_exact_FFT(3) = (2.5d0, 0.d0)
        fs_exact_FFT(5) = (0.d0, -15.d0)
        fs_exact_FFT(7) = (0.d0, 15.d0)
        fs_exact_FFT(9) = (2.5d0, 0.d0)

        fs_FFT = FFT_1D(S, "FFTW_FFT_1D")

        diff_FFT = fs_exact_FFT - fs_FFT
        IF(MAXVAL(ABS(diff_FFT)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"FFTW", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"FFTW", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test FFTW threads
        S = 1 + (COS(xlist))**2 + 3 * SIN(4 * xlist)

        fs_exact_FFT = (0.d0, 0.d0)
        fs_exact_FFT(1) = (15.d0, 0.d0)
        fs_exact_FFT(3) = (2.5d0, 0.d0)
        fs_exact_FFT(5) = (0.d0, -15.d0)
        fs_exact_FFT(7) = (0.d0, 15.d0)
        fs_exact_FFT(9) = (2.5d0, 0.d0)

        fs_FFT = FFT_1D(S, "FFTW_FFT_1D", threads=4)

        diff_FFT = fs_exact_FFT - fs_FFT
        IF(MAXVAL(ABS(diff_FFT)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"FFTW threads", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"FFTW threads", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test IFFT
        fS = (0.d0, 0.d0)
        fS(1) = Mx
        fS(2) = 2 * Mx
        fS(Mx) = 2 * Mx
        fS(3) = - im * Mx / 2
        fS(Mx-1) = im * Mx / 2

        s_exact_IFFT = [(5.d0, 0.d0), &
                    (5.19d0, 0.d0), &
                    (2.82d0, 0.d0), &
                    (-0.82d0, 0.d0), &
                    (-3.19d0, 0.d0), &
                    (-3.d0, 0.d0), &
                    (-1.29d0, 0.d0), &
                    (0.35d0, 0.d0), &
                    (1.65d0, 0.d0), &
                    (3.29d0, 0.d0)]

        s_IFFT = IFFT_1D(fS,"NAFPack_IFFT_1D")

        diff_IFFT = s_exact_IFFT - s_IFFT
        IF(MAXVAL(ABS(diff_FFT)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"IFFT", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"IFFT", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test IFFTW
        fS = (0.d0, 0.d0)
        fS(1) = Mx
        fS(2) = 2 * Mx
        fS(Mx) = 2 * Mx
        fS(3) = - im * Mx / 2
        fS(Mx-1) = im * Mx / 2

        s_exact_IFFT = [(5.d0, 0.d0), &
                    (5.19d0, 0.d0), &
                    (2.82d0, 0.d0), &
                    (-0.82d0, 0.d0), &
                    (-3.19d0, 0.d0), &
                    (-3.d0, 0.d0), &
                    (-1.29d0, 0.d0), &
                    (0.35d0, 0.d0), &
                    (1.65d0, 0.d0), &
                    (3.29d0, 0.d0)]

        s_IFFT = IFFT_1D(fS,"FFTW_IFFT_1D")

        diff_IFFT = s_exact_IFFT - s_IFFT
        IF(MAXVAL(ABS(diff_FFT)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"IFFTW", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"IFFTW", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test IFFTW threads
        fS = (0.d0, 0.d0)
        fS(1) = Mx
        fS(2) = 2 * Mx
        fS(Mx) = 2 * Mx
        fS(3) = - im * Mx / 2
        fS(Mx-1) = im * Mx / 2

        s_exact_IFFT = [(5.d0, 0.d0), &
                    (5.19d0, 0.d0), &
                    (2.82d0, 0.d0), &
                    (-0.82d0, 0.d0), &
                    (-3.19d0, 0.d0), &
                    (-3.d0, 0.d0), &
                    (-1.29d0, 0.d0), &
                    (0.35d0, 0.d0), &
                    (1.65d0, 0.d0), &
                    (3.29d0, 0.d0)]

        s_IFFT = IFFT_1D(fS,"FFTW_IFFT_1D", threads=4)

        diff_IFFT = s_exact_IFFT - s_IFFT
        IF(MAXVAL(ABS(diff_FFT)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"IFFTW threads", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"IFFTW threads", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test FFT2
        CALL meshgrid(xlist, ylist, X, Y)
        S2D = SIN(2 * X) + SIN(Y)

        fs_exact_FFT2 = (0.d0, 0.d0)
        fs_exact_FFT2(1, 3) = (0.d0, -100.d0) 
        fs_exact_FFT2(1, Mx - 1) = (0.d0, 100.d0)
        fs_exact_FFT2(2, 1) = (0d0, -100.d0)
        fs_exact_FFT2(My, 1) = (0.d0, 100.d0)

        fs_FFT2 = FFT_2D(S2D, "NAFPack_FFT_2D")

        diff_FFT2 = fs_exact_FFT2 - fs_FFT2
        IF(MAXVAL(ABS(diff_FFT2)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"FFT2", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"FFT2", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test FFTW2
        CALL meshgrid(xlist, ylist, X, Y)
        S2D = SIN(2 * X) + SIN(Y)

        fs_exact_FFT2 = (0.d0, 0.d0)
        fs_exact_FFT2(1, 3) = (0.d0, -100.d0) 
        fs_exact_FFT2(1, Mx - 1) = (0.d0, 100.d0)
        fs_exact_FFT2(2, 1) = (0d0, -100.d0)
        fs_exact_FFT2(My, 1) = (0.d0, 100.d0)

        fs_FFT2 = FFT_2D(S2D, "FFTW_FFT_2D")

        diff_FFT2 = fs_exact_FFT2 - fs_FFT2
        IF(MAXVAL(ABS(diff_FFT2)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"FFTW2", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"FFTW2", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test FFTW2 threads
        CALL meshgrid(xlist, ylist, X, Y)
        S2D = SIN(2 * X) + SIN(Y)

        fs_exact_FFT2 = (0.d0, 0.d0)
        fs_exact_FFT2(1, 3) = (0.d0, -100.d0) 
        fs_exact_FFT2(1, Mx - 1) = (0.d0, 100.d0)
        fs_exact_FFT2(2, 1) = (0d0, -100.d0)
        fs_exact_FFT2(My, 1) = (0.d0, 100.d0)

        fs_FFT2 = FFT_2D(S2D, "FFTW_FFT_2D", threads=4)

        diff_FFT2 = fs_exact_FFT2 - fs_FFT2
        IF(MAXVAL(ABS(diff_FFT2)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"FFTW2 threads", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"FFTW2 threads", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test IFFT2
        fS2S = (0.d0, 0.d0)
        fS2S(1, 1) = Mx2

        s_exact_IFFT2 = (0.16666667d0, 0.d0)

        s_IFFT2 = IFFT_2D(fS2S,"NAFPack_IFFT_2D")

        diff_IFFT2 = s_exact_IFFT2 - s_IFFT2
        IF(MAXVAL(ABS(diff_FFT2)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"IFFT2", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"IFFT2", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test IFFTW2
        fS2S = (0.d0, 0.d0)
        fS2S(1, 1) = Mx2

        s_exact_IFFT2 = (0.16666667d0, 0.d0)

        s_IFFT2 = IFFT_2D(fS2S, "FFTW_IFFT_2D")

        diff_IFFT2 = s_exact_IFFT2 - s_IFFT2
        IF(MAXVAL(ABS(diff_FFT2)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"IFFTW2", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"IFFTW2", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

        !==================================================================
        !test IFFTW2
        fS2S = (0.d0, 0.d0)
        fS2S(1, 1) = Mx2

        s_exact_IFFT2 = (0.16666667d0, 0.d0)

        s_IFFT2 = IFFT_2D(fS2S, "FFTW_IFFT_2D", threads=4)

        diff_IFFT2 = s_exact_IFFT2 - s_IFFT2
        IF(MAXVAL(ABS(diff_FFT2)) < epsi) THEN
            WRITE(*,'(A,T40,A,A)') green_color//"IFFTW2 threads", " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//"IFFTW2 threads", " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

    END SUBROUTINE test_FFT
    
END MODULE test_NAFPack_fft