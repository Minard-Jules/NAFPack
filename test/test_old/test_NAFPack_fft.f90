module test_NAFPack_fft

    use NAFPack_fft
    use NAFPack_meshgrid
    use NAFPack_constant
    use NAFPack_kinds
    use NAFPack_ANSI, only: ColorsAscii


    implicit none(type, external)

    private
    public :: test_FFT

contains

    subroutine test_FFT(stat)

        logical, intent(inout) :: stat

        integer, parameter :: Mx = 10, My = 20
        integer, parameter :: Mx2 = 5, My2 = 6
        real(dp), dimension(Mx) :: xlist
        real(dp), dimension(My) :: ylist
        complex(dp), dimension(Mx) :: S, fS
        complex(dp), dimension(Mx) :: fs_exact_DFT, fs_DFT, diff_DFT
        complex(dp), dimension(Mx) :: fs_exact_FFT, fs_FFT, diff_FFT
        complex(dp), dimension(Mx) :: s_exact_IFFT, s_IFFT, diff_IFFT
        complex(dp), dimension(My, Mx) :: fs_exact_FFT2, fs_FFT2, diff_FFT2
        complex(dp), dimension(My2, Mx2) :: s_exact_IFFT2, s_IFFT2, diff_IFFT2
        real(dp), dimension(:, :), allocatable :: X, Y
        complex(dp), dimension(My, Mx) :: S2D
        complex(dp), dimension(My2, Mx2) :: fS2S
        integer :: i
        type(ColorsAscii) :: colors

        CALL colors%init()

        do i = 1, Mx
            xlist(i) = i - 1
        end do
        xlist = xlist * 2 * pi / Mx

        do i = 1, My
            ylist(i) = i - 1
        end do
        ylist = ylist * 2 * pi / My

        !==================================================================
        !test DFT
        S = exp(-im * xlist)

        fs_exact_DFT = (0.d0, 0.d0)
        fs_exact_DFT(Mx) = (10.d0, 0.d0)

        fs_DFT = FFT_1D(S, "NAFPack_DFT")

        diff_DFT = fs_exact_DFT - fs_DFT
        if (maxval(abs(diff_DFT)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"DFT", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"DFT", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test FFT
        S = 1 + (cos(xlist))**2 + 3 * sin(4 * xlist)

        fs_exact_FFT = (0.d0, 0.d0)
        fs_exact_FFT(1) = (15.d0, 0.d0)
        fs_exact_FFT(3) = (2.5d0, 0.d0)
        fs_exact_FFT(5) = (0.d0, -15.d0)
        fs_exact_FFT(7) = (0.d0, 15.d0)
        fs_exact_FFT(9) = (2.5d0, 0.d0)

        fs_FFT = FFT_1D(S, "NAFPack_FFT_1D")

        diff_FFT = fs_exact_FFT - fs_FFT
        if (maxval(abs(diff_FFT)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"FFT", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"FFT", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test FFTW
        S = 1 + (cos(xlist))**2 + 3 * sin(4 * xlist)

        fs_exact_FFT = (0.d0, 0.d0)
        fs_exact_FFT(1) = (15.d0, 0.d0)
        fs_exact_FFT(3) = (2.5d0, 0.d0)
        fs_exact_FFT(5) = (0.d0, -15.d0)
        fs_exact_FFT(7) = (0.d0, 15.d0)
        fs_exact_FFT(9) = (2.5d0, 0.d0)

        fs_FFT = FFT_1D(S, "FFTW_FFT_1D")

        diff_FFT = fs_exact_FFT - fs_FFT
        if (maxval(abs(diff_FFT)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"FFTW", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"FFTW", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test FFTW threads
        S = 1 + (cos(xlist))**2 + 3 * sin(4 * xlist)

        fs_exact_FFT = (0.d0, 0.d0)
        fs_exact_FFT(1) = (15.d0, 0.d0)
        fs_exact_FFT(3) = (2.5d0, 0.d0)
        fs_exact_FFT(5) = (0.d0, -15.d0)
        fs_exact_FFT(7) = (0.d0, 15.d0)
        fs_exact_FFT(9) = (2.5d0, 0.d0)

        fs_FFT = FFT_1D(S, "FFTW_FFT_1D", threads=4)

        diff_FFT = fs_exact_FFT - fs_FFT
        if (maxval(abs(diff_FFT)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"FFTW threads", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"FFTW threads", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test IFFT
        fS = (0.d0, 0.d0)
        fS(1) = Mx
        fS(2) = 2 * Mx
        fS(Mx) = 2 * Mx
        fS(3) = -im * Mx / 2
        fS(Mx - 1) = im * Mx / 2

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

        s_IFFT = IFFT_1D(fS, "NAFPack_IFFT_1D")

        diff_IFFT = s_exact_IFFT - s_IFFT
        if (maxval(abs(diff_FFT)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"IFFT", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"IFFT", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test IFFTW
        fS = (0.d0, 0.d0)
        fS(1) = Mx
        fS(2) = 2 * Mx
        fS(Mx) = 2 * Mx
        fS(3) = -im * Mx / 2
        fS(Mx - 1) = im * Mx / 2

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

        s_IFFT = IFFT_1D(fS, "FFTW_IFFT_1D")

        diff_IFFT = s_exact_IFFT - s_IFFT
        if (maxval(abs(diff_FFT)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"IFFTW", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"IFFTW", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test IFFTW threads
        fS = (0.d0, 0.d0)
        fS(1) = Mx
        fS(2) = 2 * Mx
        fS(Mx) = 2 * Mx
        fS(3) = -im * Mx / 2
        fS(Mx - 1) = im * Mx / 2

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

        s_IFFT = IFFT_1D(fS, "FFTW_IFFT_1D", threads=4)

        diff_IFFT = s_exact_IFFT - s_IFFT
        if (maxval(abs(diff_FFT)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"IFFTW threads", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"IFFTW threads", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test FFT2
        call meshgrid(xlist, ylist, X, Y)

        S2D = sin(2 * X) + sin(Y)

        fs_exact_FFT2 = (0.d0, 0.d0)
        fs_exact_FFT2(1, 3) = (0.d0, -100.d0)
        fs_exact_FFT2(1, Mx - 1) = (0.d0, 100.d0)
        fs_exact_FFT2(2, 1) = (0d0, -100.d0)
        fs_exact_FFT2(My, 1) = (0.d0, 100.d0)

        fs_FFT2 = FFT_2D(S2D, "NAFPack_FFT_2D")

        diff_FFT2 = fs_exact_FFT2 - fs_FFT2
        if (maxval(abs(diff_FFT2)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"FFT2", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"FFT2", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test FFTW2
        call meshgrid(xlist, ylist, X, Y)
        S2D = sin(2 * X) + sin(Y)

        fs_exact_FFT2 = (0.d0, 0.d0)
        fs_exact_FFT2(1, 3) = (0.d0, -100.d0)
        fs_exact_FFT2(1, Mx - 1) = (0.d0, 100.d0)
        fs_exact_FFT2(2, 1) = (0d0, -100.d0)
        fs_exact_FFT2(My, 1) = (0.d0, 100.d0)

        fs_FFT2 = FFT_2D(S2D, "FFTW_FFT_2D")

        diff_FFT2 = fs_exact_FFT2 - fs_FFT2
        if (maxval(abs(diff_FFT2)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"FFTW2", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"FFTW2", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test FFTW2 threads
        call meshgrid(xlist, ylist, X, Y)
        S2D = sin(2 * X) + sin(Y)

        fs_exact_FFT2 = (0.d0, 0.d0)
        fs_exact_FFT2(1, 3) = (0.d0, -100.d0)
        fs_exact_FFT2(1, Mx - 1) = (0.d0, 100.d0)
        fs_exact_FFT2(2, 1) = (0d0, -100.d0)
        fs_exact_FFT2(My, 1) = (0.d0, 100.d0)

        fs_FFT2 = FFT_2D(S2D, "FFTW_FFT_2D", threads=4)

        diff_FFT2 = fs_exact_FFT2 - fs_FFT2
        if (maxval(abs(diff_FFT2)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"FFTW2 threads", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"FFTW2 threads", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test IFFT2
        fS2S = (0.d0, 0.d0)
        fS2S(1, 1) = Mx2

        s_exact_IFFT2 = (0.16666667d0, 0.d0)

        s_IFFT2 = IFFT_2D(fS2S, "NAFPack_IFFT_2D")

        diff_IFFT2 = s_exact_IFFT2 - s_IFFT2
        if (maxval(abs(diff_FFT2)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"IFFT2", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"IFFT2", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test IFFTW2
        fS2S = (0.d0, 0.d0)
        fS2S(1, 1) = Mx2

        s_exact_IFFT2 = (0.16666667d0, 0.d0)

        s_IFFT2 = IFFT_2D(fS2S, "FFTW_IFFT_2D")

        diff_IFFT2 = s_exact_IFFT2 - s_IFFT2
        if (maxval(abs(diff_FFT2)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"IFFTW2", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"IFFTW2", " :: ECHEC"//colors%reset
            stat = .true.
        end if

        !==================================================================
        !test IFFTW2
        fS2S = (0.d0, 0.d0)
        fS2S(1, 1) = Mx2

        s_exact_IFFT2 = (0.16666667d0, 0.d0)

        s_IFFT2 = IFFT_2D(fS2S, "FFTW_IFFT_2D", threads=4)

        diff_IFFT2 = s_exact_IFFT2 - s_IFFT2
        if (maxval(abs(diff_FFT2)) < TOL_CONVERGENCE) then
            write (*, '(A,T50,A,A)') colors%green//"IFFTW2 threads", " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//"IFFTW2 threads", " :: ECHEC"//colors%reset
            stat = .true.
        end if

    end subroutine test_FFT

end module test_NAFPack_fft
