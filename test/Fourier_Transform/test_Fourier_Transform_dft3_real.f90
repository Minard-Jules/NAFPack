module test_Fourier_Transform_dft3_real

    use NAFPack_kinds, only: dp, sp, qp
    use NAFPack_constant, only: TOL_TEST_sp, TOL_TEST_dp, TOL_TEST_qp
    use NAFPack_constant, only: pi_sp, pi_dp, pi_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_loop_method, only: LoopMethod, init_loop_method
    use NAFPack_Fourier_Transform, only: Fourier_Transform

    implicit none(type, external)

    private
    public :: collect_Fourier_Transform_dft3

contains

    subroutine collect_Fourier_Transform_dft3(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("dft real(sp)  3D", test_real_dft3_sp), &
                    new_unittest("idft real(sp) 3D", test_real_idft3_sp), &
                    new_unittest("dft real(dp)  3D", test_real_dft3_dp), &
                    new_unittest("idft real(dp) 3D", test_real_idft3_dp), &
                    new_unittest("dft real(qp)  3D", test_real_dft3_qp), &
                    new_unittest("idft real(qp) 3D", test_real_idft3_qp) &
                    ]
    end subroutine collect_Fourier_Transform_dft3

    ! ----------- SIMPLE PRECISION -----------
    subroutine test_real_dft3_sp(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: Nx = 4, Ny = 5, Nz = 6, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(sp), dimension(Nx) :: xlist
        real(sp), dimension(Ny) :: ylist
        real(sp), dimension(Nz) :: zlist
        real(sp), dimension(number_tests, Nx, Ny, Nz) :: signal
        real(sp), dimension(Nx, Ny, Nz) :: signal_used
        complex(sp), dimension(number_tests, Nx, Ny, Nz) :: fs_exact
        complex(sp), dimension(Nx, Ny, Nz) :: fs_DFT3
        integer :: i, j, k, freqx, freqy, freqz

        xlist = [(real((i - 1), kind=sp), i=1, Nx)] * 2 * pi_sp / Nx
        ylist = [(real((j - 1), kind=sp), j=1, Ny)] * 2 * pi_sp / Ny
        zlist = [(real((k - 1), kind=sp), k=1, Nz)] * 2 * pi_sp / Nz

        ! test 1: Dirac delta function
        signal(1, :, :, :) = 0._sp
        signal(1, 1, 1, 1) = 1._sp
        fs_exact(1, :, :, :) = cmplx(1._sp, 0._sp, kind=sp)

        ! test 2: onde cosinus 3D
        freqx = 2
        freqy = 1
        freqz = 3
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nz
                    signal(2, i, j, k) = cos(freqx * xlist(i) + freqy * ylist(j) + freqz * zlist(k))
                end do
            end do
        end do
        fs_exact(2, :, :, :) = cmplx(0._sp, 0._sp, kind=sp)
        fs_exact(2, freqx + 1, freqy + 1, freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=sp) / 2, 0._sp, kind=sp)
        fs_exact(2, Nx - freqx + 1, Ny - freqy + 1, Nz - freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=sp) / 2, 0._sp, kind=sp)

        ! test 3: constant signal (DC component)
        signal(3, :, :, :) = 1._sp
        fs_exact(3, :, :, :) = cmplx(0._sp, 0._sp, kind=sp)
        fs_exact(3, 1, 1, 1) = cmplx(real(Nx * Ny * Nz, kind=sp), 0._sp, kind=sp)

        do i = 1, number_tests
            signal_used = signal(i, :, :, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true., num_threads=4)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_sp)
            if (allocated(error)) return
        end do
    end subroutine test_real_dft3_sp

    subroutine test_real_idft3_sp(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: Nx = 4, Ny = 5, Nz = 6, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(sp), dimension(Nx) :: xlist
        real(sp), dimension(Ny) :: ylist
        real(sp), dimension(Nz) :: zlist
        real(sp), dimension(number_tests, Nx, Ny, Nz) :: fs
        real(sp), dimension(Nx, Ny, Nz) :: fs_used
        complex(sp), dimension(number_tests, Nx, Ny, Nz) :: signal_exact
        complex(sp), dimension(Nx, Ny, Nz) :: signal_IDFT3
        integer :: i, j, k, freqx, freqy, freqz

        xlist = [(real((i - 1), kind=sp), i=1, Nx)] * 2 * pi_sp / Nx
        ylist = [(real((j - 1), kind=sp), j=1, Ny)] * 2 * pi_sp / Ny
        zlist = [(real((k - 1), kind=sp), k=1, Nz)] * 2 * pi_sp / Nz

        ! test 1: Dirac delta function
        signal_exact(1, :, :, :) = 0._sp
        signal_exact(1, 1, 1, 1) = 1._sp
        fs(1, :, :, :) = cmplx(1._sp, 0._sp, kind=sp)

        ! test 2: onde cosinus 3D
        freqx = 2
        freqy = 1
        freqz = 3
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nz
                    signal_exact(2, i, j, k) = cos(freqx * xlist(i) + freqy * ylist(j) + freqz * zlist(k))
                end do
            end do
        end do
        fs(2, :, :, :) = cmplx(0._sp, 0._sp, kind=sp)
        fs(2, freqx + 1, freqy + 1, freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=sp) / 2, 0._sp, kind=sp)
        fs(2, Nx - freqx + 1, Ny - freqy + 1, Nz - freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=sp) / 2, 0._sp, kind=sp)

        ! test 3: constant signal (DC component)
        signal_exact(3, :, :, :) = 1._sp
        fs(3, :, :, :) = cmplx(0._sp, 0._sp, kind=sp)
        fs(3, 1, 1, 1) = cmplx(real(Nx * Ny * Nz, kind=sp), 0._sp, kind=sp)

        do i = 1, number_tests
            fs_used = fs(i, :, :, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_sp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true., num_threads=4)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_sp)
            if (allocated(error)) return
        end do
    end subroutine test_real_idft3_sp

    subroutine test_real_dft3_dp(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: Nx = 4, Ny = 5, Nz = 6, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(dp), dimension(Nx) :: xlist
        real(dp), dimension(Ny) :: ylist
        real(dp), dimension(Nz) :: zlist
        real(dp), dimension(number_tests, Nx, Ny, Nz) :: signal
        real(dp), dimension(Nx, Ny, Nz) :: signal_used
        complex(dp), dimension(number_tests, Nx, Ny, Nz) :: fs_exact
        complex(dp), dimension(Nx, Ny, Nz) :: fs_DFT3
        integer :: i, j, k, freqx, freqy, freqz

        xlist = [(real((i - 1), kind=dp), i=1, Nx)] * 2 * pi_dp / Nx
        ylist = [(real((j - 1), kind=dp), j=1, Ny)] * 2 * pi_dp / Ny
        zlist = [(real((k - 1), kind=dp), k=1, Nz)] * 2 * pi_dp / Nz

        ! test 1: Dirac delta function
        signal(1, :, :, :) = 0._dp
        signal(1, 1, 1, 1) = 1._dp
        fs_exact(1, :, :, :) = cmplx(1._dp, 0._dp, kind=dp)

        ! test 2: onde cosinus 3D
        freqx = 2
        freqy = 1
        freqz = 3
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nz
                    signal(2, i, j, k) = cos(freqx * xlist(i) + freqy * ylist(j) + freqz * zlist(k))
                end do
            end do
        end do
        fs_exact(2, :, :, :) = cmplx(0._dp, 0._dp, kind=dp)
        fs_exact(2, freqx + 1, freqy + 1, freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=dp) / 2, 0._dp, kind=dp)
        fs_exact(2, Nx - freqx + 1, Ny - freqy + 1, Nz - freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=dp) / 2, 0._dp, kind=dp)

        ! test 3: constant signal (DC component)
        signal(3, :, :, :) = 1._dp
        fs_exact(3, :, :, :) = cmplx(0._dp, 0._dp, kind=dp)
        fs_exact(3, 1, 1, 1) = cmplx(real(Nx * Ny * Nz, kind=dp), 0._dp, kind=dp)

        do i = 1, number_tests
            signal_used = signal(i, :, :, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true., num_threads=4)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_dp)
            if (allocated(error)) return
        end do
    end subroutine test_real_dft3_dp

    subroutine test_real_idft3_dp(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: Nx = 4, Ny = 5, Nz = 6, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(dp), dimension(Nx) :: xlist
        real(dp), dimension(Ny) :: ylist
        real(dp), dimension(Nz) :: zlist
        real(dp), dimension(number_tests, Nx, Ny, Nz) :: fs
        real(dp), dimension(Nx, Ny, Nz) :: fs_used
        complex(dp), dimension(number_tests, Nx, Ny, Nz) :: signal_exact
        complex(dp), dimension(Nx, Ny, Nz) :: signal_IDFT3
        integer :: i, j, k, freqx, freqy, freqz

        xlist = [(real((i - 1), kind=dp), i=1, Nx)] * 2 * pi_dp / Nx
        ylist = [(real((j - 1), kind=dp), j=1, Ny)] * 2 * pi_dp / Ny
        zlist = [(real((k - 1), kind=dp), k=1, Nz)] * 2 * pi_dp / Nz

        ! test 1: Dirac delta function
        signal_exact(1, :, :, :) = 0._dp
        signal_exact(1, 1, 1, 1) = 1._dp
        fs(1, :, :, :) = cmplx(1._dp, 0._dp, kind=dp)

        ! test 2: onde cosinus 3D
        freqx = 2
        freqy = 1
        freqz = 3
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nz
                    signal_exact(2, i, j, k) = cos(freqx * xlist(i) + freqy * ylist(j) + freqz * zlist(k))
                end do
            end do
        end do
        fs(2, :, :, :) = cmplx(0._dp, 0._dp, kind=dp)
        fs(2, freqx + 1, freqy + 1, freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=dp) / 2, 0._dp, kind=dp)
        fs(2, Nx - freqx + 1, Ny - freqy + 1, Nz - freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=dp) / 2, 0._dp, kind=dp)

        ! test 3: constant signal (DC component)
        signal_exact(3, :, :, :) = 1._dp
        fs(3, :, :, :) = cmplx(0._dp, 0._dp, kind=dp)
        fs(3, 1, 1, 1) = cmplx(real(Nx * Ny * Nz, kind=dp), 0._dp, kind=dp)

        do i = 1, number_tests
            fs_used = fs(i, :, :, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_dp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true., num_threads=4)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_dp)
            if (allocated(error)) return
        end do
    end subroutine test_real_idft3_dp

    subroutine test_real_dft3_qp(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: Nx = 4, Ny = 5, Nz = 6, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(qp), dimension(Nx) :: xlist
        real(qp), dimension(Ny) :: ylist
        real(qp), dimension(Nz) :: zlist
        real(qp), dimension(number_tests, Nx, Ny, Nz) :: signal
        real(qp), dimension(Nx, Ny, Nz) :: signal_used
        complex(qp), dimension(number_tests, Nx, Ny, Nz) :: fs_exact
        complex(qp), dimension(Nx, Ny, Nz) :: fs_DFT3
        integer :: i, j, k, freqx, freqy, freqz

        xlist = [(real((i - 1), kind=qp), i=1, Nx)] * 2 * pi_qp / Nx
        ylist = [(real((j - 1), kind=qp), j=1, Ny)] * 2 * pi_qp / Ny
        zlist = [(real((k - 1), kind=qp), k=1, Nz)] * 2 * pi_qp / Nz

        ! test 1: Dirac delta function
        signal(1, :, :, :) = 0._qp
        signal(1, 1, 1, 1) = 1._qp
        fs_exact(1, :, :, :) = cmplx(1._qp, 0._qp, kind=qp)

        ! test 2: onde cosinus 3D
        freqx = 2
        freqy = 1
        freqz = 3
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nz
                    signal(2, i, j, k) = cos(freqx * xlist(i) + freqy * ylist(j) + freqz * zlist(k))
                end do
            end do
        end do
        fs_exact(2, :, :, :) = cmplx(0._qp, 0._qp, kind=qp)
        fs_exact(2, freqx + 1, freqy + 1, freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=qp) / 2, 0._qp, kind=qp)
        fs_exact(2, Nx - freqx + 1, Ny - freqy + 1, Nz - freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=qp) / 2, 0._qp, kind=qp)

        ! test 3: constant signal (DC component)
        signal(3, :, :, :) = 1._qp
        fs_exact(3, :, :, :) = cmplx(0._qp, 0._qp, kind=qp)
        fs_exact(3, 1, 1, 1) = cmplx(real(Nx * Ny * Nz, kind=qp), 0._qp, kind=qp)

        do i = 1, number_tests
            signal_used = signal(i, :, :, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true., num_threads=4)
            fs_DFT3 = FT%dft(signal_used, loop_method=loop_method)
            call check(error, maxval(abs(fs_exact(i, :, :, :) - fs_DFT3)) < TOL_TEST_qp)
            if (allocated(error)) return
        end do
    end subroutine test_real_dft3_qp

    subroutine test_real_idft3_qp(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: Nx = 4, Ny = 5, Nz = 6, number_tests = 3
        type(Fourier_Transform) :: FT
        type(LoopMethod) :: loop_method
        real(qp), dimension(Nx) :: xlist
        real(qp), dimension(Ny) :: ylist
        real(qp), dimension(Nz) :: zlist
        real(qp), dimension(number_tests, Nx, Ny, Nz) :: fs
        real(qp), dimension(Nx, Ny, Nz) :: fs_used
        complex(qp), dimension(number_tests, Nx, Ny, Nz) :: signal_exact
        complex(qp), dimension(Nx, Ny, Nz) :: signal_IDFT3
        integer :: i, j, k, freqx, freqy, freqz

        xlist = [(real((i - 1), kind=qp), i=1, Nx)] * 2 * pi_qp / Nx
        ylist = [(real((j - 1), kind=qp), j=1, Ny)] * 2 * pi_qp / Ny
        zlist = [(real((k - 1), kind=qp), k=1, Nz)] * 2 * pi_qp / Nz

        ! test 1: Dirac delta function
        signal_exact(1, :, :, :) = 0._qp
        signal_exact(1, 1, 1, 1) = 1._qp
        fs(1, :, :, :) = cmplx(1._qp, 0._qp, kind=qp)

        ! test 2: onde cosinus 3D
        freqx = 2
        freqy = 1
        freqz = 3
        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nz
                    signal_exact(2, i, j, k) = cos(freqx * xlist(i) + freqy * ylist(j) + freqz * zlist(k))
                end do
            end do
        end do
        fs(2, :, :, :) = cmplx(0._qp, 0._qp, kind=qp)
        fs(2, freqx + 1, freqy + 1, freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=qp) / 2, 0._qp, kind=qp)
        fs(2, Nx - freqx + 1, Ny - freqy + 1, Nz - freqz + 1) = cmplx(real(Nx * Ny * Nz, kind=qp) / 2, 0._qp, kind=qp)

        ! test 3: constant signal (DC component)
        signal_exact(3, :, :, :) = 1._qp
        fs(3, :, :, :) = cmplx(0._qp, 0._qp, kind=qp)
        fs(3, 1, 1, 1) = cmplx(real(Nx * Ny * Nz, kind=qp), 0._qp, kind=qp)

        do i = 1, number_tests
            fs_used = fs(i, :, :, :)
            loop_method = init_loop_method(use_do_classic=.true.)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_vectorized=.true.)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_do_concurrent=.true.)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_qp)
            if (allocated(error)) return

            loop_method = init_loop_method(use_openmp=.true., num_threads=4)
            signal_IDFT3 = FT%idft(fs_used, loop_method=loop_method)
            call check(error, maxval(abs(signal_exact(i, :, :, :) - signal_IDFT3)) < TOL_TEST_qp)
            if (allocated(error)) return
        end do
    end subroutine test_real_idft3_qp

end module test_Fourier_Transform_dft3_real

program test
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_Fourier_Transform_dft3_real

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("Fourier Transform DFT 3D real", collect_Fourier_Transform_dft3) &
                 ]

    call init_color_output(.true.)

    do is = 1, size(testsuites)
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

    deallocate(testsuites)

end program test
