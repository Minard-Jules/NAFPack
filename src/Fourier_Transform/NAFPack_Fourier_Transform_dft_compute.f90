submodule(NAFPack_Fourier_Transform:NAFPack_Fourier_Transform_dft) NAFPack_Fourier_Transform_dft_compute

    implicit none(type, external)

contains

    !=================================================================================
    ! Compute the discrete Fourier transform of a complex signal in simple precision
    !=================================================================================

    module function compute_dft_cmplx_sp(signal, N, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        integer(isp), intent(in) :: N
        type(LoopMethod), intent(in) :: loop_method
        complex(sp), dimension(N) :: result
        real(sp), dimension(size(signal)) :: n_vec
        complex(sp) :: omega
        integer(isp) :: i

        if (N == 1) then
            result = signal
        else
            omega = exp(-2 * pi_sp * im_sp / real(N, sp))
            n_vec = [(real(i - 1, sp), i=1, N)]
            if (loop_method%use_do_classic) then
                result = compute_do_classic_cmplx_sp(signal, n_vec, omega, N)
            else if (loop_method%use_vectorized) then
                result = compute_do_vectorized_cmplx_sp(signal, n_vec, omega, N)
            else if (loop_method%use_do_concurrent) then
                result = compute_do_concurrent_cmplx_sp(signal, n_vec, omega, N)
            else if (loop_method%parallel%use_openmp) then
                result = compute_openmp_cmplx_sp(signal, n_vec, omega, N, &
                                                    loop_method%parallel%num_threads)
            end if
        end if

    end function compute_dft_cmplx_sp

    pure function compute_do_classic_cmplx_sp(signal, n_vector, omega, N) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        real(sp), dimension(:), intent(in) :: n_vector
        complex(sp), intent(in) :: omega
        integer(isp), intent(in) :: N
        complex(sp), dimension(N) :: result
        integer(isp) :: i, k

        do i = 1, N
            k = i - 1
            result(i) = sum(signal * omega**(k * n_vector))
        end do

    end function compute_do_classic_cmplx_sp

    pure function compute_do_vectorized_cmplx_sp(signal, n_vector, omega, N) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        real(sp), dimension(:), intent(in) :: n_vector
        complex(sp), intent(in) :: omega
        integer(isp), intent(in) :: N
        complex(sp), dimension(N) :: result
        real(sp), dimension(:), allocatable :: k_vector
        complex(sp), dimension(:, :), allocatable :: fourier_matrix
        integer(isp) :: stat_alloc

        allocate (fourier_matrix(N, N), stat=stat_alloc)
        if (stat_alloc /= 0) error stop "Error in allocation of fourier_matrix in dft_cmplx_sp"
        allocate (k_vector(N), stat=stat_alloc)
        if (stat_alloc /= 0) error stop "Error in allocation of k_vector in dft_cmplx_sp"

        k_vector = n_vector
        fourier_matrix = omega**reshape(spread(k_vector, 2, N) * spread(n_vector, 1, N), [N, N])
        result = matmul(fourier_matrix, signal)

        deallocate (fourier_matrix, k_vector)

    end function compute_do_vectorized_cmplx_sp

    pure function compute_do_concurrent_cmplx_sp(signal, n_vector, omega, N) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        real(sp), dimension(:), intent(in) :: n_vector
        complex(sp), intent(in) :: omega
        integer(isp), intent(in) :: N
        complex(sp), dimension(N) :: result
        integer(isp) :: i, k

        do concurrent(i=1:N)
            k = i - 1
            result(i) = sum(signal * omega**(k * n_vector))
        end do

    end function compute_do_concurrent_cmplx_sp

    function compute_openmp_cmplx_sp(signal, n_vector, omega, N, threads) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        real(sp), dimension(:), intent(in) :: n_vector
        complex(sp), intent(in) :: omega
        integer(isp), intent(in) :: N, threads
        complex(sp), dimension(N) :: result
        integer(isp) :: i, k

        !$omp parallel do default(none) private(k, i) &
        !$omp& shared(signal, n_vector, omega, result, N) &
        !$omp& num_threads(threads)
        do i = 1, N
            k = i - 1
            result(i) = sum(signal * omega**(k * n_vector))
        end do
        !$omp end parallel do

    end function compute_openmp_cmplx_sp

    !=================================================================================
    ! Compute the discrete Fourier transform of a complex signal in double precision
    !=================================================================================

    module function compute_dft_cmplx_dp(signal, N, loop_method) result(result)
        complex(dp), dimension(:), intent(in) :: signal
        integer(isp), intent(in) :: N
        type(LoopMethod), intent(in) :: loop_method
        complex(dp), dimension(N) :: result
        real(dp), dimension(size(signal)) :: n_vec
        complex(dp) :: omega
        integer(isp) :: i

        if (N == 1) then
            result = signal
        else
            omega = exp(-2 * pi_dp * im_dp / real(N, dp))
            n_vec = [(real(i - 1, dp), i=1, N)]
            if (loop_method%use_do_classic) then
                result = compute_do_classic_cmplx_dp(signal, n_vec, omega, N)
            else if (loop_method%use_vectorized) then
                result = compute_do_vectorized_cmplx_dp(signal, n_vec, omega, N)
            else if (loop_method%use_do_concurrent) then
                result = compute_do_concurrent_cmplx_dp(signal, n_vec, omega, N)
            else if (loop_method%parallel%use_openmp) then
                result = compute_openmp_cmplx_dp(signal, n_vec, omega, N, &
                                                    loop_method%parallel%num_threads)
            end if
        end if

    end function compute_dft_cmplx_dp

    module pure function compute_do_classic_cmplx_dp(signal, n_vector, omega, N) result(result)
        complex(dp), dimension(:), intent(in) :: signal
        real(dp), dimension(:), intent(in) :: n_vector
        complex(dp), intent(in) :: omega
        integer(isp), intent(in) :: N
        complex(dp), dimension(N) :: result
        integer(isp) :: i, k

        do i = 1, N
            k = i - 1
            result(i) = sum(signal * omega**(k * n_vector))
        end do

    end function compute_do_classic_cmplx_dp

    module pure function compute_do_vectorized_cmplx_dp(signal, n_vector, omega, N) result(result)
        complex(dp), dimension(:), intent(in) :: signal
        real(dp), dimension(:), intent(in) :: n_vector
        complex(dp), intent(in) :: omega
        integer(isp), intent(in) :: N
        complex(dp), dimension(N) :: result
        real(dp), dimension(:), allocatable :: k_vector
        complex(dp), dimension(:, :), allocatable :: fourier_matrix
        integer(isp) :: stat_alloc

        allocate (fourier_matrix(N, N), stat=stat_alloc)
        if (stat_alloc /= 0) error stop "Error in allocation of fourier_matrix in dft_cmplx_sp"
        allocate (k_vector(N), stat=stat_alloc)
        if (stat_alloc /= 0) error stop "Error in allocation of k_vector in dft_cmplx_sp"

        k_vector = n_vector
        fourier_matrix = omega**reshape(spread(k_vector, 2, N) * spread(n_vector, 1, N), [N, N])
        result = matmul(fourier_matrix, signal)

        deallocate (fourier_matrix, k_vector)

    end function compute_do_vectorized_cmplx_dp

    module pure function compute_do_concurrent_cmplx_dp(signal, n_vector, omega, N) result(result)
        complex(dp), dimension(:), intent(in) :: signal
        real(dp), dimension(:), intent(in) :: n_vector
        complex(dp), intent(in) :: omega
        integer(isp), intent(in) :: N
        complex(dp), dimension(N) :: result
        integer(isp) :: i, k

        do concurrent(i=1:N)
            k = i - 1
            result(i) = sum(signal * omega**(k * n_vector))
        end do

    end function compute_do_concurrent_cmplx_dp

    module function compute_openmp_cmplx_dp(signal, n_vector, omega, N, threads) result(result)
        complex(dp), dimension(:), intent(in) :: signal
        real(dp), dimension(:), intent(in) :: n_vector
        complex(dp), intent(in) :: omega
        integer(isp), intent(in) :: N, threads
        complex(dp), dimension(N) :: result
        integer(isp) :: i, k

        !$omp parallel do default(none) private(k, i) &
        !$omp& shared(signal, n_vector, omega, result, N) &
        !$omp& num_threads(threads)
        do i = 1, N
            k = i - 1
            result(i) = sum(signal * omega**(k * n_vector))
        end do
        !$omp end parallel do

    end function compute_openmp_cmplx_dp

    !=================================================================================
    ! Compute the discrete Fourier transform of a complex signal in quadruple precision
    !=================================================================================

    module function compute_dft_cmplx_qp(signal, N, loop_method) result(result)
        complex(qp), dimension(:), intent(in) :: signal
        integer(isp), intent(in) :: N
        type(LoopMethod), intent(in) :: loop_method
        complex(qp), dimension(N) :: result
        real(qp), dimension(size(signal)) :: n_vec
        complex(qp) :: omega
        integer(isp) :: i

        if (N == 1) then
            result = signal
        else
            omega = exp(-2 * pi_qp * im_qp / real(N, qp))
            n_vec = [(real(i - 1, qp), i=1, N)]
            if (loop_method%use_do_classic) then
                result = compute_do_classic_cmplx_qp(signal, n_vec, omega, N)
            else if (loop_method%use_vectorized) then
                result = compute_do_vectorized_cmplx_qp(signal, n_vec, omega, N)
            else if (loop_method%use_do_concurrent) then
                result = compute_do_concurrent_cmplx_qp(signal, n_vec, omega, N)
            else if (loop_method%parallel%use_openmp) then
                result = compute_openmp_cmplx_qp(signal, n_vec, omega, N, &
                                                    loop_method%parallel%num_threads)
            end if
        end if

    end function compute_dft_cmplx_qp

    module pure function compute_do_classic_cmplx_qp(signal, n_vector, omega, N) result(result)
        complex(qp), dimension(:), intent(in) :: signal
        real(qp), dimension(:), intent(in) :: n_vector
        complex(qp), intent(in) :: omega
        integer(isp), intent(in) :: N
        complex(qp), dimension(N) :: result
        integer(isp) :: i, k

        do i = 1, N
            k = i - 1
            result(i) = sum(signal * omega**(k * n_vector))
        end do

    end function compute_do_classic_cmplx_qp

    module pure function compute_do_vectorized_cmplx_qp(signal, n_vector, omega, N) result(result)
        complex(qp), dimension(:), intent(in) :: signal
        real(qp), dimension(:), intent(in) :: n_vector
        complex(qp), intent(in) :: omega
        integer(isp), intent(in) :: N
        complex(qp), dimension(N) :: result
        real(qp), dimension(:), allocatable :: k_vector
        complex(qp), dimension(:, :), allocatable :: fourier_matrix
        integer(isp) :: stat_alloc

        allocate (fourier_matrix(N, N), stat=stat_alloc)
        if (stat_alloc /= 0) error stop "Error in allocation of fourier_matrix in dft_cmplx_sp"
        allocate (k_vector(N), stat=stat_alloc)
        if (stat_alloc /= 0) error stop "Error in allocation of k_vector in dft_cmplx_sp"

        k_vector = n_vector
        fourier_matrix = omega**reshape(spread(k_vector, 2, N) * spread(n_vector, 1, N), [N, N])
        result = matmul(fourier_matrix, signal)

        deallocate (fourier_matrix, k_vector)

    end function compute_do_vectorized_cmplx_qp

    module pure function compute_do_concurrent_cmplx_qp(signal, n_vector, omega, N) result(result)
        complex(qp), dimension(:), intent(in) :: signal
        real(qp), dimension(:), intent(in) :: n_vector
        complex(qp), intent(in) :: omega
        integer(isp), intent(in) :: N
        complex(qp), dimension(N) :: result
        integer(isp) :: i, k

        do concurrent(i=1:N)
            k = i - 1
            result(i) = sum(signal * omega**(k * n_vector))
        end do

    end function compute_do_concurrent_cmplx_qp

    module function compute_openmp_cmplx_qp(signal, n_vector, omega, N, threads) result(result)
        complex(qp), dimension(:), intent(in) :: signal
        real(qp), dimension(:), intent(in) :: n_vector
        complex(qp), intent(in) :: omega
        integer(isp), intent(in) :: N, threads
        complex(qp), dimension(N) :: result
        integer(isp) :: i, k

        !$omp parallel do default(none) private(k, i) &
        !$omp& shared(signal, n_vector, omega, result, N) &
        !$omp& num_threads(threads)
        do i = 1, N
            k = i - 1
            result(i) = sum(signal * omega**(k * n_vector))
        end do
        !$omp end parallel do

    end function compute_openmp_cmplx_qp

end submodule NAFPack_Fourier_Transform_dft_compute
