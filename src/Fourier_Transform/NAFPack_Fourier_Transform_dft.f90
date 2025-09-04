submodule(NAFPack_Fourier_Transform) NAFPack_Fourier_Transform_dft

contains

    module function dft_cmplx_sp(signal, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(size(signal)) :: result
        real(sp), dimension(size(signal)) :: n_normalised
        real(sp), dimension(:), allocatable :: k_vector
        complex(sp), dimension(:, :), allocatable :: exp_matrix
        type(LoopMethod) :: loop_method_used
        integer :: N, i, k, j, stat_alloc

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        N = size(signal)

        if (N == 1) then
            result = signal
        else
            n_normalised = [(real(j - 1, sp), j=1, N)] / real(N, sp)
            if (loop_method_used%use_do_classic) then
                do i = 1, N
                    k = i - 1
                    result(i) = sum(signal * exp(im_sp * (-2 * pi_sp * k * n_normalised)))
                end do
            else if (loop_method_used%use_vectorized) then
                allocate (exp_matrix(N, N), stat=stat_alloc)
                if (stat_alloc /= 0) error stop "Error in allocation of exp_matrix in dft_cmplx_dp"

                allocate (k_vector(N), stat=stat_alloc)
                if (stat_alloc /= 0) error stop "Error in allocation of k_vector in dft_cmplx_dp"

                k_vector = [(real(i - 1, sp), i=1, N)]
                exp_matrix = exp(im_sp * reshape(-2 * pi_sp * spread(k_vector, 2, N) * spread(n_normalised, 1, N), [N, N]))
                result = matmul(exp_matrix, signal)

                deallocate (exp_matrix, k_vector)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(i=1:N)
                    k = i - 1
                    result(i) = sum(signal * exp(im_sp * (-2 * pi_sp * k * n_normalised)))
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel do default(none) private(k, i) &
                !$omp& shared(signal, n_normalised, result, N) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)
                do i = 1, N
                    k = i - 1
                    result(i) = sum(signal * exp(im_sp * (-2 * pi_sp * k * n_normalised)))
                end do
                !$omp end parallel do
            end if
        end if

    end function dft_cmplx_sp

    module function dft_cmplx_dp(signal, loop_method) result(result)
        complex(dp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(size(signal)) :: result
        real(dp), dimension(size(signal)) :: n_normalised
        real(dp), dimension(:), allocatable :: k_vector
        complex(dp), dimension(:, :), allocatable :: exp_matrix
        type(LoopMethod) :: loop_method_used
        integer :: N, i, k, j, stat_alloc

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        N = size(signal)

        if (N == 1) then
            result = signal
        else
            n_normalised = [(real(j - 1, dp), j=1, N)] / real(N, dp)
            if (loop_method_used%use_do_classic) then
                do i = 1, N
                    k = i - 1
                    result(i) = sum(signal * exp(im_dp * (-2 * pi_dp * k * n_normalised)))
                end do
            else if (loop_method_used%use_vectorized) then
                allocate (exp_matrix(N, N), stat=stat_alloc)
                if (stat_alloc /= 0) error stop "Error in allocation of exp_matrix in dft_cmplx_dp"

                allocate (k_vector(N), stat=stat_alloc)
                if (stat_alloc /= 0) error stop "Error in allocation of k_vector in dft_cmplx_dp"

                k_vector = [(real(i - 1, dp), i=1, N)]
                exp_matrix = exp(im_dp * reshape(-2 * pi_dp * spread(k_vector, 2, N) * spread(n_normalised, 1, N), [N, N]))
                result = matmul(exp_matrix, signal)

                deallocate (exp_matrix, k_vector)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(i=1:N)
                    k = i - 1
                    result(i) = sum(signal * exp(im_dp * (-2 * pi_dp * k * n_normalised)))
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel do default(none) private(k, i) &
                !$omp& shared(signal, n_normalised, result, N) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)
                do i = 1, N
                    k = i - 1
                    result(i) = sum(signal * exp(im_dp * (-2 * pi_dp * k * n_normalised)))
                end do
                !$omp end parallel do
            end if
        end if

    end function dft_cmplx_dp

    module function dft_cmplx_qp(signal, loop_method) result(result)
        complex(qp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(size(signal)) :: result
        real(qp), dimension(size(signal)) :: n_normalised
        real(qp), dimension(:), allocatable :: k_vector
        complex(qp), dimension(:, :), allocatable :: exp_matrix
        type(LoopMethod) :: loop_method_used
        integer :: N, i, k, j, stat_alloc

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        N = size(signal)

        if (N == 1) then
            result = signal
        else
            n_normalised = [(real(j - 1, qp), j=1, N)] / real(N, qp)
            if (loop_method_used%use_do_classic) then
                do i = 1, N
                    k = i - 1
                    result(i) = sum(signal * exp(im_qp * (-2 * pi_qp * k * n_normalised)))
                end do
            else if (loop_method_used%use_vectorized) then
                allocate (exp_matrix(N, N), stat=stat_alloc)
                if (stat_alloc /= 0) error stop "Error in allocation of exp_matrix in dft_cmplx_dp"

                allocate (k_vector(N), stat=stat_alloc)
                if (stat_alloc /= 0) error stop "Error in allocation of k_vector in dft_cmplx_dp"

                k_vector = [(real(i - 1, qp), i=1, N)]
                exp_matrix = exp(im_qp * reshape(-2 * pi_qp * spread(k_vector, 2, N) * spread(n_normalised, 1, N), [N, N]))
                result = matmul(exp_matrix, signal)

                deallocate (exp_matrix, k_vector)
            else if (loop_method_used%use_do_concurrent) then
                do concurrent(i=1:N)
                    k = i - 1
                    result(i) = sum(signal * exp(im_qp * (-2 * pi_qp * k * n_normalised)))
                end do
            else if (loop_method_used%parallel%use_openmp) then
                !$omp parallel do default(none) private(k, i) &
                !$omp& shared(signal, n_normalised, result, N) &
                !$omp& num_threads(loop_method_used%parallel%num_threads)
                do i = 1, N
                    k = i - 1
                    result(i) = sum(signal * exp(im_qp * (-2 * pi_qp * k * n_normalised)))
                end do
                !$omp end parallel do
            end if
        end if

    end function dft_cmplx_qp

    module function dft_real_sp(signal, loop_method) result(result)
        real(sp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(size(signal)) :: result
        complex(sp), allocatable :: signal_cmplx(:)

        signal_cmplx = cmplx(signal, 0.0_dp, kind=sp)

        if(present(loop_method)) then
            result = dft_cmplx_sp(signal_cmplx, loop_method)
        else
            result = dft_cmplx_sp(signal_cmplx)
        end if
    end function dft_real_sp

    module function dft_real_dp(signal, loop_method) result(result)
        real(dp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(dp), dimension(size(signal)) :: result
        complex(dp), allocatable :: signal_cmplx(:)

        signal_cmplx = cmplx(signal, 0.0_dp, kind=dp)
        
        if(present(loop_method)) then
            result = dft_cmplx_dp(signal_cmplx, loop_method)
        else
            result = dft_cmplx_dp(signal_cmplx)
        end if
    end function dft_real_dp

    module function dft_real_qp(signal, loop_method) result(result)
        real(qp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(qp), dimension(size(signal)) :: result
        complex(qp), allocatable :: signal_cmplx(:)

        signal_cmplx = cmplx(signal, 0.0_qp, kind=qp)
        
        if(present(loop_method)) then
            result = dft_cmplx_qp(signal_cmplx, loop_method)
        else
            result = dft_cmplx_qp(signal_cmplx)
        end if
    end function dft_real_qp

end submodule NAFPack_Fourier_Transform_dft
