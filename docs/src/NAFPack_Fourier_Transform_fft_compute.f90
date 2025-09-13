submodule(NAFPack_Fourier_Transform:NAFPack_Fourier_Transform_fft) NAFPack_Fourier_Transform_fft_compute

    implicit none(type, external)

contains

    !============================================================================
    ! mixed-radix FFT computation 
    !============================================================================

    module function compute_fft_cmplx_sp(&
        signal, plan, stage, current_block_size, block_size, nb_blocks, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        integer(isp), intent(in) :: stage, current_block_size, block_size, nb_blocks
        type(LoopMethod), intent(in) :: loop_method
        complex(sp), dimension(plan%N) :: result

        if (plan%is_pure_radix2) then
            if (loop_method%use_do_classic) then
                result = compute_do_classic_radix2_cmplx_sp(&
                signal, plan, stage, current_block_size, block_size, nb_blocks)
            else if (loop_method%use_vectorized) then
                result = compute_vectorized_radix2_cmplx_sp(&
                    signal, plan, stage, current_block_size, block_size, nb_blocks)
            else if (loop_method%use_do_concurrent) then
                result = compute_do_concurrent_radix2_cmplx_sp(&
                    signal, plan, stage, current_block_size, block_size, nb_blocks)
            else if (loop_method%parallel%use_openmp) then
                result = compute_openmp_radix2_cmplx_sp(&
                    signal, plan, stage, current_block_size, block_size, nb_blocks, &
                    loop_method%parallel%num_threads)
            end if
        else
            if (loop_method%use_do_classic) then
                result = compute_do_classic_mixed_radix_cmplx_sp(&
                    signal, plan, stage, current_block_size, block_size, nb_blocks)
            else if (loop_method%use_vectorized) then
            else if (loop_method%use_do_concurrent) then
            else if (loop_method%parallel%use_openmp) then
            end if
        end if

    end function compute_fft_cmplx_sp

    pure function compute_do_classic_mixed_radix_cmplx_sp(&
    signal, plan, stage, current_block_size, block_size, nb_blocks) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        integer(isp), intent(in) :: stage, current_block_size, block_size, nb_blocks
        complex(sp), dimension(plan%N) :: result
        complex(sp), dimension(plan%twiddles(stage)%radix) :: result_temp
        integer(isp) :: block_index, base, j, r
        complex(sp) :: wj

        do block_index = 0, nb_blocks - 1
            base = block_index * current_block_size
            do j = 0, block_size - 1

                wj = plan%twiddles(stage)%twiddles_factor(j + 1)
                do r = 0, plan%twiddles(stage)%radix - 1
                    result_temp(r + 1) = signal(base + j + r * block_size + 1) * wj**r
                end do

                result_temp = small_dft_kernel_sp(result_temp, plan, stage)

                do r = 0, plan%twiddles(stage)%radix - 1
                    result(base + j + r * block_size + 1) = result_temp(r + 1)
                end do
            end do
        end do

    end function compute_do_classic_mixed_radix_cmplx_sp

    pure function small_dft_kernel_sp(x, plan, stage) result(y)
        complex(sp), dimension(:), intent(in) :: x
        type(FFTPlan), intent(in) :: plan
        integer(isp), intent(in) :: stage
        complex(sp), dimension(plan%twiddles(stage)%radix) :: y
        complex(sp) :: W3, W4, W5, w_kn
        integer(isp) :: k, n

        W3 = exp(-2.0_sp * pi_sp * im_sp / 3.0_sp)
        W4 = im_sp
        W5 = exp(-2.0_sp * pi_sp * im_sp / 5.0_sp)

        select case (plan%twiddles(stage)%radix)
        case (2)
            y(1) = x(1) + x(2)
            y(2) = x(1) - x(2)
        case (3)
            y(1) = x(1) +         x(2) +         x(3)
            y(2) = x(1) + W3 *    x(2) + W3**2 * x(3)
            y(3) = x(1) + W3**2 * x(2) + W3 *    x(3)
        case (4)
            y(1) = x(1) +      x(2) + x(3) +      x(4)
            y(2) = x(1) - W4 * x(2) - x(3) + W4 * x(4)
            y(3) = x(1) -      x(2) + x(3) -      x(4)
            y(4) = x(1) + W4 * x(2) - x(3) - W4 * x(4)
        case (5)
            y(1) = x(1) +         x(2) +         x(3) +         x(4) +         x(5)
            y(2) = x(1) + W5 *    x(2) + W5**2 * x(3) + W5**3 * x(4) + W5**4 * x(5)
            y(3) = x(1) + W5**2 * x(2) + W5**4 * x(3) + W5    * x(4) + W5**3 * x(5)
            y(4) = x(1) + W5**3 * x(2) + W5    * x(3) + W5**4 * x(4) + W5**2 * x(5)
            y(5) = x(1) + W5**4 * x(2) + W5**3 * x(3) + W5**2 * x(4) + W5    * x(5)
        case default
            y = (0.0_sp, 0.0_sp)
            do k = 1, size(x)
                do n = 1, size(x)
                    w_kn = exp(-2.0_sp * pi_sp * im_sp * real((n-1)*(k-1), sp) / real(size(x), sp))
                    y(k) = y(k) + x(n) * w_kn
                end do
            end do
        end select

    end function small_dft_kernel_sp

    !=============================================================================
    ! radix-2 FFT computation with different loop methods
    !=============================================================================

    pure function compute_do_classic_radix2_cmplx_sp(&
    signal, plan, stage, current_block_size, block_size, nb_blocks) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        integer(isp), intent(in) :: stage, current_block_size, block_size, nb_blocks
        complex(sp), dimension(plan%N) :: result
        integer(isp) :: j, block_index, base
        complex(sp) :: even, odd, twiddle

        do block_index = 0, nb_blocks - 1
            base = block_index * current_block_size
            do j = 0, block_size - 1
                even = signal(base + j + 1)
                odd = signal(base + j + block_size + 1)

                twiddle = plan%twiddles(stage)%twiddles_factor(j + 1)

                result(base + j + 1) = even + twiddle * odd
                result(base + j + block_size + 1) = even - twiddle * odd
            end do
        end do

    end function compute_do_classic_radix2_cmplx_sp

    pure function compute_vectorized_radix2_cmplx_sp(&
    signal, plan, stage, current_block_size, block_size, nb_blocks) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        integer(isp), intent(in) :: stage, current_block_size, block_size, nb_blocks
        complex(sp), dimension(plan%N) :: result
        complex(sp), dimension(:), allocatable :: even, odd, twiddles
        integer(sp), dimension(:), allocatable :: idx
        integer(isp) :: j, block_index, base


        allocate (idx(block_size), twiddles(block_size), odd(block_size), even(block_size))
        idx = [(j, j=0, block_size - 1)]
        twiddles = plan%twiddles(stage)%twiddles_factor(idx + 1)

        do block_index = 0, nb_blocks - 1
            base = block_index * current_block_size

            even = signal(base + idx + 1)
            odd = signal(base + idx + block_size + 1)

            result(base + idx + 1) = even + twiddles * odd
            result(base + idx + block_size + 1) = even - twiddles * odd
        end do

        deallocate (idx, twiddles, odd, even)

    end function compute_vectorized_radix2_cmplx_sp

    pure function compute_do_concurrent_radix2_cmplx_sp(&
    signal, plan, stage, current_block_size, block_size, nb_blocks) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        integer(isp), intent(in) :: stage, current_block_size, block_size, nb_blocks
        complex(sp), dimension(plan%N) :: result
        integer(isp) :: j, base, block_index
        complex(sp) :: even, odd, twiddle

        do concurrent(block_index = 0:nb_blocks - 1, j = 0:block_size - 1)
            base = block_index * current_block_size

            even = signal(base + j + 1)
            odd  = signal(base + j + block_size + 1)
            twiddle = plan%twiddles(stage)%twiddles_factor(j + 1)

            result(base + j + 1) = even + twiddle * odd
            result(base + j + block_size + 1) = even - twiddle * odd
        end do

    end function compute_do_concurrent_radix2_cmplx_sp

    function compute_openmp_radix2_cmplx_sp(&
        signal, plan, stage, current_block_size, block_size, nb_blocks, threads) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        integer(isp), intent(in) :: stage, current_block_size, block_size, nb_blocks, threads
        complex(sp), dimension(plan%N) :: result
        complex(sp) :: even, odd, twiddle
        integer(isp) :: j, block_index, base

        !$omp parallel do default(none) private(block_index, base, j, even, odd, twiddle) &
        !$omp& shared(result, plan, block_size, stage, current_block_size, signal, nb_blocks) &
        !$omp& num_threads(threads)
        do block_index = 0, nb_blocks - 1
            base = block_index * current_block_size
            do j = 0, block_size - 1
                even = signal(base + j + 1)
                odd = signal(base + j + block_size + 1)

                twiddle = plan%twiddles(stage)%twiddles_factor(j + 1)

                result(base + j + 1) = even + twiddle * odd
                result(base + j + block_size + 1) = even - twiddle * odd
            end do
        end do
        !$omp end parallel do

    end function compute_openmp_radix2_cmplx_sp

end submodule NAFPack_Fourier_Transform_fft_compute
