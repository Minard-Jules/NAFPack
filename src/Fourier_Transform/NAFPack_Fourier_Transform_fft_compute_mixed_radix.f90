submodule(NAFPack_Fourier_Transform:NAFPack_Fourier_Transform_fft) NAFPack_Fourier_Transform_fft_compute_mixed_radix

    implicit none(type, external)

contains

    module function compute_fft_mixed_radix_cmplx_sp( &
        signal, plan, stage_params, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params
        type(LoopMethod), intent(in) :: loop_method
        complex(sp), dimension(plan%N) :: result

        if (loop_method%use_do_classic) then
            result = compute_do_classic_cmplx_sp( &
                    signal, plan, stage_params)
        else if (loop_method%use_vectorized) then
            result = compute_vectorized_cmplx_sp( &
                    signal, plan, stage_params)
        else if (loop_method%use_do_concurrent) then
            result = compute_do_concurrent_cmplx_sp( &
                    signal, plan, stage_params)
        else if (loop_method%parallel%use_openmp) then
            result = compute_openmp_cmplx_sp( &
                    signal, plan, stage_params, &
                    loop_method%parallel%num_threads)
        end if

    end function compute_fft_mixed_radix_cmplx_sp

    pure function compute_do_classic_cmplx_sp( &
        signal, plan, stage_params) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params
        complex(sp), dimension(plan%N) :: result
        complex(sp), dimension(plan%twiddles(stage_params%stage)%radix) :: result_temp
        integer(isp) :: block_index, base, j, r, radix
        complex(sp) :: wj

        result = signal

        radix = plan%twiddles(stage_params%stage)%radix

        do block_index = 0, stage_params%nb_blocks - 1
            base = block_index * stage_params%current_block_size
            do j = 0, stage_params%block_size - 1
                select case (plan%algorithm%decimation_method%id)
                case (DIT%id)
                    wj = plan%twiddles(stage_params%stage)%twiddles_factor(j + 1)
                    do r = 0, radix - 1
                        result_temp(r + 1) = result(base + j + r * stage_params%block_size + 1) * wj**r
                    end do

            result_temp = reshape(small_dft_kernel_sp(reshape(result_temp, [radix, 1]), plan, stage_params%stage, radix, 1), [radix])

                    do r = 0, radix - 1
                        result(base + j + r * stage_params%block_size + 1) = result_temp(r + 1)
                    end do
                case (DIF%id)
                    do r = 0, radix - 1
                        result_temp(r + 1) = result(base + j + r * stage_params%block_size + 1)
                    end do

            result_temp = reshape(small_dft_kernel_sp(reshape(result_temp, [radix, 1]), plan, stage_params%stage, radix, 1), [radix])

                    ! * apply twiddle factors after the small DFT (exept for r=0)
                    wj = plan%twiddles(stage_params%stage)%twiddles_factor(j + 1)
                    do r = 1, radix - 1
                        result_temp(r + 1) = result_temp(r + 1) * wj**r
                    end do

                    do r = 0, radix - 1
                        result(base + j + r * stage_params%block_size + 1) = result_temp(r + 1)
                    end do
                end select
            end do
        end do

    end function compute_do_classic_cmplx_sp

    pure function compute_vectorized_cmplx_sp( &
        signal, plan, stage_params) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params
        complex(sp), dimension(plan%N) :: result
        complex(sp), dimension(plan%twiddles(stage_params%stage)%radix, stage_params%block_size) :: result_temp
        complex(sp), dimension(plan%twiddles(stage_params%stage)%radix, stage_params%block_size) :: fourier_matrix
        complex(sp), dimension(stage_params%block_size) :: twiddles
        integer(isp), dimension(plan%twiddles(stage_params%stage)%radix*stage_params%block_size) :: indices
        integer(isp), dimension(plan%twiddles(stage_params%stage)%radix) :: r_idx
        integer(isp), dimension(stage_params%block_size) :: idx
        integer(isp) :: j, block_index, base

        result = signal

        idx = [(j, j=0, stage_params%block_size - 1)]
        r_idx = [(j, j=0, plan%twiddles(stage_params%stage)%radix - 1)]
        indices = reshape(spread(idx, dim=1, ncopies=plan%twiddles(stage_params%stage)%radix) + &
                        spread(r_idx * stage_params%block_size, dim=2, ncopies=stage_params%block_size), &
                        [plan%twiddles(stage_params%stage)%radix * stage_params%block_size])

        twiddles = plan%twiddles(stage_params%stage)%twiddles_factor(idx + 1)
        fourier_matrix = spread(twiddles, 1, plan%twiddles(stage_params%stage)%radix)**spread(r_idx, 2, stage_params%block_size)

        do block_index = 0, stage_params%nb_blocks - 1

            select case (plan%algorithm%decimation_method%id)
            case (DIT%id)
                base = block_index * stage_params%current_block_size

                result_temp = reshape(result(base + indices + 1), [plan%twiddles(stage_params%stage)%radix, stage_params%block_size])
                result_temp = fourier_matrix * result_temp

                    result_temp = small_dft_kernel_sp(result_temp, plan, stage_params%stage, plan%twiddles(stage_params%stage)%radix, stage_params%block_size)

                result(base + indices + 1) = reshape(result_temp, [plan%twiddles(stage_params%stage)%radix * stage_params%block_size])
            case (DIF%id)
                base = block_index * stage_params%current_block_size

                result_temp = reshape(result(base + indices + 1), [plan%twiddles(stage_params%stage)%radix, stage_params%block_size])

                    result_temp = small_dft_kernel_sp(result_temp, plan, stage_params%stage, plan%twiddles(stage_params%stage)%radix, stage_params%block_size)

                ! * apply twiddle factors after the small DFT (exept for r=0)
                fourier_matrix(1, :) = 1.0_sp
                result_temp = fourier_matrix * result_temp

                result(base + indices + 1) = reshape(result_temp, [plan%twiddles(stage_params%stage)%radix * stage_params%block_size])
            end select
        end do

    end function compute_vectorized_cmplx_sp

    pure function compute_do_concurrent_cmplx_sp( &
        signal, plan, stage_params) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params
        complex(sp), dimension(plan%N) :: result
        complex(sp), dimension(plan%twiddles(stage_params%stage)%radix) :: result_temp
        integer(isp) :: block_index, base, j, r, radix
        complex(sp) :: wj

        result = signal

        radix = plan%twiddles(stage_params%stage)%radix

        do concurrent(block_index=0:stage_params%nb_blocks - 1, j=0:stage_params%block_size - 1)
            base = block_index * stage_params%current_block_size

            select case (plan%algorithm%decimation_method%id)
            case (DIT%id)
                wj = plan%twiddles(stage_params%stage)%twiddles_factor(j + 1)
                do r = 0, radix - 1
                    result_temp(r + 1) = result(base + j + r * stage_params%block_size + 1) * wj**r
                end do

            result_temp = reshape(small_dft_kernel_sp(reshape(result_temp, [radix, 1]), plan, stage_params%stage, radix, 1), [radix])

                do r = 0, radix - 1
                    result(base + j + r * stage_params%block_size + 1) = result_temp(r + 1)
                end do
            case (DIF%id)
                do r = 0, radix - 1
                    result_temp(r + 1) = result(base + j + r * stage_params%block_size + 1)
                end do

            result_temp = reshape(small_dft_kernel_sp(reshape(result_temp, [radix, 1]), plan, stage_params%stage, radix, 1), [radix])

                ! * apply twiddle factors after the small DFT (exept for r=0)
                wj = plan%twiddles(stage_params%stage)%twiddles_factor(j + 1)
                do r = 1, radix - 1
                    result_temp(r + 1) = result_temp(r + 1) * wj**r
                end do

                do r = 0, radix - 1
                    result(base + j + r * stage_params%block_size + 1) = result_temp(r + 1)
                end do
            end select
        end do

    end function compute_do_concurrent_cmplx_sp

    function compute_openmp_cmplx_sp( &
        signal, plan, stage_params, threads) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params
        integer(isp), intent(in) :: threads
        complex(sp), dimension(plan%N) :: result
        complex(sp), dimension(plan%twiddles(stage_params%stage)%radix) :: result_temp
        integer(isp) :: block_index, base, j, r, radix
        complex(sp) :: wj

        result = signal

        radix = plan%twiddles(stage_params%stage)%radix

        !$omp parallel do default(none) private(block_index, base, j, r, wj, result_temp) &
        !$omp& shared(result, plan, stage_params, radix) &
        !$omp& num_threads(threads)
        do block_index = 0, stage_params%nb_blocks - 1
            base = block_index * stage_params%current_block_size
            do j = 0, stage_params%block_size - 1
                select case (plan%algorithm%decimation_method%id)
                case (DIT%id)
                    wj = plan%twiddles(stage_params%stage)%twiddles_factor(j + 1)
                    do r = 0, radix - 1
                        result_temp(r + 1) = result(base + j + r * stage_params%block_size + 1) * wj**r
                    end do

            result_temp = reshape(small_dft_kernel_sp(reshape(result_temp, [radix, 1]), plan, stage_params%stage, radix, 1), [radix])

                    do r = 0, radix - 1
                        result(base + j + r * stage_params%block_size + 1) = result_temp(r + 1)
                    end do
                case (DIF%id)
                    do r = 0, radix - 1
                        result_temp(r + 1) = result(base + j + r * stage_params%block_size + 1)
                    end do

            result_temp = reshape(small_dft_kernel_sp(reshape(result_temp, [radix, 1]), plan, stage_params%stage, radix, 1), [radix])

                    ! * apply twiddle factors after the small DFT (exept for r=0)
                    wj = plan%twiddles(stage_params%stage)%twiddles_factor(j + 1)
                    do r = 1, radix - 1
                        result_temp(r + 1) = result_temp(r + 1) * wj**r
                    end do

                    do r = 0, radix - 1
                        result(base + j + r * stage_params%block_size + 1) = result_temp(r + 1)
                    end do
                end select
            end do
        end do
        !$omp end parallel do

    end function compute_openmp_cmplx_sp

    pure function small_dft_kernel_sp(x, plan, stage, radix, block_size) result(y)
        complex(sp), dimension(:, :), intent(in) :: x
        type(FFTPlan), intent(in) :: plan
        integer(isp), intent(in) :: stage, radix, block_size
        complex(sp), dimension(radix, block_size) :: y
        complex(sp), dimension(:, :), allocatable :: w_kn
        real(sp), dimension(:, :), allocatable :: kvec, nvec
        complex(sp) :: W3, W4, W5
        integer(isp) :: k, n

        W3 = exp(-2.0_sp * pi_sp * im_sp / 3.0_sp)
        W4 = im_sp
        W5 = exp(-2.0_sp * pi_sp * im_sp / 5.0_sp)

        select case (plan%twiddles(stage)%radix)
        case (2)
            y(1, :) = x(1, :) + x(2, :)
            y(2, :) = x(1, :) - x(2, :)
        case (3)
            y(1, :) = x(1, :) + x(2, :) + x(3, :)
            y(2, :) = x(1, :) + W3 * x(2, :) + W3**2 * x(3, :)
            y(3, :) = x(1, :) + W3**2 * x(2, :) + W3 * x(3, :)
        case (4)
            y(1, :) = x(1, :) + x(2, :) + x(3, :) + x(4, :)
            y(2, :) = x(1, :) - W4 * x(2, :) - x(3, :) + W4 * x(4, :)
            y(3, :) = x(1, :) - x(2, :) + x(3, :) - x(4, :)
            y(4, :) = x(1, :) + W4 * x(2, :) - x(3, :) - W4 * x(4, :)
        case (5)
            y(1, :) = x(1, :) + x(2, :) + x(3, :) + x(4, :) + x(5, :)
            y(2, :) = x(1, :) + W5 * x(2, :) + W5**2 * x(3, :) + W5**3 * x(4, :) + W5**4 * x(5, :)
            y(3, :) = x(1, :) + W5**2 * x(2, :) + W5**4 * x(3, :) + W5 * x(4, :) + W5**3 * x(5, :)
            y(4, :) = x(1, :) + W5**3 * x(2, :) + W5 * x(3, :) + W5**4 * x(4, :) + W5**2 * x(5, :)
            y(5, :) = x(1, :) + W5**4 * x(2, :) + W5**3 * x(3, :) + W5**2 * x(4, :) + W5 * x(5, :)
        case default
            allocate (w_kn(radix, radix), kvec(radix, radix), nvec(radix, radix))
            kvec = spread([(k - 1, k=1, radix)], dim=1, ncopies=radix)
            nvec = spread([(n - 1, n=1, radix)], dim=2, ncopies=radix)
            w_kn = exp(-2.0_sp * pi_sp * im_sp * kvec * nvec / real(radix, sp))
            y = matmul(w_kn, x)
        end select

    end function small_dft_kernel_sp

end submodule NAFPack_Fourier_Transform_fft_compute_mixed_radix
