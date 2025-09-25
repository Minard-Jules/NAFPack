submodule(NAFPack_Fourier_Transform:NAFPack_Fourier_Transform_fft) NAFPack_Fourier_Transform_fft_compute_radix2

implicit none(type, external)

contains

recursive module function compute_fft_radix2_recursive_cmplx_sp( &
    signal, plan, stage) result(result)
    complex(sp), dimension(:), intent(in) :: signal
    type(FFTPlan), intent(in) :: plan
    integer(isp), intent(in) :: stage
    complex(sp), dimension(size(signal)) :: result
    complex(sp), dimension(size(signal)/2) :: twiddle, even, odd
    integer(isp) :: N, i0, i1

    N = size(signal)
    i0 = N / 2
    i1 = i0 + 1

    if (N == 1) then
        result = signal
    else if (N == 2) then
        result(1) = signal(1) + signal(2)
        result(2) = signal(1) - signal(2)
        return
    end if

    if (iand(N, N - 1) /= 0) error stop "N must be power of 2"

    twiddle = plan%twiddles(stage)%twiddles_factor
    select case (plan%algorithm%decimation_method%id)
    case (DIT%id)
        even = compute_fft_radix2_recursive_cmplx_sp(signal(1:N:2), plan, stage - 1)
        odd = compute_fft_radix2_recursive_cmplx_sp(signal(2:N:2), plan, stage - 1)

        result(:i0) = even + odd * twiddle
        result(i1:) = even - odd * twiddle
    case (DIF%id)
        even = signal(:i0) + signal(i1:)
        odd = (signal(:i0) - signal(i1:)) * twiddle

        even = compute_fft_radix2_recursive_cmplx_sp(even, plan, stage + 1)
        odd = compute_fft_radix2_recursive_cmplx_sp(odd, plan, stage + 1)

        result(1:N:2) = even
        result(2:N:2) = odd
    end select

end function compute_fft_radix2_recursive_cmplx_sp

module function compute_fft_radix2_iterative_cmplx_sp( &
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

end function compute_fft_radix2_iterative_cmplx_sp

pure function compute_do_classic_cmplx_sp( &
    signal, plan, stage_params) result(result)
    complex(sp), dimension(:), intent(in) :: signal
    type(FFTPlan), intent(in) :: plan
    type(FFTStageParams), intent(in) :: stage_params
    complex(sp), dimension(plan%N) :: result
    complex(sp) :: even, odd, twiddle
    integer(isp) :: j, block_index, base
    integer(isp) :: i0, i1

    do block_index = 0, stage_params%nb_blocks - 1
        base = block_index * stage_params%current_block_size
        do j = 0, stage_params%block_size - 1
            i0 = base + j + 1
            i1 = i0 + stage_params%block_size
            even = signal(i0)
            odd = signal(i1)

            twiddle = plan%twiddles(stage_params%stage)%twiddles_factor(j + 1)
            select case (plan%algorithm%decimation_method%id)
            case (DIT%id)
                result(i0) = even + twiddle * odd
                result(i1) = even - twiddle * odd
            case (DIF%id)
                result(i0) = even + odd
                result(i1) = (even - odd) * twiddle
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
    complex(sp), dimension(stage_params%block_size) :: even, odd, twiddles
    integer(sp), dimension(stage_params%block_size) :: idx, i0, i1
    integer(isp) :: j, block_index, base

    result = signal

    idx = [(j, j=0, stage_params%block_size - 1)]
    twiddles = plan%twiddles(stage_params%stage)%twiddles_factor(idx + 1)

    do block_index = 0, stage_params%nb_blocks - 1
        base = block_index * stage_params%current_block_size
        i0 = base + idx + 1
        i1 = i0 + stage_params%block_size

        even = result(i0)
        odd = result(i1)

        select case (plan%algorithm%decimation_method%id)
        case (DIT%id)
            result(i0) = even + twiddles * odd
            result(i1) = even - twiddles * odd
        case (DIF%id)
            result(i0) = even + odd
            result(i1) = (even - odd) * twiddles
        end select
    end do

end function compute_vectorized_cmplx_sp

pure function compute_do_concurrent_cmplx_sp( &
    signal, plan, stage_params) result(result)
    complex(sp), dimension(:), intent(in) :: signal
    type(FFTPlan), intent(in) :: plan
    type(FFTStageParams), intent(in) :: stage_params
    complex(sp), dimension(plan%N) :: result
    complex(sp) :: even, odd, twiddle
    integer(isp) :: j, base, block_index
    integer(isp) :: i0, i1

    result = signal

    do concurrent(block_index=0:stage_params%nb_blocks - 1, j=0:stage_params%block_size - 1)
        base = block_index * stage_params%current_block_size
        i0 = base + j + 1
        i1 = i0 + stage_params%block_size

        even = result(i0)
        odd = result(i1)
        twiddle = plan%twiddles(stage_params%stage)%twiddles_factor(j + 1)

        select case (plan%algorithm%decimation_method%id)
        case (DIT%id)
            result(i0) = even + twiddle * odd
            result(i1) = even - twiddle * odd
        case (DIF%id)
            result(i0) = even + odd
            result(i1) = (even - odd) * twiddle
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
    complex(sp) :: even, odd, twiddle
    integer(isp) :: j, block_index, base
    integer(isp) :: i0, i1

    result = signal

    !$omp parallel do default(none) private(block_index, base, j, even, odd, twiddle, i0, i1) &
    !$omp& shared(result, plan, stage_params) &
    !$omp& num_threads(threads)
    do block_index = 0, stage_params%nb_blocks - 1
        base = block_index * stage_params%current_block_size
        do j = 0, stage_params%block_size - 1
            i0 = base + j + 1
            i1 = i0 + stage_params%block_size
            even = result(i0)
            odd = result(i1)

            twiddle = plan%twiddles(stage_params%stage)%twiddles_factor(j + 1)
            select case (plan%algorithm%decimation_method%id)
            case (DIT%id)
                result(i0) = even + twiddle * odd
                result(i1) = even - twiddle * odd
            case (DIF%id)
                result(i0) = even + odd
                result(i1) = (even - odd) * twiddle
            end select
        end do
    end do
    !$omp end parallel do

end function compute_openmp_cmplx_sp

end submodule NAFPack_Fourier_Transform_fft_compute_radix2
