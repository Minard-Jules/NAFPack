submodule(NAFPack_Fourier_Transform:NAFPack_Fourier_Transform_fft) NAFPack_Fourier_Transform_fft_compute_split_radix

    implicit none(type, external)

contains

    recursive module function compute_fft_split_radix_recursive_cmplx_sp( &
        signal, plan, stage) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        integer(isp), intent(in) :: stage
        complex(sp), dimension(size(signal)) :: result
        complex(sp), dimension(size(signal)/2) :: even
        complex(sp), dimension(size(signal)/4) :: twiddle_wk, twiddle_w3k, odd1, odd2, even1, even2
        integer(isp) :: N, i0, i1, i2

        N = size(signal)
        i0 = N/4
        i1 = i0 + N/4
        i2 = i1 + N/4

        if (N == 1) then
            result = signal
            return
        else if (N == 2) then
            result(1) = signal(1) + signal(2)
            result(2) = signal(1) - signal(2)
            return
        end if

        if (iand(N, N-1) /= 0) error stop "N must be power of 2"

        twiddle_wk = plan%split_radix_twiddles(stage)%twiddles_wk
        twiddle_w3k = plan%split_radix_twiddles(stage)%twiddles_w3k
        select case (plan%algorithm%decimation_method%id)
        case (DIT%id)
            even = compute_fft_split_radix_recursive_cmplx_sp(signal(1:N:2), plan, stage-1)
            odd1 = compute_fft_split_radix_recursive_cmplx_sp(signal(2:N:4), plan, stage-2)
            odd2 = compute_fft_split_radix_recursive_cmplx_sp(signal(4:N:4), plan, stage-2)
            even1 = even(1:i0)
            even2 = even(i0+1:)

            ! radix-2 part
            result( 1:i0) = even1 + (odd1 * twiddle_wk + odd2 * twiddle_w3k)
            result(i1+1:i2) = even1 - (odd1 * twiddle_wk + odd2 * twiddle_w3k)

            ! radix-4 part
            result(i0+1:i1) = even2 - im_sp * (odd1 * twiddle_wk - odd2 * twiddle_w3k)
            result(i2+1:) = even2 + im_sp * (odd1 * twiddle_wk - odd2 * twiddle_w3k)
        case (DIF%id)
            ! radix-2 part
            even(1:i0) = signal(1:i0) + signal(i1+1:i2)
            even(i0+1:) = signal(i0+1:i1) + signal(i2+1:)

            ! radix-4 part
            odd1 = ((signal(1:i0) - signal(i1+1:i2)) - im_sp * (signal(i0+1:i1) - signal(i2+1:))) * twiddle_wk
            odd2 = ((signal(1:i0) - signal(i1+1:i2)) + im_sp * (signal(i0+1:i1) - signal(i2+1:))) * twiddle_w3k

            even = compute_fft_split_radix_recursive_cmplx_sp(even, plan, stage+1)
            odd1 = compute_fft_split_radix_recursive_cmplx_sp(odd1, plan, stage+2)
            odd2 = compute_fft_split_radix_recursive_cmplx_sp(odd2, plan, stage+2)

            result(1:N:2) = even
            result(2:N:4) = odd1
            result(4:N:4) = odd2
        end select

    end function compute_fft_split_radix_recursive_cmplx_sp

    module function compute_fft_split_radix_cmplx_sp( &
        signal, plan, stage_params, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params
        type(LoopMethod), intent(in) :: loop_method
        complex(sp), dimension(plan%N) :: result

        if (loop_method%use_do_classic) then
            result = compute_do_classic_cmplx_sp( &
                    signal, plan, stage_params)
            ! else if (loop_method%use_vectorized) then
            !     result = compute_vectorized_cmplx_sp( &
            !              signal, plan, stage_params)
            ! else if (loop_method%use_do_concurrent) then
            !     result = compute_do_concurrent_cmplx_sp( &
            !              signal, plan, stage_params)
            ! else if (loop_method%parallel%use_openmp) then
            !     result = compute_openmp_cmplx_sp( &
            !              signal, plan, stage_params, &
            !              loop_method%parallel%num_threads)
        end if

    end function compute_fft_split_radix_cmplx_sp

    pure function compute_do_classic_cmplx_sp( &
        signal, plan, stage_params) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params
        complex(sp), dimension(plan%N) :: result

        result = signal

        if (stage_params%block_size < 1) then
            call apply_radix2_stage(result, plan%N)
        else
            call apply_split_radix_stage(result, plan, stage_params)
        end if

    end function compute_do_classic_cmplx_sp

    pure subroutine apply_radix2_stage(result, N)
        complex(sp), dimension(:), intent(inout) :: result
        integer(isp), intent(in) :: N

        integer(isp) :: max_iter_radix2, j, start_index, stride

        max_iter_radix2 = nint(log(real(N)) / log(4.0))
        stride = 4
        start_index = 1

        do j = 1, max_iter_radix2
            call process_radix2_butterflies(result, start_index, stride, N)
            start_index = 2 * stride - 1
            stride = 4 * stride
        end do

    end subroutine apply_radix2_stage

    pure subroutine process_radix2_butterflies(result, start_index, stride, N)
        complex(sp), dimension(:), intent(inout) :: result
        integer(isp), intent(in) :: start_index, stride, N

        integer(isp) :: i0, i1
        complex(sp) :: even, odd

        do i0 = start_index, N, stride
            i1 = i0 + 1
            if (i1 > N) cycle

            even = result(i0)
            odd = result(i1)

            result(i0) = even + odd
            result(i1) = even - odd
        end do

    end subroutine process_radix2_butterflies

    pure subroutine apply_split_radix_stage(result, plan, stage_params)
        complex(sp), dimension(:), intent(inout) :: result
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params

        integer(isp) :: j
        complex(sp) :: twiddles_wk, twiddles_w3k

        do j = 1, stage_params%block_size
            ! load twiddle factors for this j
            twiddles_wk = plan%split_radix_twiddles(stage_params%stage)%twiddles_wk(j)
            twiddles_w3k = plan%split_radix_twiddles(stage_params%stage)%twiddles_w3k(j)

            ! process all internal stages for this j
            call process_split_radix_stages_for_j(result, plan, stage_params, j, &
                                                twiddles_wk, twiddles_w3k)
        end do

    end subroutine apply_split_radix_stage

    pure subroutine process_split_radix_stages_for_j(result, plan, stage_params, j, &
                                                    twiddles_wk, twiddles_w3k)
        complex(sp), dimension(:), intent(inout) :: result
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params
        integer(isp), intent(in) :: j
        complex(sp), intent(in) :: twiddles_wk, twiddles_w3k

        integer(isp) :: inner_stage, start_index, stride

        do inner_stage = 0, plan%split_radix_twiddles(stage_params%stage)%indices(j)%num_stages
            ! Get start index and stride for this inner stage
            start_index = plan%split_radix_twiddles(stage_params%stage)%indices(j)%start_indices(inner_stage)
            stride = plan%split_radix_twiddles(stage_params%stage)%indices(j)%strides(inner_stage)

            ! Process split-radix butterflies (L shaped) for this stage
            call process_L_butterflies(result, plan, stage_params, &
                                    start_index, stride, &
                                    twiddles_wk, twiddles_w3k)
        end do

    end subroutine process_split_radix_stages_for_j

    pure subroutine process_L_butterflies(result, plan, stage_params, &
                                        start_index, stride, &
                                        twiddles_wk, twiddles_w3k)
        complex(sp), dimension(:), intent(inout) :: result
        type(FFTPlan), intent(in) :: plan
        type(FFTStageParams), intent(in) :: stage_params
        integer(isp), intent(in) :: start_index, stride
        complex(sp), intent(in) :: twiddles_wk, twiddles_w3k
        integer(isp) :: i1, i2, i3
        complex(sp) :: even1, even2, odd1, odd2
        complex(sp) :: temp_radix2, temp_radix4

        integer(isp) :: i0

        do i0 = start_index, plan%N - 1, stride
            i1 = i0 + stage_params%block_size
            i2 = i1 + stage_params%block_size
            i3 = i2 + stage_params%block_size
            if (i3 > plan%N) return

            even1 = result(i0)
            even2 = result(i1)
            odd1 = result(i2)
            odd2 = result(i3)

            select case (plan%algorithm%decimation_method%id)
            case (DIT%id)
                ! radix-2
                temp_radix2 = odd1 * twiddles_wk + odd2 * twiddles_w3k
                result(i0) = even1 + temp_radix2
                result(i2) = even1 - temp_radix2

                ! radix-4
                temp_radix4 = odd1 * twiddles_wk - odd2 * twiddles_w3k
                result(i1) = even2 - im_sp * temp_radix4
                result(i3) = even2 + im_sp * temp_radix4
            case (DIF%id)
                ! radix-2 part
                result(i0) = even1 + odd1
                result(i1) = even2 + odd2
                ! radix-4 part
                result(i2) = ((even1 - odd1) - im_sp * (even2 - odd2)) * twiddles_wk
                result(i3) = ((even1 - odd1) + im_sp * (even2 - odd2)) * twiddles_w3k
            end select

        end do

    end subroutine process_L_butterflies
    
end submodule NAFPack_Fourier_Transform_fft_compute_split_radix
