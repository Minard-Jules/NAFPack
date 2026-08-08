submodule(NAFPack_Fourier_Transform) NAFPack_plan_execute

    use NAFPack_butterfly, only: execute_butterfly_radix2_recursive, execute_butterfly_radix2, &
                                 execute_butterfly_radix3_recursive, execute_butterfly_radix3, &
                                 execute_butterfly_radix4_recursive, execute_butterfly_radix4, &
                                 execute_butterfly_radix5_recursive, execute_butterfly_radix5, &
                                 execute_butterfly_dft, execute_butterfly_mixed_radix_recursive, &
                                 execute_butterfly_split_radix_recursive, execute_butterfly_split_radix

    use NAFPack_utils, only: present_arg

    implicit none(type, external)

contains

    module subroutine execute_fft_plan_sp(this, signal_re, signal_im, result_re, result_im, loop_method, implementation_type)
        class(Fourier_Transform), intent(in) :: this
        real(sp), dimension(:), intent(in) :: signal_re, signal_im
        real(sp), dimension(:), intent(out) :: result_re, result_im
        type(LoopMethod), optional, intent(in) :: loop_method
        type(ImplementationType), optional, intent(in) :: implementation_type
        type(ImplementationType) :: implementation_type_val

        if (.not. this%plan%is_initialized) &
            error stop "Error in execute_fft_plan_sp: FFT plan is not initialized"

        if (this%plan%algorithm%id == ALG_NONE%id) &
            error stop "Error in execute_fft_plan_sp: FFT plan algorithm is not set"

        result_re = signal_re
        result_im = signal_im

        implementation_type_val = present_arg(implementation_type, DEFAULT_IMPLEMENTATION_TYPE)

        select case (implementation_type_val%id)
        case (recursive%id)
            call execute_recursive_fft_plan_sp(this%plan, result_re, result_im)
        case (ITERATIVE%id)
            call execute_iterative_fft_plan_sp(this%plan, result_re, result_im, present_arg(loop_method, default_loop_method))
        case default
            error stop "Error in execute_fft_plan_sp: unknown implementation type"
        end select

    end subroutine execute_fft_plan_sp

    pure subroutine execute_recursive_fft_plan_sp(plan, signal_re, signal_im)
        type(FFTPlan), intent(in) :: plan
        real(sp), dimension(:), intent(inout) :: signal_re, signal_im
        integer(isp) :: num_stages

        select case (plan%algorithm%decimation_method%id)
        case (DIT%id)
            num_stages = plan%num_stages
        case (DIF%id)
            num_stages = 1
        case default
            error stop "Error in execute_recursive_fft_plan_sp: unknown decimation method"
        end select

        if (plan%use_radix2) then
            call execute_butterfly_radix2_recursive(plan%stages, plan%N, signal_re, signal_im, num_stages, plan%algorithm%decimation_method)
        else if (plan%use_radix3) then
            call execute_butterfly_radix3_recursive(plan%stages, plan%N, signal_re, signal_im, num_stages, plan%algorithm%decimation_method)
        else if (plan%use_radix4) then
            call execute_butterfly_radix4_recursive(plan%stages, plan%N, signal_re, signal_im, num_stages, plan%algorithm%decimation_method)
        else if (plan%use_radix5) then
            call execute_butterfly_radix5_recursive(plan%stages, plan%N, signal_re, signal_im, num_stages, plan%algorithm%decimation_method)
        else if (plan%use_mixed_radix) then
            call execute_butterfly_mixed_radix_recursive(plan%stages, plan%N, signal_re, signal_im, num_stages, plan%algorithm%decimation_method)
        else if (plan%use_split_radix) then
            call execute_butterfly_split_radix_recursive(plan%stages, plan%N, signal_re, signal_im, num_stages, plan%algorithm%decimation_method)
        else
            error stop "Error in execute_recursive_fft_plan_sp: unsupported radix"
        end if

    end subroutine execute_recursive_fft_plan_sp

    subroutine execute_iterative_fft_plan_sp(plan, signal_re, signal_im, loop_method)
        type(FFTPlan), intent(in) :: plan
        real(sp), dimension(:), intent(inout) :: signal_re, signal_im
        type(LoopMethod), intent(in) :: loop_method

        if (plan%algorithm%decimation_method%id == DIT%id) then
            if (plan%use_radix2 .or. plan%use_split_radix) then
                call bit_reversal(signal_re, signal_im, plan%N)
            else
                call digit_reversal(signal_re, signal_im, plan%N, plan%radix_plan)
            end if
            call execute_stages(plan, signal_re, signal_im, loop_method)
        else
            call execute_stages(plan, signal_re, signal_im, loop_method)
            if (plan%use_radix2 .or. plan%use_split_radix) then
                call bit_reversal(signal_re, signal_im, plan%N)
            else
                call digit_reversal(signal_re, signal_im, plan%N, plan%radix_plan)
            end if
        end if

    end subroutine execute_iterative_fft_plan_sp

    pure subroutine bit_reversal(re, im, N)
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re
        real(sp), dimension(N), intent(inout) :: im
        integer(isp) :: i, j, k
        real(sp) :: temp_re, temp_im

        j = 0
        do i = 0, N - 2
            if (i < j) then
                temp_re = re(i + 1)
                re(i + 1) = re(j + 1)
                re(j + 1) = temp_re

                temp_im = im(i + 1)
                im(i + 1) = im(j + 1)
                im(j + 1) = temp_im
            end if

            k = N / 2
            do while (k <= j)
                j = j - k
                k = k / 2
            end do
            j = j + k
        end do
    end subroutine bit_reversal

    pure subroutine digit_reversal(re, im, N, radix_plan)
        integer(isp), intent(in) :: N
        integer(isp), dimension(:), intent(in) :: radix_plan
        real(sp), dimension(N), intent(inout) :: re
        real(sp), dimension(N), intent(inout) :: im
        integer(isp) :: i, num_stages, curr, nxt
        real(sp) :: tmp_re, tmp_im
        logical :: done(N)

        num_stages = size(radix_plan)
        done = .false.

        do i = 0, N - 1
            if (done(i + 1)) cycle

            nxt = digit_rev_idx(i, num_stages, radix_plan)

            if (nxt /= i) then
                tmp_re = re(i + 1)
                tmp_im = im(i + 1)
                curr = i
                do while (nxt /= i)
                    re(curr + 1) = re(nxt + 1)
                    im(curr + 1) = im(nxt + 1)
                    done(curr + 1) = .true.
                    curr = nxt
                    nxt = digit_rev_idx(curr, num_stages, radix_plan)
                end do
                re(curr + 1) = tmp_re
                im(curr + 1) = tmp_im
                done(curr + 1) = .true.
            end if
            done(i + 1) = .true.
        end do

    end subroutine digit_reversal

    pure function digit_rev_idx(idx, num_stages, radix_plan) result(rev_idx)
        integer(isp), intent(in) :: idx
        integer(isp), intent(in) :: num_stages
        integer(isp), dimension(:), intent(in) :: radix_plan
        integer(isp) :: rev_idx, n, s, weight
        integer(isp) :: digits(num_stages)

        n = idx
        do s = 1, num_stages
            digits(s) = mod(n, radix_plan(s))
            n = n / radix_plan(s)
        end do

        rev_idx = 0
        weight = 1
        do s = num_stages, 1, -1
            rev_idx = rev_idx + digits(s) * weight
            if (s > 1) weight = weight * radix_plan(s)
        end do
    end function digit_rev_idx

    subroutine execute_stages(plan, data_re, data_im, loop_method)
        type(FFTPlan), intent(in) :: plan
        real(sp), dimension(:), intent(inout) :: data_re
        real(sp), dimension(:), intent(inout) :: data_im
        type(LoopMethod), intent(in) :: loop_method
        integer(isp) :: s

        do s = 1, plan%num_stages

            if (plan%use_split_radix) then
                call execute_butterfly_split_radix( &
                    plan%stages(s), plan%N, &
                    data_re, data_im, &
                    plan%algorithm%decimation_method)
                cycle
            end if

            select case (plan%stages(s)%radix)
            case (2)
                call execute_butterfly_radix2( &
                    plan%stages(s), plan%N, &
                    data_re, data_im, &
                    loop_method, plan%algorithm%decimation_method)
            case (3)
                call execute_butterfly_radix3( &
                    plan%stages(s), plan%N, &
                    data_re, data_im, &
                    plan%algorithm%decimation_method)
            case (4)
                call execute_butterfly_radix4( &
                    plan%stages(s), plan%N, &
                    data_re, data_im, &
                    plan%algorithm%decimation_method)
            case (5)
                call execute_butterfly_radix5( &
                    plan%stages(s), plan%N, &
                    data_re, data_im, &
                    plan%algorithm%decimation_method)
            case default
                ! print *, "Warning: stage ", s, " has unsupported radix ", plan%stages(s)%radix, ". Falling back to DFT."
                call execute_butterfly_dft( &
                    plan%stages(s), plan%N, &
                    data_re, data_im, &
                    plan%algorithm%decimation_method)
            end select
        end do
    end subroutine execute_stages

end submodule NAFPack_plan_execute
