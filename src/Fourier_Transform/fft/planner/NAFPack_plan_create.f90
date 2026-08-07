submodule(NAFPack_Fourier_Transform) NAFPack_plan_create

    use NAFPack_math_utils, only: sieve_of_eratosthenes, is_power_of_two, is_power_of_p, power_of_p_exponent
    use NAFPack_utils, only: present_arg

    implicit none(type, external)

contains

    pure module subroutine create_fft_plan_sp(this, N, algorithm, decimation_method)
        class(Fourier_Transform), intent(inout) :: this
        integer(isp), intent(in) :: N
        type(FFTAlgorithm), optional, intent(in) :: algorithm
        type(DecimationMethod), optional, intent(in) :: decimation_method

        this%plan%N = N
        this%plan%algorithm = present_arg(algorithm, ALG_AUTO)

        if (this%plan%algorithm%id == ALG_AUTO%id) call resolve_auto_algorithm(this%plan, decimation_method)

        call build_radix_plan(this%plan)

        call build_twiddles(this%plan)

        this%plan%is_initialized = .true.

    end subroutine create_fft_plan_sp

    pure subroutine resolve_auto_algorithm(plan, decimation_method)
        type(FFTPlan), intent(inout) :: plan
        type(DecimationMethod), optional, intent(in) :: decimation_method
        type(DecimationMethod) :: decimation_method_use

        decimation_method_use = present_arg(decimation_method, DIT)

        if (is_power_of_two(plan%N)) then
            plan%algorithm = merge(ALG_RADIX2_DIT, ALG_RADIX2_DIF, decimation_method_use%id == DIT%id)
        else if (is_power_of_p(plan%N, 3)) then
            plan%algorithm = merge(ALG_RADIX3_DIT, ALG_RADIX3_DIF, decimation_method_use%id == DIT%id)
        else if (is_power_of_p(plan%N, 4)) then
            plan%algorithm = merge(ALG_RADIX4_DIT, ALG_RADIX4_DIF, decimation_method_use%id == DIT%id)
        else if (is_power_of_p(plan%N, 5)) then
            plan%algorithm = merge(ALG_RADIX5_DIT, ALG_RADIX5_DIF, decimation_method_use%id == DIT%id)
        else
            plan%algorithm = merge(ALG_MIXED_DIT, ALG_MIXED_DIF, decimation_method_use%id == DIT%id)
        end if
    end subroutine resolve_auto_algorithm

    pure subroutine build_radix_plan(plan)
        type(FFTPlan), intent(inout) :: plan

        select case (plan%algorithm%id)
        case (ALG_RADIX2_DIT%id, ALG_RADIX2_DIF%id)
            call build_radix2_plan(plan)
        case (ALG_RADIX3_DIT%id, ALG_RADIX3_DIF%id)
            call build_radix3_plan(plan)
        case (ALG_RADIX4_DIT%id, ALG_RADIX4_DIF%id)
            call build_radix4_plan(plan)
        case (ALG_RADIX5_DIT%id, ALG_RADIX5_DIF%id)
            call build_radix5_plan(plan)
        case (ALG_MIXED_DIT%id, ALG_MIXED_DIF%id)
            call build_mixed_radix_plan(plan)
        case (ALG_SPLIT_DIT%id, ALG_SPLIT_DIF%id)
            call build_split_radix_plan(plan)
        case default
            error stop "Error in build_radix_plan: unknown fft algorithm"
        end select

        if (plan%algorithm%decimation_method%id == DIF%id) &
            plan%radix_plan = plan%radix_plan(plan%num_stages:1:-1)

    end subroutine build_radix_plan

    pure subroutine build_split_radix_plan(plan)
        type(FFTPlan), intent(inout) :: plan

        if (.not. is_power_of_two(plan%N)) &
            error stop "Error in build_split_radix_plan: N must be a power of two for split-radix algorithm"

        plan%num_stages = power_of_p_exponent(plan%N, 2)
        allocate (plan%radix_plan(plan%num_stages))
        plan%radix_plan = 2
        plan%use_split_radix = .true.

    end subroutine build_split_radix_plan

    pure subroutine build_mixed_radix_plan(plan)
        type(FFTPlan), intent(inout) :: plan
        integer(isp), dimension(:), allocatable :: primes
        integer(isp), dimension(:), allocatable :: radix_exponents
        integer(isp) :: N_primes, i, remaining_N

        primes = sieve_of_eratosthenes(plan%N)
        N_primes = size(primes)
        allocate (radix_exponents(N_primes))
        radix_exponents = 0

        remaining_N = plan%N
        do i = 1, N_primes
            radix_exponents(i) = power_of_p_exponent(remaining_N, primes(i))
            remaining_N = remaining_N / primes(i)**radix_exponents(i)
        end do

        if (remaining_N /= 1) &
            error stop "Error in build_mixed_radix_plan: N has a prime factor larger than the largest prime in the sieve"

        plan%num_stages = sum(radix_exponents)
        allocate (plan%radix_plan(plan%num_stages))
        plan%radix_plan = 0
        do i = 1, N_primes
            plan%radix_plan(sum(radix_exponents(1:i - 1)) + 1:sum(radix_exponents(1:i))) = primes(i)
        end do

        plan%use_mixed_radix = .true.

        deallocate (radix_exponents, primes)

    end subroutine build_mixed_radix_plan

    pure subroutine build_radix2_plan(plan)
        type(FFTPlan), intent(inout) :: plan

        if (.not. is_power_of_two(plan%N) .and. plan%use_radix2) then
            error stop "Error in build_radix2_plan: N must be a power of two for radix-2 algorithm"
        else if (.not. is_power_of_two(plan%N) .and. plan%use_split_radix) then
            error stop "Error in build_radix2_plan: N must be a power of two for split radix algorithm"
        end if

        plan%num_stages = power_of_p_exponent(plan%N, 2)
        allocate (plan%radix_plan(plan%num_stages))
        plan%radix_plan = 2
        plan%use_radix2 = .true.

    end subroutine build_radix2_plan

    pure subroutine build_radix3_plan(plan)
        type(FFTPlan), intent(inout) :: plan

        if (.not. is_power_of_p(plan%N, 3)) &
            error stop "Error in build_radix3_plan: N must be a power of three for radix-3 algorithm"

        plan%num_stages = power_of_p_exponent(plan%N, 3)
        allocate (plan%radix_plan(plan%num_stages))
        plan%radix_plan = 3
        plan%use_radix3 = .true.

    end subroutine build_radix3_plan

    pure subroutine build_radix4_plan(plan)
        type(FFTPlan), intent(inout) :: plan

        if (.not. is_power_of_p(plan%N, 4)) &
            error stop "Error in build_radix4_plan: N must be a power of four for radix-4 algorithm"

        plan%num_stages = power_of_p_exponent(plan%N, 4)
        allocate (plan%radix_plan(plan%num_stages))
        plan%radix_plan = 4
        plan%use_radix4 = .true.

    end subroutine build_radix4_plan

    pure subroutine build_radix5_plan(plan)
        type(FFTPlan), intent(inout) :: plan

        if (.not. is_power_of_p(plan%N, 5)) &
            error stop "Error in build_radix5_plan: N must be a power of five for radix-5 algorithm"

        plan%num_stages = power_of_p_exponent(plan%N, 5)
        allocate (plan%radix_plan(plan%num_stages))
        plan%radix_plan = 5
        plan%use_radix5 = .true.

    end subroutine build_radix5_plan

    pure subroutine build_split_radix_twiddles(plan)
        type(FFTPlan), intent(inout) :: plan
        integer(isp) :: s, j, n_groups, butterfly_size, stride, offset, n_inner, two_stride
        integer(isp) :: start_index, group_stride
        integer(isp) :: max_depth
        integer(isp), dimension(:), allocatable :: tmp_base, tmp_stride
        logical :: is_dit

        if (plan%N < 4) &
            error stop "Error in build_split_radix_twiddles: N must be at least 4 for split-radix"

        max_depth = plan%num_stages / 2 + 1
        allocate (tmp_base(0:max_depth), tmp_stride(0:max_depth))

        if (.not. associated(plan%twiddles_factors%real)) allocate (plan%twiddles_factors%real(plan%N - 2))
        if (.not. associated(plan%twiddles_factors%imag)) allocate (plan%twiddles_factors%imag(plan%N - 2))
        if (.not. allocated(plan%stages)) allocate (plan%stages(plan%num_stages))

        is_dit = plan%algorithm%decimation_method%id == DIT%id
        two_stride = merge(1, plan%N / 2, is_dit)

        offset = 1
        do s = 1, plan%num_stages
            stride = two_stride / 2
            butterfly_size = 2 * two_stride

            plan%stages(s)%radix = 2
            plan%stages(s)%offset = offset
            plan%stages(s)%stride = stride
            plan%stages(s)%butterfly_size = butterfly_size

            plan%stages(s)%twiddles_real => &
                plan%twiddles_factors%real(offset:offset + 2 * stride - 1)
            plan%stages(s)%twiddles_imag => &
                plan%twiddles_factors%imag(offset:offset + 2 * stride - 1)

            if (stride >= 1) then
                ! Wk = W^k
                call compute_twiddles_for_stage_sp( &
                    plan%stages(s)%twiddles_real(1:stride), &
                    plan%stages(s)%twiddles_imag(1:stride), &
                    butterfly_size, stride, plan%inverse, power=1)

                ! W3k = W^3k
                call compute_twiddles_for_stage_sp( &
                    plan%stages(s)%twiddles_real(stride + 1:2 * stride), &
                    plan%stages(s)%twiddles_imag(stride + 1:2 * stride), &
                    butterfly_size, stride, plan%inverse, power=3)
            end if

            n_inner = 0
            tmp_base(0) = 0
            tmp_stride(0) = 2 * butterfly_size
            do while (tmp_base(n_inner) + 1 < plan%N)
                n_inner = n_inner + 1
                tmp_base(n_inner) = 2 * tmp_stride(n_inner - 1) - butterfly_size
                tmp_stride(n_inner) = 4 * tmp_stride(n_inner - 1)
            end do

            allocate (plan%stages(s)%groups%start_index(0:n_inner))
            allocate (plan%stages(s)%groups%group_stride(0:n_inner))
            plan%stages(s)%groups%start_index = tmp_base(0:n_inner)
            plan%stages(s)%groups%group_stride = tmp_stride(0:n_inner)

            offset = offset + 2 * stride
            two_stride = merge(two_stride * 2, two_stride / 2, is_dit)
        end do

    end subroutine build_split_radix_twiddles

    pure subroutine build_twiddles(plan)
        type(FFTPlan), intent(inout) :: plan
        integer(isp) :: s, p, offset, stride, current_radix
        logical :: is_dit

        if (plan%use_split_radix) then
            call build_split_radix_twiddles(plan)
            return
        end if

        if (.not. associated(plan%twiddles_factors%real)) allocate (plan%twiddles_factors%real(plan%N - 1))
        if (.not. associated(plan%twiddles_factors%imag)) allocate (plan%twiddles_factors%imag(plan%N - 1))
        if (.not. allocated(plan%stages)) allocate (plan%stages(plan%num_stages))

        is_dit = plan%algorithm%decimation_method%id == DIT%id
        stride = merge(1, plan%N / plan%radix_plan(1), is_dit)

        offset = 1
        do s = 1, plan%num_stages
            current_radix = plan%radix_plan(s)

            plan%stages(s)%radix = current_radix
            plan%stages(s)%offset = offset
            plan%stages(s)%stride = stride
            plan%stages(s)%butterfly_size = current_radix * stride

            plan%stages(s)%twiddles_real => &
                plan%twiddles_factors%real(offset:offset + stride * (current_radix - 1) - 1)
            plan%stages(s)%twiddles_imag => &
                plan%twiddles_factors%imag(offset:offset + stride * (current_radix - 1) - 1)

            do p = 1, current_radix - 1
                call compute_twiddles_for_stage_sp( &
                    plan%stages(s)%twiddles_real((p - 1) * stride + 1:p * stride), &
                    plan%stages(s)%twiddles_imag((p - 1) * stride + 1:p * stride), &
                    plan%stages(s)%butterfly_size, stride, plan%inverse, power=p)
            end do

            offset = offset + stride * (current_radix - 1)

            if (s < plan%num_stages) stride = merge(stride * current_radix, stride / plan%radix_plan(s + 1), is_dit)
        end do

    end subroutine build_twiddles

    pure subroutine compute_twiddles_for_stage_sp(twiddles_real, twiddles_imag, butterfly_size, number_of_twiddles, inverse, power)
        real(sp), dimension(:), intent(inout) :: twiddles_real, twiddles_imag
        integer(isp), intent(in) :: butterfly_size, number_of_twiddles
        logical, intent(in) :: inverse
        integer(isp), optional, intent(in) :: power
        integer(isp), dimension(number_of_twiddles) :: k
        integer(isp) :: i
        real(sp) :: sign, angle

        sign = merge(1.0_sp, -1.0_sp, inverse)
        angle = sign * 2.0_sp * pi_sp * real(present_arg(power, 1), sp) / real(butterfly_size, sp)

        k = [(i, i=0, number_of_twiddles - 1)]

        twiddles_real = cos(angle * real(k, sp))
        twiddles_imag = sin(angle * real(k, sp))

    end subroutine compute_twiddles_for_stage_sp

end submodule NAFPack_plan_create
