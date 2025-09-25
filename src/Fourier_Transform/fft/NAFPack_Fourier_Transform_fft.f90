submodule(NAFPack_Fourier_Transform) NAFPack_Fourier_Transform_fft

    use NAFPack_memory_management, only: realloc

    implicit none(type, external)

    type FFTStageParams
        integer(isp) :: stage
        integer(isp) :: current_block_size
        integer(isp) :: block_size
        integer(isp) :: nb_blocks
    end type FFTStageParams

    interface
        recursive module function compute_fft_radix2_recursive_cmplx_sp( &
            signal, plan, stage) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            type(FFTPlan), intent(in) :: plan
            integer(isp), intent(in) :: stage
            complex(sp), dimension(plan%N) :: result
        end function compute_fft_radix2_recursive_cmplx_sp

        module function compute_fft_radix2_iterative_cmplx_sp( &
            signal, plan, stage_params, loop_method) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            type(FFTPlan), intent(in) :: plan
            type(FFTStageParams), intent(in) :: stage_params
            type(LoopMethod), intent(in) :: loop_method
            complex(sp), dimension(plan%N) :: result
        end function compute_fft_radix2_iterative_cmplx_sp
    end interface

    interface
        recursive module function compute_fft_mixed_radix_recursive_cmplx_sp( &
            signal, plan, stage, radix) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            type(FFTPlan), intent(in) :: plan
            integer(isp), intent(in) :: stage, radix
            complex(sp), dimension(size(signal)) :: result
        end function compute_fft_mixed_radix_recursive_cmplx_sp

        module function compute_fft_mixed_radix_cmplx_sp( &
            signal, plan, stage_params, loop_method) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            type(FFTPlan), intent(in) :: plan
            type(FFTStageParams), intent(in) :: stage_params
            type(LoopMethod), intent(in) :: loop_method
            complex(sp), dimension(plan%N) :: result
        end function compute_fft_mixed_radix_cmplx_sp
    end interface

    interface
        recursive module function compute_fft_split_radix_recursive_cmplx_sp( &
            signal, plan, stage) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            type(FFTPlan), intent(in) :: plan
            integer(isp), intent(in) :: stage
            complex(sp), dimension(size(signal)) :: result
        end function compute_fft_split_radix_recursive_cmplx_sp

        module function compute_fft_split_radix_cmplx_sp( &
            signal, plan, stage_params, loop_method) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            type(FFTPlan), intent(in) :: plan
            type(FFTStageParams), intent(in) :: stage_params
            type(LoopMethod), intent(in) :: loop_method
            complex(sp), dimension(plan%N) :: result
        end function compute_fft_split_radix_cmplx_sp
    end interface

contains

    ! module function fft_plan_create(N) result(plan)
    !     integer(isp), intent(in) :: N
    !     complex(sp), dimension(:), allocatable :: plan
    !     integer(isp), dimension(:), allocatable :: radix_fft

    !     call get_radix(N, radix_fft)

    ! end function fft_plan_create

    module subroutine init_fft_plan_sp(this, N, algorithm, decimation_method)
        class(Fourier_Transform), intent(inout) :: this
        integer(isp), intent(in) :: N
        type(FFTAlgorithm), optional, intent(in) :: algorithm
        type(DecimationMethod), optional, intent(in) :: decimation_method
        integer(isp) :: N_radix

        this%fft_plan%N = N
        this%fft_plan%is_initialized = .true.

        if (present(algorithm)) then
            this%fft_plan%algorithm = algorithm
        else
            this%fft_plan%algorithm = ALG_AUTO
        end if

        select case (this%fft_plan%algorithm%id)
        case (ALG_RADIX2_DIT%id, ALG_RADIX2_DIF%id)
            if (.not. is_power_of_two(N)) &
                error stop "Error in init_fft_plan_sp: N must be a power of two for radix-2 DIT FFT."
            call get_radix_2(this%fft_plan)
            this%fft_plan%use_pure_radix2 = .true.
        case (ALG_SPLIT_DIF%id, ALG_SPLIT_DIT%id)
            if (.not. is_power_of_two(N)) &
                error stop "Error in init_fft_plan_sp: N must be a power of two for split FFT."
            call get_radix_2(this%fft_plan)
            this%fft_plan%use_split_radix = .true.
        case (ALG_MIXED_DIT%id, ALG_MIXED_DIF%id)
            call get_mixed_radix(N, this%fft_plan)
            this%fft_plan%use_mixed_radix = .true.
        case (ALG_AUTO%id)

            if (present(decimation_method)) then
                this%fft_plan%algorithm%decimation_method = decimation_method
            else
                this%fft_plan%algorithm%decimation_method = DIT
            end if

            if (is_power_of_two(N)) then
                call get_radix_2(this%fft_plan)
                this%fft_plan%use_pure_radix2 = .true.
            else
                call get_mixed_radix(N, this%fft_plan)
                this%fft_plan%use_mixed_radix = .true.
            end if

        case (ALG_NONE%id)
            print*,"Error in init_fft_plan_sp: FFT algorithm must be specified."
            error stop
        case default
            print*,"Error in init_fft_plan_sp: Unknown FFT algorithm."
            error stop
        end select

        N_radix = size(this%fft_plan%radix_plan)

        if (this%fft_plan%algorithm%decimation_method%id == DIF%id) then
            this%fft_plan%radix_plan = this%fft_plan%radix_plan(N_radix:1:-1)
        end if

        if (this%fft_plan%use_split_radix) then
            call generate_split_twiddles_sp(this%fft_plan, N_radix)
        else
            call generate_twiddles_sp(this%fft_plan, N_radix)
        end if

        if (this%fft_plan%algorithm%decimation_method%id == DIF%id) then
            if (this%fft_plan%use_split_radix) then
                this%fft_plan%split_radix_twiddles = this%fft_plan%split_radix_twiddles(N_radix:1:-1)
            else
                this%fft_plan%twiddles = this%fft_plan%twiddles(N_radix:1:-1)
            end if
        end if

    end subroutine init_fft_plan_sp

    pure subroutine get_radix_2(plan)
        type(FFTPlan), intent(inout) :: plan
        integer(isp) :: N_radix

        N_radix = power_of_p_exponent(plan%N, 2)
        allocate (plan%radix_plan(N_radix))
        plan%radix_plan = 2

    end subroutine get_radix_2

    pure subroutine get_mixed_radix(N, plan, group_radices)
        integer(isp), intent(in) :: N
        logical, optional, intent(in) :: group_radices
        type(FFTPlan), intent(inout) :: plan
        integer(isp), dimension(:), allocatable :: radix_exponents
        integer(isp), dimension(:), allocatable :: prime_number
        integer(isp) :: N_prime, i, N_temp, N_radix, idx, idx_old
        logical :: group_radices_used

        group_radices_used = .false.
        if (present(group_radices)) group_radices_used = group_radices

        N_temp = N
        prime_number = sieve_of_eratosthenes(N)
        N_prime = size(prime_number, 1)
        allocate (radix_exponents(N_prime))
        radix_exponents = 0

        do i = 1, size(prime_number)
            do while (mod(N_temp, prime_number(i)) == 0)
                radix_exponents(i) = radix_exponents(i) + 1
                N_temp = N_temp / prime_number(i)
            end do
        end do

        if (group_radices_used) then
            N_radix = count(radix_exponents /= 0)
        else
            N_radix = sum(radix_exponents)
        end if
        if (N_temp > 1) N_radix = N_radix + 1

        allocate (plan%radix_plan(N_radix))
        idx = 0
        idx_old = 0

        do i = 1, N_prime
            if (radix_exponents(i) /= 0) then
                if (group_radices_used) then
                    idx = idx + 1
                    plan%radix_plan(idx) = prime_number(i)**radix_exponents(i)
                else
                    idx = idx_old + radix_exponents(i)
                    plan%radix_plan((idx_old) + 1:idx) = prime_number(i)
                    idx_old = idx
                end if
            end if
        end do

        if (N_temp > 1) plan%radix_plan(N_radix) = N_temp

        deallocate (radix_exponents, prime_number)
    end subroutine get_mixed_radix

    pure subroutine generate_twiddles_sp(plan, N_radix)
        type(FFTPlan), intent(inout) :: plan
        integer(isp), intent(in) :: N_radix
        integer(isp), dimension(:), allocatable :: k
        integer(isp) :: i, j, r, block_size, current_block_size
        real(sp) :: sign_factor

        if (plan%algorithm%decimation_method%id == DIF%id) then
            sign_factor = -1.0_sp ! Pour DIF
        else
            sign_factor = -1.0_sp ! Pour DIT (défaut)
        end if

        allocate (plan%twiddles(N_radix))

        block_size = 1
        do i = 1, N_radix
            r = plan%radix_plan(i)
            current_block_size = block_size * r

            plan%twiddles(i)%radix = r
            plan%twiddles(i)%block_size = block_size
            plan%twiddles(i)%current_block_size = current_block_size

            allocate (plan%twiddles(i)%twiddles_factor(block_size))
            allocate (k(0:block_size - 1))
            k = [(j, j=0, block_size - 1)]
            plan%twiddles(i)%twiddles_factor = exp(sign_factor * 2.0_sp * im_sp * pi_sp * k / current_block_size)
            deallocate (k)

            block_size = current_block_size
        end do
    end subroutine generate_twiddles_sp

    subroutine generate_split_twiddles_sp(plan, N_radix)
        type(FFTPlan), intent(inout) :: plan
        integer(isp), intent(in) :: N_radix
        integer(isp), parameter :: MAX_ITER = 20
        integer(isp), dimension(:), allocatable :: tmp_start_indices, tmp_strides
        integer(isp), dimension(:), allocatable :: k
        integer(isp) :: i, j, r, block_size, current_block_size
        real(sp) :: sign_factor
        integer(isp) :: start_index, stride, num_stages, tmp_size

        if (plan%algorithm%decimation_method%id == DIF%id) then
            sign_factor = 1.0_sp
        else
            sign_factor = -1.0_sp
        end if

        allocate (plan%split_radix_twiddles(N_radix))

        current_block_size = 2
        do i = 1, N_radix
            r = plan%radix_plan(i)
            block_size = current_block_size / (2 * r)

            ! Initialize stage parameters
            plan%split_radix_twiddles(i)%radix = r
            plan%split_radix_twiddles(i)%block_size = block_size
            plan%split_radix_twiddles(i)%current_block_size = current_block_size

            ! Allocation and computation of twiddles_Wk and twiddles_W3k
            if (block_size < 1) then
                allocate (plan%split_radix_twiddles(i)%twiddles_Wk(1))
                allocate (plan%split_radix_twiddles(i)%twiddles_W3k(1))
                plan%split_radix_twiddles(i)%twiddles_Wk = (1._sp, 0._sp)
                plan%split_radix_twiddles(i)%twiddles_W3k = (1._sp, 0._sp)
            else
                allocate (plan%split_radix_twiddles(i)%twiddles_Wk(block_size))
                allocate (plan%split_radix_twiddles(i)%twiddles_W3k(block_size))

                allocate (k(0:block_size - 1))
                k = [(j, j=0, block_size - 1)]
                plan%split_radix_twiddles(i)%twiddles_Wk = exp(-2.0_sp * im_sp * pi_sp * k / current_block_size)
                plan%split_radix_twiddles(i)%twiddles_W3k = exp(-2.0_sp * im_sp * pi_sp * 3 * k / current_block_size)
                deallocate (k)
            end if

            ! Allocation and computation of stages split radix
            allocate (plan%split_radix_twiddles(i)%indices(block_size))
            do j = 1, block_size
                allocate (tmp_start_indices(0:MAX_ITER), tmp_strides(0:MAX_ITER))

                num_stages = 0
                start_index = j
                stride = 2 * current_block_size

                tmp_start_indices(num_stages) = start_index
                tmp_strides(num_stages) = stride

                do while (start_index < plan%N)
                    num_stages = num_stages + 1
                    tmp_size = size(tmp_start_indices)
                    if (num_stages > size(tmp_start_indices)) then
                        call realloc(tmp_start_indices, 2 * tmp_size)
                        call realloc(tmp_strides, 2 * tmp_size)
                    end if
                    start_index = 2 * stride - current_block_size + j
                    stride = 4 * stride
                    tmp_start_indices(num_stages) = start_index
                    tmp_strides(num_stages) = stride
                end do
                plan%split_radix_twiddles(i)%indices(j)%num_stages = num_stages

                allocate (plan%split_radix_twiddles(i)%indices(j)%start_indices(0:num_stages))
                allocate (plan%split_radix_twiddles(i)%indices(j)%strides(0:num_stages))
                plan%split_radix_twiddles(i)%indices(j)%start_indices = tmp_start_indices(0:num_stages)
                plan%split_radix_twiddles(i)%indices(j)%strides = tmp_strides(0:num_stages)
                deallocate (tmp_start_indices, tmp_strides)
            end do

            current_block_size = current_block_size * r
        end do
    end subroutine generate_split_twiddles_sp

    module function fft_cmplx_sp(this, signal, loop_method, implementation_type) result(result)
        class(Fourier_Transform), intent(inout) :: this
        complex(sp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        type(ImplementationType), optional, intent(in) :: implementation_type
        complex(sp), dimension(:), allocatable :: result
        type(ImplementationType) :: implementation_type_used
        type(LoopMethod) :: loop_method_used

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        if (present(implementation_type)) then
            implementation_type_used = implementation_type
        else
            implementation_type_used = ITERATIVE
        end if

        if (.not. this%fft_plan%is_initialized) then
            call init_fft_plan_sp(this, size(signal))
        else if (this%fft_plan%N /= size(signal)) then
            print*,"Error in fft_cmplx_sp: FFT plan size does not match signal size."
            error stop
        end if

        select case (implementation_type_used%id)
        case (ITERATIVE%id)
            result = fft_sp_iterative(signal, this%fft_plan, loop_method_used)
        case (recursive%id)
            result = fft_sp_recursive(signal, this%fft_plan)
        end select

    end function fft_cmplx_sp

    function fft_sp_iterative(signal, plan, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(LoopMethod), intent(in) :: loop_method
        complex(sp), dimension(:), allocatable :: result
        complex(sp), dimension(:), allocatable :: signal_reversed
        integer(isp) :: N_radix

        N_radix = size(plan%radix_plan)

        select case (plan%algorithm%decimation_method%id)
        case (DIT%id)
            if (plan%use_pure_radix2 .or. plan%use_split_radix) then
                signal_reversed = bit_reverse(signal, plan%N)
            else if (plan%use_mixed_radix) then
                signal_reversed = digit_reverse(signal, plan%radix_plan)
            else
                print*,"Error in fft_cmplx_sp: Unknown FFT plan type."
                error stop
            end if
            result = compute_fft_sp_iterative(signal_reversed, plan, loop_method)
        case (DIF%id)
            result = compute_fft_sp_iterative(signal, plan, loop_method)
            if (plan%use_pure_radix2 .or. plan%use_split_radix) then
                result = bit_reverse(result, plan%N)
            else if (plan%use_mixed_radix) then
                result = digit_reverse(result, plan%radix_plan(N_radix:1:-1))
            else
                print*,"Error in fft_cmplx_sp: Unknown FFT plan type."
                error stop
            end if
        end select

    end function fft_sp_iterative

    pure function bit_reverse(x, N) result(y)
        complex(sp), dimension(N), intent(in) :: x
        integer(isp), intent(in) :: N
        complex(sp), dimension(N) :: y
        integer(isp) :: i, j, m
        complex(sp) :: tmp

        y = x

        j = 0
        do i = 0, N - 1
            if (i < j) then
                tmp = y(i + 1)
                y(i + 1) = y(j + 1)
                y(j + 1) = tmp
            end if

            m = N / 2
            do while (m > 0 .and. j >= m)
                j = j - m
                m = m / 2
            end do
            j = j + m
        end do

    end function bit_reverse

    function digit_reverse(x, radix_plan) result(y)
        complex(sp), dimension(:), intent(in) :: x
        integer, dimension(:), intent(in) :: radix_plan ! [r1, r2, ..., rm]
        complex(sp), dimension(size(x)) :: y
        integer(sp), dimension(size(radix_plan)) :: d
        integer :: N, m, i, j, idx, rev, factor

        N = size(x)
        m = size(radix_plan)
        do i = 0, N - 1

            idx = i
            do j = m, 1, -1
                d(j) = mod(idx, radix_plan(j))
                idx = idx / radix_plan(j)
            end do

            rev = 0
            factor = 1
            do j = 1, m
                rev = rev + d(j) * factor
                factor = factor * radix_plan(j)
            end do

            y(rev + 1) = x(i + 1)
        end do

    end function digit_reverse

    function compute_fft_sp_iterative(signal, plan, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(LoopMethod), intent(in) :: loop_method
        complex(sp), dimension(:), allocatable :: result
        integer(isp) :: N, num_stages, stage, r
        type(FFTStageParams) :: stage_params

        N = size(signal)
        num_stages = size(plan%radix_plan)
        allocate (result(N))
        result = signal

        do stage = 1, num_stages

            stage_params%stage = stage
            if (allocated(plan%twiddles)) then
                r = plan%twiddles(stage)%radix
                stage_params%block_size = plan%twiddles(stage)%block_size
                stage_params%current_block_size = plan%twiddles(stage)%current_block_size
                stage_params%nb_blocks = N / stage_params%current_block_size
            else if (allocated(plan%split_radix_twiddles)) then
                r = plan%split_radix_twiddles(stage)%radix
                stage_params%block_size = plan%split_radix_twiddles(stage)%block_size
                stage_params%current_block_size = plan%split_radix_twiddles(stage)%current_block_size
                stage_params%nb_blocks = -1
            else
                print*,"Error in compute_fft_sp: Twiddles not allocated."
                error stop
            end if

            if (plan%use_pure_radix2) then
                result = compute_fft_radix2_iterative_cmplx_sp( &
                         result, plan, stage_params, loop_method)
            else if (plan%use_mixed_radix) then
                result = compute_fft_mixed_radix_cmplx_sp( &
                         result, plan, stage_params, loop_method)
            else if (plan%use_split_radix) then
                result = compute_fft_split_radix_cmplx_sp( &
                         result, plan, stage_params, loop_method)
            else
                print*,"Error in compute_fft_sp: Unknown FFT plan type."
                error stop
            end if
        end do

    end function compute_fft_sp_iterative

    function fft_sp_recursive(signal, plan) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        complex(sp), dimension(:), allocatable :: result
        integer(isp) :: N, num_stages, radix

        N = size(signal)

        select case (plan%algorithm%decimation_method%id)
        case (DIT%id)
            num_stages = size(plan%radix_plan)
        case (DIF%id)
            num_stages = 1
        end select
        radix = plan%radix_plan(size(plan%radix_plan))

        if (plan%use_pure_radix2) then
            result = compute_fft_radix2_recursive_cmplx_sp( &
                     signal, plan, num_stages)
        else if (plan%use_mixed_radix) then
            result = compute_fft_mixed_radix_recursive_cmplx_sp( &
                     signal, plan, num_stages, radix)
        else if (plan%use_split_radix) then
            result = compute_fft_split_radix_recursive_cmplx_sp( &
                     signal, plan, num_stages)
        else
            print*,"Error in compute_fft_sp: Unknown FFT plan type."
            error stop
        end if

    end function fft_sp_recursive

    pure module subroutine destroy_fft_plan_sp(this)
        class(Fourier_Transform), intent(inout) :: this
        integer(isp) :: i

        if (this%fft_plan%is_initialized) then

            if (allocated(this%fft_plan%twiddles)) then
                do i = 1, size(this%fft_plan%twiddles)
                    if (allocated(this%fft_plan%twiddles(i)%twiddles_factor)) then
                        deallocate (this%fft_plan%twiddles(i)%twiddles_factor)
                    end if
                end do
                deallocate (this%fft_plan%twiddles)
            end if

            if (allocated(this%fft_plan%split_radix_twiddles)) then
                do i = 1, size(this%fft_plan%split_radix_twiddles)
                    if (allocated(this%fft_plan%split_radix_twiddles(i)%twiddles_Wk)) then
                        deallocate (this%fft_plan%split_radix_twiddles(i)%twiddles_Wk)
                    end if
                    if (allocated(this%fft_plan%split_radix_twiddles(i)%twiddles_W3k)) then
                        deallocate (this%fft_plan%split_radix_twiddles(i)%twiddles_W3k)
                    end if
                end do
                deallocate (this%fft_plan%split_radix_twiddles)
            end if

            if (allocated(this%fft_plan%radix_plan)) then
                deallocate (this%fft_plan%radix_plan)
            end if

            this%fft_plan%N = 0
            this%fft_plan%is_initialized = .false.
            this%fft_plan%use_pure_radix2 = .false.
            this%fft_plan%use_split_radix = .false.
            this%fft_plan%use_mixed_radix = .false.
            this%fft_plan%algorithm = ALG_NONE
        end if
    end subroutine destroy_fft_plan_sp

end submodule NAFPack_Fourier_Transform_fft
