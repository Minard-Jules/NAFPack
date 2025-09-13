submodule(NAFPack_Fourier_Transform) NAFPack_Fourier_Transform_fft

    implicit none(type, external)

    interface
        module function compute_fft_cmplx_sp( &
            signal, plan, stage, current_block_size, block_size, nb_blocks, loop_method) result(result)
            complex(sp), dimension(:), intent(in) :: signal
            type(FFTPlan), intent(in) :: plan
            integer(isp), intent(in) :: stage, current_block_size, block_size, nb_blocks
            type(LoopMethod), intent(in) :: loop_method
            complex(sp), dimension(plan%N) :: result
        end function compute_fft_cmplx_sp
    end interface

contains

    ! module function fft_plan_create(N) result(plan)
    !     integer(isp), intent(in) :: N
    !     complex(sp), dimension(:), allocatable :: plan
    !     integer(isp), dimension(:), allocatable :: radix_fft

    !     call get_radix(N, radix_fft)

    ! end function fft_plan_create

    module subroutine init_fft_plan_sp(this, N)
        class(Fourier_Transform), intent(inout) :: this
        integer(isp), intent(in) :: N
        integer(isp) :: N_radix

        this%fft_plan%N = N
        this%fft_plan%is_initialized = .true.

        if (is_power_of_two(N)) then
            call get_radix_2_sp(this%fft_plan)
        else
            call get_radix(N, this%fft_plan%radix_plan)
        end if

        N_radix = size(this%fft_plan%radix_plan)

        allocate (this%fft_plan%twiddles(0:N_radix))
        this%fft_plan%twiddles(0)%block_size = this%fft_plan%N

        call generate_twiddles_sp(this%fft_plan, N_radix)

    end subroutine init_fft_plan_sp

    subroutine get_radix_2_sp(plan)
        type(FFTPlan), intent(inout) :: plan
        integer(isp) :: N_radix

        plan%is_pure_radix2 = .true.
        N_radix = power_of_p_exponent(plan%N, 2)
        allocate (plan%radix_plan(N_radix))
        plan%radix_plan = 2

    end subroutine get_radix_2_sp

    subroutine get_radix(N, radix_plan, group_radices)
        integer(isp), intent(in) :: N
        logical, optional, intent(in) :: group_radices
        integer(isp), dimension(:), allocatable, intent(out) :: radix_plan
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

        allocate (radix_plan(N_radix))
        idx = 0
        idx_old = 0

        do i = 1, N_prime
            if (radix_exponents(i) /= 0) then
                if (group_radices_used) then
                    idx = idx + 1
                    radix_plan(idx) = prime_number(i)**radix_exponents(i)
                else
                    idx = idx_old + radix_exponents(i)
                    radix_plan((idx_old) + 1:idx) = prime_number(i)
                    idx_old = idx
                end if
            end if
        end do

        if (N_temp > 1) radix_plan(N_radix) = N_temp

        deallocate (radix_exponents, prime_number)
    end subroutine get_radix

    subroutine generate_twiddles_sp(plan, N_radix)
        type(FFTPlan), intent(inout) :: plan
        integer(isp), intent(in) :: N_radix
        integer(isp), dimension(:), allocatable :: k
        integer(isp) :: i, j, r, block_size, current_block_size

        block_size = 1
        do i = 1, N_radix
            r = plan%radix_plan(i)
            current_block_size = block_size * r

            plan%twiddles(i)%radix = r
            plan%twiddles(i)%block_size = block_size

            allocate (plan%twiddles(i)%twiddles_factor(block_size))

            allocate (k(0:block_size - 1))
            k = [(j, j=0, block_size - 1)]
            plan%twiddles(i)%twiddles_factor = exp(-2.0_sp * im_sp * pi_sp * k / current_block_size)
            deallocate (k)

            block_size = current_block_size
        end do
    end subroutine generate_twiddles_sp

    module function fft_cmplx_sp(this, signal, loop_method) result(result)
        class(Fourier_Transform), intent(inout) :: this
        complex(sp), dimension(:), intent(in) :: signal
        type(LoopMethod), optional, intent(in) :: loop_method
        complex(sp), dimension(:), allocatable :: result
        complex(sp), dimension(:), allocatable :: signal_reversed
        type(LoopMethod) :: loop_method_used

        if (present(loop_method)) then
            loop_method_used = check_loop_method(loop_method)
        else
            loop_method_used = default_loop_method
        end if

        if (.not. this%fft_plan%is_initialized) then
            call init_fft_plan_sp(this, size(signal))
        else if (this%fft_plan%N /= size(signal)) then
            print*,"Error in fft_cmplx_sp: FFT plan size does not match signal size."
            error stop
        end if

        signal_reversed = bit_reverse(signal, this%fft_plan%N)

        result = compute_fft_radix2_sp(signal_reversed, this%fft_plan, loop_method_used)

    end function fft_cmplx_sp

    function bit_reverse(x, N) result(y)
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

    function compute_fft_radix2_sp(signal, plan, loop_method) result(result)
        complex(sp), dimension(:), intent(in) :: signal
        type(FFTPlan), intent(in) :: plan
        type(LoopMethod), intent(in) :: loop_method
        complex(sp), dimension(:), allocatable :: result
        integer(isp) :: N, num_stages, stage, block_size, r, nb_blocks, current_block_size

        N = size(signal)
        num_stages = size(plan%radix_plan)
        allocate (result(N))
        result = signal

        do stage = 1, num_stages
            r = plan%twiddles(stage)%radix
            block_size = plan%twiddles(stage)%block_size
            current_block_size = block_size * r
            nb_blocks = N / current_block_size

            result = compute_fft_cmplx_sp(&
            result, plan, stage, current_block_size, block_size, nb_blocks, loop_method)
        end do

    end function compute_fft_radix2_sp

    module subroutine destroy_fft_plan_sp(this)
        class(Fourier_Transform), intent(inout) :: this
        integer(isp) :: i

        if (this%fft_plan%is_initialized) then
            do i = 0, size(this%fft_plan%twiddles) - 1
                if (allocated(this%fft_plan%twiddles(i)%twiddles_factor)) then
                    deallocate (this%fft_plan%twiddles(i)%twiddles_factor)
                end if
            end do

            if (allocated(this%fft_plan%twiddles)) then
                deallocate (this%fft_plan%twiddles)
            end if

            if (allocated(this%fft_plan%radix_plan)) then
                deallocate (this%fft_plan%radix_plan)
            end if

            this%fft_plan%is_initialized = .false.
            this%fft_plan%is_pure_radix2 = .false.
        end if
    end subroutine destroy_fft_plan_sp

end submodule NAFPack_Fourier_Transform_fft
