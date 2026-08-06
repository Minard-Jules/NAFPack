submodule(NAFPack_Fourier_Transform) NAFPack_plan_destroy

    implicit none(type, external)

contains

    pure module subroutine destroy_fft_plan_sp(this)
        class(Fourier_Transform), intent(inout) :: this

        if (.not. this%plan%is_initialized) return

        call deallocate_fft_plan(this%plan)

        call zero_out_plan(this%plan)

    end subroutine destroy_fft_plan_sp

    pure subroutine deallocate_fft_plan(plan)
        type(FFTPlan), intent(inout) :: plan
        integer(isp) :: i

        if (allocated(plan%radix_plan)) deallocate (plan%radix_plan)

        if (allocated(plan%stages)) then
            do i = 1, plan%num_stages
                plan%stages(i)%twiddles_real => null()
                plan%stages(i)%twiddles_imag => null()
            end do
            deallocate (plan%stages)
        end if

        if (associated(plan%twiddles_factors%real)) then
            deallocate (plan%twiddles_factors%real)
            plan%twiddles_factors%real => null()
        end if

        if (associated(plan%twiddles_factors%imag)) then
            deallocate (plan%twiddles_factors%imag)
            plan%twiddles_factors%imag => null()
        end if

    end subroutine deallocate_fft_plan

    pure subroutine zero_out_plan(plan)
        type(FFTPlan), intent(inout) :: plan

        plan%N = 0
        plan%num_stages = 0
        plan%algorithm = ALG_NONE
        plan%inverse = .false.
        plan%use_radix2 = .false.
        plan%use_radix3 = .false.
        plan%use_radix4 = .false.
        plan%use_radix5 = .false.
        plan%use_split_radix = .false.
        plan%use_mixed_radix = .false.
        plan%is_initialized = .false.

    end subroutine zero_out_plan

end submodule NAFPack_plan_destroy
