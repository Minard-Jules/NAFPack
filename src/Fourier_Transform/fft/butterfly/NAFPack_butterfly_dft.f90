submodule(NAFPack_butterfly) NAFPack_butterfly_dft

    use NAFPack_constant, only: pi_sp

    implicit none(type, external)

contains

    pure module subroutine execute_butterfly_dft(stage, N, re, im, decimation_method)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method
        integer(isp) :: m, k, num_groups, step, radix
        real(sp), dimension(stage%radix, stage%radix) :: dft_re, dft_im

        radix = stage%radix
        num_groups = N / stage%butterfly_size
        step = radix * stage%stride

        do k = 1, radix
            do m = 1, radix
                dft_re(m, k) = cos(-2.0_sp * pi_sp * real((m - 1) * (k - 1), sp) / real(radix, sp))
                dft_im(m, k) = sin(-2.0_sp * pi_sp * real((m - 1) * (k - 1), sp) / real(radix, sp))
            end do
        end do

        select case (decimation_method%id)
        case (DIT%id)
            call kernel_butterfly_dft_dit(stage, N, re, im, num_groups, step, radix, dft_re, dft_im)
        case (DIF%id)
            call kernel_butterfly_dft_dif(stage, N, re, im, num_groups, step, radix, dft_re, dft_im)
        case default
            error stop "Error in execute_butterfly_dft: unknown decimation method"
        end select

    end subroutine execute_butterfly_dft

    pure subroutine kernel_butterfly_dft_dit(stage, N, re, im, num_groups, step, radix, dft_re, dft_im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp), intent(in) :: num_groups, step, radix
        real(sp), dimension(stage%radix, stage%radix), intent(in) :: dft_re, dft_im
        integer(isp) :: group_idx, twiddle_idx, j, m, base
        real(sp), dimension(stage%radix) :: x_re, x_im, out_re, out_im
        real(sp) :: wr, wi, ar, ai

        do group_idx = 0, num_groups - 1
            do twiddle_idx = 1, stage%stride
                base = group_idx * step + twiddle_idx

                x_re(1) = re(base)
                x_im(1) = im(base)

                do j = 2, radix
                    ar = re(base + (j - 1) * stage%stride)
                    ai = im(base + (j - 1) * stage%stride)
                    wr = stage%twiddles_real(twiddle_idx + (j - 2) * stage%stride)
                    wi = stage%twiddles_imag(twiddle_idx + (j - 2) * stage%stride)
                    x_re(j) = ar * wr - ai * wi
                    x_im(j) = ar * wi + ai * wr
                end do

                do m = 1, radix
                    out_re(m) = dot_product(dft_re(:, m), x_re) - dot_product(dft_im(:, m), x_im)
                    out_im(m) = dot_product(dft_re(:, m), x_im) + dot_product(dft_im(:, m), x_re)
                end do
                do m = 1, radix
                    re(base + (m - 1) * stage%stride) = out_re(m)
                    im(base + (m - 1) * stage%stride) = out_im(m)
                end do
            end do
        end do

    end subroutine kernel_butterfly_dft_dit

    pure subroutine kernel_butterfly_dft_dif(stage, N, re, im, num_groups, step, radix, dft_re, dft_im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp), intent(in) :: num_groups, step, radix
        real(sp), dimension(stage%radix, stage%radix), intent(in) :: dft_re, dft_im
        integer(isp) :: group_idx, twiddle_idx, j, m, base
        real(sp), dimension(stage%radix) :: x_re, x_im, out_re, out_im
        real(sp) :: wr, wi

        do group_idx = 0, num_groups - 1
            do twiddle_idx = 1, stage%stride
                base = group_idx * step + twiddle_idx

                do j = 1, radix
                    x_re(j) = re(base + (j - 1) * stage%stride)
                    x_im(j) = im(base + (j - 1) * stage%stride)
                end do

                do m = 1, radix
                    out_re(m) = dot_product(dft_re(:, m), x_re) - dot_product(dft_im(:, m), x_im)
                    out_im(m) = dot_product(dft_re(:, m), x_im) + dot_product(dft_im(:, m), x_re)
                end do

                re(base) = out_re(1)
                im(base) = out_im(1)

                do m = 2, radix
                    wr = stage%twiddles_real(twiddle_idx + (m - 2) * stage%stride)
                    wi = stage%twiddles_imag(twiddle_idx + (m - 2) * stage%stride)
                    re(base + (m - 1) * stage%stride) = out_re(m) * wr - out_im(m) * wi
                    im(base + (m - 1) * stage%stride) = out_re(m) * wi + out_im(m) * wr
                end do
            end do
        end do
    end subroutine kernel_butterfly_dft_dif

end submodule NAFPack_butterfly_dft
