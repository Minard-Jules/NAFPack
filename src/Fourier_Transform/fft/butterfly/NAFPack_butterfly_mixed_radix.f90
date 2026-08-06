submodule(NAFPack_butterfly) NAFPack_butterfly_mixed_radix

    use NAFPack_constant, only: pi_sp

    implicit none(type, external)

contains

    pure module subroutine execute_butterfly_mixed_radix_recursive(stages, N, re, im, stage_idx, decimation_method)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        if (stage_idx < 1 .or. stage_idx > size(stages)) &
            error stop "stage_idx out of bounds"

        select case (decimation_method%id)
            case (DIT%id)
                call execute_butterfly_mixed_radix_recursive_dit(stages, N, re, im, stage_idx)
            case (DIF%id)
                call execute_butterfly_mixed_radix_recursive_dif(stages, N, re, im, stage_idx)
            case default
                error stop "Error in execute_butterfly_mixed_radix_recursive: unknown decimation method"
        end select

    end subroutine execute_butterfly_mixed_radix_recursive

    pure recursive subroutine execute_butterfly_mixed_radix_recursive_dit(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: radix, nr, r, l

        if (N == 1) return

        radix   = stages(stage_idx)%radix
        nr = N / radix

        block
            real(sp), dimension(nr, radix) :: sub_re, sub_im
            real(sp), dimension(nr, radix) :: tmp_re_vec,   tmp_im_vec

            real(sp), dimension(radix, radix)   :: w_re,   w_im
            integer(isp) :: m, k
            real(sp)     :: angle

            do k = 0, radix - 1
                do m = 0, radix - 1
                    angle        = -2.0_sp * pi_sp * real(m * k, sp) / real(radix, sp)
                    w_re(m+1, k+1) = cos(angle)
                    w_im(m+1, k+1) = sin(angle)
                end do
            end do

            do r = 1, radix
                sub_re(:, r) = re(r:N:radix)
                sub_im(:, r) = im(r:N:radix)
            end do

            if (stage_idx > 1) then
                do r = 1, radix
                    call execute_butterfly_mixed_radix_recursive_dit( &
                        stages, nr, sub_re(:, r), sub_im(:, r), stage_idx - 1)
                end do
            end if

            tmp_re_vec(:, 1) = sub_re(:, 1)
            tmp_im_vec(:, 1) = sub_im(:, 1)
            do r = 2, radix
                associate( &
                    tw_r => stages(stage_idx)%twiddles_real((r-2)*nr+1 : (r-1)*nr), &
                    tw_i => stages(stage_idx)%twiddles_imag((r-2)*nr+1 : (r-1)*nr))
                    tmp_re_vec(:, r) = sub_re(:, r) * tw_r - sub_im(:, r) * tw_i
                    tmp_im_vec(:, r) = sub_re(:, r) * tw_i + sub_im(:, r) * tw_r
                end associate
            end do

            do l = 0, radix - 1
                re(l*nr+1 : (l+1)*nr) = 0.0_sp
                im(l*nr+1 : (l+1)*nr) = 0.0_sp
                do r = 1, radix
                    re(l*nr+1 : (l+1)*nr) = re(l*nr+1 : (l+1)*nr) &
                        + tmp_re_vec(:, r) * w_re(r, l+1) - tmp_im_vec(:, r) * w_im(r, l+1)
                    im(l*nr+1 : (l+1)*nr) = im(l*nr+1 : (l+1)*nr) &
                        + tmp_re_vec(:, r) * w_im(r, l+1) + tmp_im_vec(:, r) * w_re(r, l+1)
                end do
            end do
        end block

    end subroutine execute_butterfly_mixed_radix_recursive_dit

    pure recursive subroutine execute_butterfly_mixed_radix_recursive_dif(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: radix, nr, r, l, num_stages

        if (N == 1) return

        radix          = stages(stage_idx)%radix
        nr        = N / radix
        num_stages = size(stages)

        block
            real(sp), dimension(nr, radix) :: blk_re, blk_im
            real(sp), dimension(nr, radix) :: Y_re,   Y_im
            real(sp), dimension(radix, radix)   :: w_re,   w_im
            integer(isp) :: k, m
            real(sp)     :: angle

            do k = 0, radix - 1
                do m = 0, radix - 1
                    angle            = -2.0_sp * pi_sp * real(k * m, sp) / real(radix, sp)
                    w_re(m+1, k+1) = cos(angle)
                    w_im(m+1, k+1) = sin(angle)
                end do
            end do

            do r = 1, radix
                blk_re(:, r) = re((r-1)*nr+1 : r*nr)
                blk_im(:, r) = im((r-1)*nr+1 : r*nr)
            end do

            do l = 0, radix - 1
                Y_re(:, l+1) = 0.0_sp
                Y_im(:, l+1) = 0.0_sp
                do r = 1, radix
                    Y_re(:, l+1) = Y_re(:, l+1) &
                        + blk_re(:, r) * w_re(r, l+1) - blk_im(:, r) * w_im(r, l+1)
                    Y_im(:, l+1) = Y_im(:, l+1) &
                        + blk_re(:, r) * w_im(r, l+1) + blk_im(:, r) * w_re(r, l+1)
                end do
            end do

            blk_re(:, 1) = Y_re(:, 1)
            blk_im(:, 1) = Y_im(:, 1)
            do r = 2, radix
                associate( &
                    tw_r => stages(stage_idx)%twiddles_real((r-2)*nr+1 : (r-1)*nr), &
                    tw_i => stages(stage_idx)%twiddles_imag((r-2)*nr+1 : (r-1)*nr))
                    blk_re(:, r) = Y_re(:, r) * tw_r - Y_im(:, r) * tw_i
                    blk_im(:, r) = Y_re(:, r) * tw_i + Y_im(:, r) * tw_r
                end associate
            end do

            if (stage_idx < num_stages) then
                do r = 1, radix
                    call execute_butterfly_mixed_radix_recursive_dif( &
                        stages, nr, blk_re(:, r), blk_im(:, r), stage_idx + 1)
                end do
            end if

            do r = 1, radix
                re(r:N:radix) = blk_re(:, r)
                im(r:N:radix) = blk_im(:, r)
            end do
        end block

    end subroutine execute_butterfly_mixed_radix_recursive_dif

end submodule NAFPack_butterfly_mixed_radix