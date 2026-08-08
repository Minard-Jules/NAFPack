submodule(NAFPack_butterfly) NAFPack_butterfly_radix3

    use NAFPack_constant, only: SQRT3_HALF_sp, SQRT3_HALF_dp, SQRT3_HALF_qp

    implicit none(type, external)

contains

    pure module subroutine execute_butterfly_radix3_recursive(stages, N, re, im, stage_idx, decimation_method)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        if (stage_idx < 1 .or. stage_idx > size(stages)) &
            error stop "stage_idx out of bounds"

        select case (decimation_method%id)
        case (DIT%id)
            call execute_butterfly_radix3_recursive_dit(stages, N, re, im, stage_idx)
        case (DIF%id)
            call execute_butterfly_radix3_recursive_dif(stages, N, re, im, stage_idx)
        case default
            error stop "Error in execute_butterfly_radix3_recursive: unknown decimation method"
        end select

    end subroutine execute_butterfly_radix3_recursive

    pure recursive subroutine execute_butterfly_radix3_recursive_dit(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/3) :: re1, re2, re3, im1, im2, im3
        real(sp), dimension(N/3) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im
        real(sp), dimension(N/3) :: u1_re, u1_im, v1_re, v1_im
        real(sp), dimension(N/3) :: p_re, p_im, q_re, q_im
        integer(isp) :: n3

        if (N == 1) then
            return
        end if

        n3 = N / 3

        re1 = re(1:N:3)
        im1 = im(1:N:3)

        re2 = re(2:N:3)
        im2 = im(2:N:3)

        re3 = re(3:N:3)
        im3 = im(3:N:3)

        call execute_butterfly_radix3_recursive_dit(stages, n3, re1, im1, stage_idx - 1)
        call execute_butterfly_radix3_recursive_dit(stages, n3, re2, im2, stage_idx - 1)
        call execute_butterfly_radix3_recursive_dit(stages, n3, re3, im3, stage_idx - 1)

        associate (tw1_re => stages(stage_idx)%twiddles_real(1:n3), &
                   tw1_im => stages(stage_idx)%twiddles_imag(1:n3), &
                   tw2_re => stages(stage_idx)%twiddles_real(n3 + 1:2 * n3), &
                   tw2_im => stages(stage_idx)%twiddles_imag(n3 + 1:2 * n3))

            tmp1_re = re2 * tw1_re - im2 * tw1_im
            tmp1_im = re2 * tw1_im + im2 * tw1_re

            tmp2_re = re3 * tw2_re - im3 * tw2_im
            tmp2_im = re3 * tw2_im + im3 * tw2_re
        end associate

        u1_re = tmp1_re + tmp2_re; u1_im = tmp1_im + tmp2_im
        v1_re = tmp1_re - tmp2_re; v1_im = tmp1_im - tmp2_im

        p_re = re1 - 0.5_sp * u1_re; p_im = im1 - 0.5_sp * u1_im
        q_re = SQRT3_HALF_sp * v1_re; q_im = SQRT3_HALF_sp * v1_im

        re(1:n3) = re1 + u1_re
        im(1:n3) = im1 + u1_im

        re(n3 + 1:2 * n3) = p_re + q_im
        im(n3 + 1:2 * n3) = p_im - q_re

        re(2 * n3 + 1:N) = p_re - q_im
        im(2 * n3 + 1:N) = p_im + q_re

    end subroutine execute_butterfly_radix3_recursive_dit

    pure recursive subroutine execute_butterfly_radix3_recursive_dif(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/3) :: re1, re2, re3, im1, im2, im3
        real(sp), dimension(N/3) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im
        real(sp), dimension(N/3) :: u1_re, u1_im, v1_re, v1_im
        real(sp), dimension(N/3) :: p_re, p_im, q_re, q_im
        integer(isp) :: n3

        if (N == 1) then
            return
        end if

        n3 = N / 3

        re1 = re(1:n3)
        im1 = im(1:n3)

        re2 = re(n3 + 1:2 * n3)
        im2 = im(n3 + 1:2 * n3)

        re3 = re(2 * n3 + 1:N)
        im3 = im(2 * n3 + 1:N)

        u1_re = re2 + re3; u1_im = im2 + im3
        v1_re = re2 - re3; v1_im = im2 - im3

        p_re = re1 - 0.5_sp * u1_re; p_im = im1 - 0.5_sp * u1_im
        q_re = SQRT3_HALF_sp * v1_re; q_im = SQRT3_HALF_sp * v1_im

        tmp1_re = p_re + q_im; tmp1_im = p_im - q_re
        tmp2_re = p_re - q_im; tmp2_im = p_im + q_re

        re1 = re1 + u1_re
        im1 = im1 + u1_im

        associate (tw1_re => stages(stage_idx)%twiddles_real(1:n3), &
                   tw1_im => stages(stage_idx)%twiddles_imag(1:n3), &
                   tw2_re => stages(stage_idx)%twiddles_real(n3 + 1:2 * n3), &
                   tw2_im => stages(stage_idx)%twiddles_imag(n3 + 1:2 * n3))

            re2 = tmp1_re * tw1_re - tmp1_im * tw1_im
            im2 = tmp1_re * tw1_im + tmp1_im * tw1_re

            re3 = tmp2_re * tw2_re - tmp2_im * tw2_im
            im3 = tmp2_re * tw2_im + tmp2_im * tw2_re
        end associate

        call execute_butterfly_radix3_recursive_dif(stages, n3, re1, im1, stage_idx + 1)
        call execute_butterfly_radix3_recursive_dif(stages, n3, re2, im2, stage_idx + 1)
        call execute_butterfly_radix3_recursive_dif(stages, n3, re3, im3, stage_idx + 1)

        re(1:N:3) = re1
        im(1:N:3) = im1

        re(2:N:3) = re2
        im(2:N:3) = im2

        re(3:N:3) = re3
        im(3:N:3) = im3

    end subroutine execute_butterfly_radix3_recursive_dif

    pure module subroutine execute_butterfly_radix3(stage, N, re, im, decimation_method)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        !TODO : check loop method and decimation method validity
        !TODO : check power of three for radix 3
        select case (decimation_method%id)
        case (DIT%id)
            call execute_do_classique_dit(stage, N, re, im)
        case (DIF%id)
            call execute_do_classique_dif(stage, N, re, im)
        case default
            error stop "Error in execute_butterfly_radix3: unknown decimation method"
        end select

    end subroutine execute_butterfly_radix3

    pure module subroutine execute_do_classique_dit(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, idx3, num_groups, offset, step
        real(sp) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im
        real(sp) :: u1_re, u1_im, v1_re, v1_im
        real(sp) :: p_re, p_im, q_re, q_im
        real(sp) :: tw1_re, tw1_im, tw2_re, tw2_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride
                idx3 = idx2 + stage%stride

                tw1_re = stage%twiddles_real(twiddle_idx)
                tw1_im = stage%twiddles_imag(twiddle_idx)
                tw2_re = stage%twiddles_real(twiddle_idx + stage%stride)
                tw2_im = stage%twiddles_imag(twiddle_idx + stage%stride)

                tmp1_re = re(idx2) * tw1_re - im(idx2) * tw1_im
                tmp1_im = re(idx2) * tw1_im + im(idx2) * tw1_re

                tmp2_re = re(idx3) * tw2_re - im(idx3) * tw2_im
                tmp2_im = re(idx3) * tw2_im + im(idx3) * tw2_re

                u1_re = tmp1_re + tmp2_re; u1_im = tmp1_im + tmp2_im
                v1_re = tmp1_re - tmp2_re; v1_im = tmp1_im - tmp2_im

                p_re = re(idx1) - 0.5_sp * u1_re
                p_im = im(idx1) - 0.5_sp * u1_im

                q_re = SQRT3_HALF_sp * v1_re
                q_im = SQRT3_HALF_sp * v1_im

                re(idx1) = re(idx1) + u1_re
                im(idx1) = im(idx1) + u1_im

                re(idx2) = p_re + q_im
                im(idx2) = p_im - q_re

                re(idx3) = p_re - q_im
                im(idx3) = p_im + q_re
            end do
        end do

    end subroutine execute_do_classique_dit

    pure module subroutine execute_do_classique_dif(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, idx3, num_groups, offset, step
        real(sp) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im
        real(sp) :: u1_re, u1_im, v1_re, v1_im
        real(sp) :: p_re, p_im, q_re, q_im
        real(sp) :: tw1_re, tw1_im, tw2_re, tw2_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride
                idx3 = idx2 + stage%stride

                tw1_re = stage%twiddles_real(twiddle_idx)
                tw1_im = stage%twiddles_imag(twiddle_idx)
                tw2_re = stage%twiddles_real(twiddle_idx + stage%stride)
                tw2_im = stage%twiddles_imag(twiddle_idx + stage%stride)

                u1_re = re(idx2) + re(idx3)
                u1_im = im(idx2) + im(idx3)
                v1_re = re(idx2) - re(idx3)
                v1_im = im(idx2) - im(idx3)

                p_re = re(idx1) - 0.5_sp * u1_re
                p_im = im(idx1) - 0.5_sp * u1_im

                q_re = SQRT3_HALF_sp * v1_re
                q_im = SQRT3_HALF_sp * v1_im

                tmp1_re = p_re + q_im
                tmp1_im = p_im - q_re

                tmp2_re = p_re - q_im
                tmp2_im = p_im + q_re

                re(idx1) = re(idx1) + u1_re
                im(idx1) = im(idx1) + u1_im

                re(idx2) = tmp1_re * tw1_re - tmp1_im * tw1_im
                im(idx2) = tmp1_re * tw1_im + tmp1_im * tw1_re

                re(idx3) = tmp2_re * tw2_re - tmp2_im * tw2_im
                im(idx3) = tmp2_re * tw2_im + tmp2_im * tw2_re
            end do
        end do

    end subroutine execute_do_classique_dif

end submodule NAFPack_butterfly_radix3
