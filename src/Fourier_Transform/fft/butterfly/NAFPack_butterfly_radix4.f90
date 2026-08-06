submodule(NAFPack_butterfly) NAFPack_butterfly_radix4

    implicit none(type, external)

contains

    pure module subroutine execute_butterfly_radix4_recursive(stages, N, re, im, stage_idx, decimation_method)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        if (stage_idx < 1 .or. stage_idx > size(stages)) &
            error stop "stage_idx out of bounds"

        select case (decimation_method%id)
        case (DIT%id)
            call execute_butterfly_radix4_recursive_dit(stages, N, re, im, stage_idx)
        case (DIF%id)
            call execute_butterfly_radix4_recursive_dif(stages, N, re, im, stage_idx)
        case default
            error stop "Error in execute_butterfly_radix4_recursive: unknown decimation method"
        end select

    end subroutine execute_butterfly_radix4_recursive

    pure recursive subroutine execute_butterfly_radix4_recursive_dit(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/4) :: re1, re2, re3, re4, im1, im2, im3, im4
        real(sp), dimension(N/4) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im, tmp3_re, tmp3_im
        real(sp), dimension(N/4) :: u1_re, u1_im, u2_re, u2_im, v1_re, v1_im, v2_re, v2_im
        integer(isp) :: n4

        if (N == 1) then
            return
        end if

        n4 = N / 4

        re1 = re(1:N:4); im1 = im(1:N:4)
        re2 = re(2:N:4); im2 = im(2:N:4)
        re3 = re(3:N:4); im3 = im(3:N:4)
        re4 = re(4:N:4); im4 = im(4:N:4)

        call execute_butterfly_radix4_recursive_dit(stages, n4, re1, im1, stage_idx - 1)
        call execute_butterfly_radix4_recursive_dit(stages, n4, re2, im2, stage_idx - 1)
        call execute_butterfly_radix4_recursive_dit(stages, n4, re3, im3, stage_idx - 1)
        call execute_butterfly_radix4_recursive_dit(stages, n4, re4, im4, stage_idx - 1)

        associate (tw1_re => stages(stage_idx)%twiddles_real(1:n4), &
                   tw1_im => stages(stage_idx)%twiddles_imag(1:n4), &
                   tw2_re => stages(stage_idx)%twiddles_real(n4 + 1:2 * n4), &
                   tw2_im => stages(stage_idx)%twiddles_imag(n4 + 1:2 * n4), &
                   tw3_re => stages(stage_idx)%twiddles_real(2 * n4 + 1:3 * n4), &
                   tw3_im => stages(stage_idx)%twiddles_imag(2 * n4 + 1:3 * n4))

            tmp1_re = re2 * tw1_re - im2 * tw1_im
            tmp1_im = re2 * tw1_im + im2 * tw1_re

            tmp2_re = re3 * tw2_re - im3 * tw2_im
            tmp2_im = re3 * tw2_im + im3 * tw2_re

            tmp3_re = re4 * tw3_re - im4 * tw3_im
            tmp3_im = re4 * tw3_im + im4 * tw3_re
        end associate

        u1_re = re1 + tmp2_re; u1_im = im1 + tmp2_im
        u2_re = re1 - tmp2_re; u2_im = im1 - tmp2_im
        v1_re = tmp1_re + tmp3_re; v1_im = tmp1_im + tmp3_im
        v2_re = tmp1_re - tmp3_re; v2_im = tmp1_im - tmp3_im

        re(1:n4) = u1_re + v1_re
        im(1:n4) = u1_im + v1_im

        re(n4 + 1:2 * n4) = u2_re + v2_im
        im(n4 + 1:2 * n4) = u2_im - v2_re

        re(2 * n4 + 1:3 * n4) = u1_re - v1_re
        im(2 * n4 + 1:3 * n4) = u1_im - v1_im

        re(3 * n4 + 1:4 * n4) = u2_re - v2_im
        im(3 * n4 + 1:4 * n4) = u2_im + v2_re

    end subroutine execute_butterfly_radix4_recursive_dit

    pure recursive subroutine execute_butterfly_radix4_recursive_dif(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/4) :: re1, re2, re3, re4, im1, im2, im3, im4
        real(sp), dimension(N/4) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im, tmp3_re, tmp3_im
        real(sp), dimension(N/4) :: u1_re, u1_im, u2_re, u2_im, v1_re, v1_im, v2_re, v2_im
        integer(isp) :: n4

        if (N == 1) then
            return
        end if

        n4 = N / 4

        re1 = re(1:n4); im1 = im(1:n4)
        re2 = re(n4 + 1:2 * n4); im2 = im(n4 + 1:2 * n4)
        re3 = re(2 * n4 + 1:3 * n4); im3 = im(2 * n4 + 1:3 * n4)
        re4 = re(3 * n4 + 1:4 * n4); im4 = im(3 * n4 + 1:4 * n4)

        u1_re = re1 + re3; u1_im = im1 + im3
        u2_re = re1 - re3; u2_im = im1 - im3
        v1_re = re2 + re4; v1_im = im2 + im4
        v2_re = re2 - re4; v2_im = im2 - im4

        tmp1_re = u2_re + v2_im; tmp1_im = u2_im - v2_re
        tmp2_re = u1_re - v1_re; tmp2_im = u1_im - v1_im
        tmp3_re = u2_re - v2_im; tmp3_im = u2_im + v2_re

        re1 = u1_re + v1_re
        im1 = u1_im + v1_im

        associate (tw1_re => stages(stage_idx)%twiddles_real(1:n4), &
                   tw1_im => stages(stage_idx)%twiddles_imag(1:n4), &
                   tw2_re => stages(stage_idx)%twiddles_real(n4 + 1:2 * n4), &
                   tw2_im => stages(stage_idx)%twiddles_imag(n4 + 1:2 * n4), &
                   tw3_re => stages(stage_idx)%twiddles_real(2 * n4 + 1:3 * n4), &
                   tw3_im => stages(stage_idx)%twiddles_imag(2 * n4 + 1:3 * n4))

            re2 = tmp1_re * tw1_re - tmp1_im * tw1_im
            im2 = tmp1_re * tw1_im + tmp1_im * tw1_re

            re3 = tmp2_re * tw2_re - tmp2_im * tw2_im
            im3 = tmp2_re * tw2_im + tmp2_im * tw2_re

            re4 = tmp3_re * tw3_re - tmp3_im * tw3_im
            im4 = tmp3_re * tw3_im + tmp3_im * tw3_re
        end associate

        call execute_butterfly_radix4_recursive_dif(stages, n4, re1, im1, stage_idx + 1)
        call execute_butterfly_radix4_recursive_dif(stages, n4, re2, im2, stage_idx + 1)
        call execute_butterfly_radix4_recursive_dif(stages, n4, re3, im3, stage_idx + 1)
        call execute_butterfly_radix4_recursive_dif(stages, n4, re4, im4, stage_idx + 1)

        re(1:N:4) = re1
        im(1:N:4) = im1

        re(2:N:4) = re2
        im(2:N:4) = im2

        re(3:N:4) = re3
        im(3:N:4) = im3

        re(4:N:4) = re4
        im(4:N:4) = im4

    end subroutine execute_butterfly_radix4_recursive_dif

    pure module subroutine execute_butterfly_radix4(stage, N, re, im, decimation_method)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        !TODO : check loop method and decimation method validity
        !TODO : check power of four for radix 4
        select case (decimation_method%id)
        case (DIT%id)
            call execute_do_classique_dit(stage, N, re, im)
        case (DIF%id)
            call execute_do_classique_dif(stage, N, re, im)
        case default
            error stop "Error in execute_butterfly_radix4: unknown decimation method"
        end select

    end subroutine execute_butterfly_radix4

    pure module subroutine execute_do_classique_dit(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, idx3, idx4, num_groups, offset, step
        real(sp) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im, tmp3_re, tmp3_im
        real(sp) :: u1_re, u1_im, u2_re, u2_im, v1_re, v1_im, v2_re, v2_im
        real(sp) :: tw1_re, tw1_im, tw2_re, tw2_im, tw3_re, tw3_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride
                idx3 = idx2 + stage%stride
                idx4 = idx3 + stage%stride

                tw1_re = stage%twiddles_real(twiddle_idx)
                tw1_im = stage%twiddles_imag(twiddle_idx)
                tw2_re = stage%twiddles_real(twiddle_idx + stage%stride)
                tw2_im = stage%twiddles_imag(twiddle_idx + stage%stride)
                tw3_re = stage%twiddles_real(twiddle_idx + 2 * stage%stride)
                tw3_im = stage%twiddles_imag(twiddle_idx + 2 * stage%stride)

                tmp1_re = re(idx2) * tw1_re - im(idx2) * tw1_im
                tmp1_im = re(idx2) * tw1_im + im(idx2) * tw1_re

                tmp2_re = re(idx3) * tw2_re - im(idx3) * tw2_im
                tmp2_im = re(idx3) * tw2_im + im(idx3) * tw2_re

                tmp3_re = re(idx4) * tw3_re - im(idx4) * tw3_im
                tmp3_im = re(idx4) * tw3_im + im(idx4) * tw3_re

                u1_re = re(idx1) + tmp2_re; u1_im = im(idx1) + tmp2_im
                u2_re = re(idx1) - tmp2_re; u2_im = im(idx1) - tmp2_im
                v1_re = tmp1_re + tmp3_re; v1_im = tmp1_im + tmp3_im
                v2_re = tmp1_re - tmp3_re; v2_im = tmp1_im - tmp3_im

                re(idx1) = u1_re + v1_re
                im(idx1) = u1_im + v1_im

                re(idx2) = u2_re + v2_im
                im(idx2) = u2_im - v2_re

                re(idx3) = u1_re - v1_re
                im(idx3) = u1_im - v1_im

                re(idx4) = u2_re - v2_im
                im(idx4) = u2_im + v2_re
            end do
        end do

    end subroutine execute_do_classique_dit

    pure module subroutine execute_do_classique_dif(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, idx3, idx4, num_groups, offset, step
        real(sp) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im, tmp3_re, tmp3_im
        real(sp) :: u1_re, u1_im, u2_re, u2_im, v1_re, v1_im, v2_re, v2_im
        real(sp) :: tw1_re, tw1_im, tw2_re, tw2_im, tw3_re, tw3_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride
                idx3 = idx2 + stage%stride
                idx4 = idx3 + stage%stride

                tw1_re = stage%twiddles_real(twiddle_idx)
                tw1_im = stage%twiddles_imag(twiddle_idx)
                tw2_re = stage%twiddles_real(twiddle_idx + stage%stride)
                tw2_im = stage%twiddles_imag(twiddle_idx + stage%stride)
                tw3_re = stage%twiddles_real(twiddle_idx + 2 * stage%stride)
                tw3_im = stage%twiddles_imag(twiddle_idx + 2 * stage%stride)

                u1_re = re(idx1) + re(idx3); u1_im = im(idx1) + im(idx3)
                u2_re = re(idx1) - re(idx3); u2_im = im(idx1) - im(idx3)
                v1_re = re(idx2) + re(idx4); v1_im = im(idx2) + im(idx4)
                v2_re = re(idx2) - re(idx4); v2_im = im(idx2) - im(idx4)

                tmp1_re = u2_re + v2_im; tmp1_im = u2_im - v2_re
                tmp2_re = u1_re - v1_re; tmp2_im = u1_im - v1_im
                tmp3_re = u2_re - v2_im; tmp3_im = u2_im + v2_re

                re(idx1) = u1_re + v1_re
                im(idx1) = u1_im + v1_im

                re(idx2) = tmp1_re * tw1_re - tmp1_im * tw1_im
                im(idx2) = tmp1_re * tw1_im + tmp1_im * tw1_re

                re(idx3) = tmp2_re * tw2_re - tmp2_im * tw2_im
                im(idx3) = tmp2_re * tw2_im + tmp2_im * tw2_re

                re(idx4) = tmp3_re * tw3_re - tmp3_im * tw3_im
                im(idx4) = tmp3_re * tw3_im + tmp3_im * tw3_re
            end do
        end do

    end subroutine execute_do_classique_dif

end submodule NAFPack_butterfly_radix4
