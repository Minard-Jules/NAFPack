submodule(NAFPack_butterfly) NAFPack_butterfly_radix5

    implicit none(type, external)

    real(sp), parameter :: C1 = 0.30901699437494745_sp ! cos(2*pi/5)
    real(sp), parameter :: C2 = 0.80901699437494745_sp ! cos(4*pi/5)
    real(sp), parameter :: S1 = 0.9510565162951535_sp ! sin(2*pi/5)
    real(sp), parameter :: S2 = 0.5877852522924731_sp ! sin(4*pi/5)

contains

    pure module subroutine execute_butterfly_radix5_recursive(stages, N, re, im, stage_idx, decimation_method)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        if (stage_idx < 1 .or. stage_idx > size(stages)) &
            error stop "stage_idx out of bounds"

        select case (decimation_method%id)
        case (DIT%id)
            call execute_butterfly_radix5_recursive_dit(stages, N, re, im, stage_idx)
        case (DIF%id)
            call execute_butterfly_radix5_recursive_dif(stages, N, re, im, stage_idx)
        case default
            error stop "Error in execute_butterfly_radix5_recursive: unknown decimation method"
        end select

    end subroutine execute_butterfly_radix5_recursive

    pure recursive subroutine execute_butterfly_radix5_recursive_dit(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/5) :: re1, re2, re3, re4, re5, im1, im2, im3, im4, im5
        real(sp), dimension(N/5) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im
        real(sp), dimension(N/5) :: tmp3_re, tmp3_im, tmp4_re, tmp4_im
        real(sp), dimension(N/5) :: u1_re, u1_im, u2_re, u2_im, v1_re, v1_im, v2_re, v2_im
        real(sp), dimension(N/5) :: p_re, p_im, q_re, q_im
        integer(isp) :: n5

        if (N == 1) return

        n5 = N / 5

        re1 = re(1:N:5); im1 = im(1:N:5)
        re2 = re(2:N:5); im2 = im(2:N:5)
        re3 = re(3:N:5); im3 = im(3:N:5)
        re4 = re(4:N:5); im4 = im(4:N:5)
        re5 = re(5:N:5); im5 = im(5:N:5)

        call execute_butterfly_radix5_recursive_dit(stages, n5, re1, im1, stage_idx - 1)
        call execute_butterfly_radix5_recursive_dit(stages, n5, re2, im2, stage_idx - 1)
        call execute_butterfly_radix5_recursive_dit(stages, n5, re3, im3, stage_idx - 1)
        call execute_butterfly_radix5_recursive_dit(stages, n5, re4, im4, stage_idx - 1)
        call execute_butterfly_radix5_recursive_dit(stages, n5, re5, im5, stage_idx - 1)

        associate (tw1_re => stages(stage_idx)%twiddles_real(1:n5), &
                   tw1_im => stages(stage_idx)%twiddles_imag(1:n5), &
                   tw2_re => stages(stage_idx)%twiddles_real(n5 + 1:2 * n5), &
                   tw2_im => stages(stage_idx)%twiddles_imag(n5 + 1:2 * n5), &
                   tw3_re => stages(stage_idx)%twiddles_real(2 * n5 + 1:3 * n5), &
                   tw3_im => stages(stage_idx)%twiddles_imag(2 * n5 + 1:3 * n5), &
                   tw4_re => stages(stage_idx)%twiddles_real(3 * n5 + 1:4 * n5), &
                   tw4_im => stages(stage_idx)%twiddles_imag(3 * n5 + 1:4 * n5))

            tmp1_re = re2 * tw1_re - im2 * tw1_im
            tmp1_im = re2 * tw1_im + im2 * tw1_re

            tmp2_re = re3 * tw2_re - im3 * tw2_im
            tmp2_im = re3 * tw2_im + im3 * tw2_re

            tmp3_re = re4 * tw3_re - im4 * tw3_im
            tmp3_im = re4 * tw3_im + im4 * tw3_re

            tmp4_re = re5 * tw4_re - im5 * tw4_im
            tmp4_im = re5 * tw4_im + im5 * tw4_re
        end associate

        u1_re = tmp1_re + tmp4_re; u1_im = tmp1_im + tmp4_im
        u2_re = tmp2_re + tmp3_re; u2_im = tmp2_im + tmp3_im
        v1_re = tmp1_re - tmp4_re; v1_im = tmp1_im - tmp4_im
        v2_re = tmp2_re - tmp3_re; v2_im = tmp2_im - tmp3_im

        p_re = C1 * u1_re - C2 * u2_re; p_im = C1 * u1_im - C2 * u2_im
        q_re = C1 * u2_re - C2 * u1_re; q_im = C1 * u2_im - C2 * u1_im

        re(1:n5) = re1 + u1_re + u2_re
        im(1:n5) = im1 + u1_im + u2_im

        re(n5 + 1:2 * n5) = re1 + p_re + S1 * v1_im + S2 * v2_im
        im(n5 + 1:2 * n5) = im1 + p_im - S1 * v1_re - S2 * v2_re

        re(2 * n5 + 1:3 * n5) = re1 + q_re + S2 * v1_im - S1 * v2_im
        im(2 * n5 + 1:3 * n5) = im1 + q_im - S2 * v1_re + S1 * v2_re

        re(3 * n5 + 1:4 * n5) = re1 + q_re - S2 * v1_im + S1 * v2_im
        im(3 * n5 + 1:4 * n5) = im1 + q_im + S2 * v1_re - S1 * v2_re

        re(4 * n5 + 1:N) = re1 + p_re - S1 * v1_im - S2 * v2_im
        im(4 * n5 + 1:N) = im1 + p_im + S1 * v1_re + S2 * v2_re

    end subroutine execute_butterfly_radix5_recursive_dit

    pure recursive subroutine execute_butterfly_radix5_recursive_dif(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/5) :: re1, re2, re3, re4, re5, im1, im2, im3, im4, im5
        real(sp), dimension(N/5) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im
        real(sp), dimension(N/5) :: tmp3_re, tmp3_im, tmp4_re, tmp4_im
        real(sp), dimension(N/5) :: u1_re, u1_im, u2_re, u2_im, v1_re, v1_im, v2_re, v2_im
        real(sp), dimension(N/5) :: p_re, p_im, q_re, q_im
        integer(isp) :: n5

        if (N == 1) return

        n5 = N / 5

        re1 = re(1:n5); im1 = im(1:n5)
        re2 = re(n5 + 1:2 * n5); im2 = im(n5 + 1:2 * n5)
        re3 = re(2 * n5 + 1:3 * n5); im3 = im(2 * n5 + 1:3 * n5)
        re4 = re(3 * n5 + 1:4 * n5); im4 = im(3 * n5 + 1:4 * n5)
        re5 = re(4 * n5 + 1:N); im5 = im(4 * n5 + 1:N)

        u1_re = re2 + re5; u1_im = im2 + im5
        u2_re = re3 + re4; u2_im = im3 + im4
        v1_re = re2 - re5; v1_im = im2 - im5
        v2_re = re3 - re4; v2_im = im3 - im4

        p_re = C1 * u1_re - C2 * u2_re; p_im = C1 * u1_im - C2 * u2_im
        q_re = C1 * u2_re - C2 * u1_re; q_im = C1 * u2_im - C2 * u1_im

        tmp1_re = re1 + p_re + S1 * v1_im + S2 * v2_im
        tmp1_im = im1 + p_im - S1 * v1_re - S2 * v2_re

        tmp2_re = re1 + q_re + S2 * v1_im - S1 * v2_im
        tmp2_im = im1 + q_im - S2 * v1_re + S1 * v2_re

        tmp3_re = re1 + q_re - S2 * v1_im + S1 * v2_im
        tmp3_im = im1 + q_im + S2 * v1_re - S1 * v2_re

        tmp4_re = re1 + p_re - S1 * v1_im - S2 * v2_im
        tmp4_im = im1 + p_im + S1 * v1_re + S2 * v2_re

        re1 = re1 + u1_re + u2_re
        im1 = im1 + u1_im + u2_im

        associate (tw1_re => stages(stage_idx)%twiddles_real(1:n5), &
                   tw1_im => stages(stage_idx)%twiddles_imag(1:n5), &
                   tw2_re => stages(stage_idx)%twiddles_real(n5 + 1:2 * n5), &
                   tw2_im => stages(stage_idx)%twiddles_imag(n5 + 1:2 * n5), &
                   tw3_re => stages(stage_idx)%twiddles_real(2 * n5 + 1:3 * n5), &
                   tw3_im => stages(stage_idx)%twiddles_imag(2 * n5 + 1:3 * n5), &
                   tw4_re => stages(stage_idx)%twiddles_real(3 * n5 + 1:4 * n5), &
                   tw4_im => stages(stage_idx)%twiddles_imag(3 * n5 + 1:4 * n5))

            re2 = tmp1_re * tw1_re - tmp1_im * tw1_im
            im2 = tmp1_re * tw1_im + tmp1_im * tw1_re

            re3 = tmp2_re * tw2_re - tmp2_im * tw2_im
            im3 = tmp2_re * tw2_im + tmp2_im * tw2_re

            re4 = tmp3_re * tw3_re - tmp3_im * tw3_im
            im4 = tmp3_re * tw3_im + tmp3_im * tw3_re

            re5 = tmp4_re * tw4_re - tmp4_im * tw4_im
            im5 = tmp4_re * tw4_im + tmp4_im * tw4_re
        end associate

        call execute_butterfly_radix5_recursive_dif(stages, n5, re1, im1, stage_idx + 1)
        call execute_butterfly_radix5_recursive_dif(stages, n5, re2, im2, stage_idx + 1)
        call execute_butterfly_radix5_recursive_dif(stages, n5, re3, im3, stage_idx + 1)
        call execute_butterfly_radix5_recursive_dif(stages, n5, re4, im4, stage_idx + 1)
        call execute_butterfly_radix5_recursive_dif(stages, n5, re5, im5, stage_idx + 1)

        re(1:N:5) = re1; im(1:N:5) = im1
        re(2:N:5) = re2; im(2:N:5) = im2
        re(3:N:5) = re3; im(3:N:5) = im3
        re(4:N:5) = re4; im(4:N:5) = im4
        re(5:N:5) = re5; im(5:N:5) = im5

    end subroutine execute_butterfly_radix5_recursive_dif

    pure module subroutine execute_butterfly_radix5(stage, N, re, im, decimation_method)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        !TODO : check loop method and decimation method validity
        !TODO : check power of five for radix 5
        select case (decimation_method%id)
        case (DIT%id)
            call execute_do_classique_dit(stage, N, re, im)
        case (DIF%id)
            call execute_do_classique_dif(stage, N, re, im)
        case default
            error stop "Error in execute_butterfly_radix5: unknown decimation method"
        end select

    end subroutine execute_butterfly_radix5

    pure module subroutine execute_do_classique_dit(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, idx3, idx4, idx5, num_groups, offset, step
        real(sp) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im, tmp3_re, tmp3_im, tmp4_re, tmp4_im
        real(sp) :: u1_re, u1_im, u2_re, u2_im, v1_re, v1_im, v2_re, v2_im
        real(sp) :: p_re, p_im, q_re, q_im
        real(sp) :: tw1_re, tw1_im, tw2_re, tw2_im, tw3_re, tw3_im, tw4_re, tw4_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride
                idx3 = idx2 + stage%stride
                idx4 = idx3 + stage%stride
                idx5 = idx4 + stage%stride

                tw1_re = stage%twiddles_real(twiddle_idx)
                tw1_im = stage%twiddles_imag(twiddle_idx)
                tw2_re = stage%twiddles_real(twiddle_idx + stage%stride)
                tw2_im = stage%twiddles_imag(twiddle_idx + stage%stride)
                tw3_re = stage%twiddles_real(twiddle_idx + 2 * stage%stride)
                tw3_im = stage%twiddles_imag(twiddle_idx + 2 * stage%stride)
                tw4_re = stage%twiddles_real(twiddle_idx + 3 * stage%stride)
                tw4_im = stage%twiddles_imag(twiddle_idx + 3 * stage%stride)

                tmp1_re = re(idx2) * tw1_re - im(idx2) * tw1_im
                tmp1_im = re(idx2) * tw1_im + im(idx2) * tw1_re

                tmp2_re = re(idx3) * tw2_re - im(idx3) * tw2_im
                tmp2_im = re(idx3) * tw2_im + im(idx3) * tw2_re

                tmp3_re = re(idx4) * tw3_re - im(idx4) * tw3_im
                tmp3_im = re(idx4) * tw3_im + im(idx4) * tw3_re

                tmp4_re = re(idx5) * tw4_re - im(idx5) * tw4_im
                tmp4_im = re(idx5) * tw4_im + im(idx5) * tw4_re

                u1_re = tmp1_re + tmp4_re; u1_im = tmp1_im + tmp4_im
                u2_re = tmp2_re + tmp3_re; u2_im = tmp2_im + tmp3_im
                v1_re = tmp1_re - tmp4_re; v1_im = tmp1_im - tmp4_im
                v2_re = tmp2_re - tmp3_re; v2_im = tmp2_im - tmp3_im

                p_re = re(idx1) + C1 * u1_re - C2 * u2_re; p_im = im(idx1) + C1 * u1_im - C2 * u2_im
                q_re = re(idx1) + C1 * u2_re - C2 * u1_re; q_im = im(idx1) + C1 * u2_im - C2 * u1_im

                re(idx1) = re(idx1) + u1_re + u2_re
                im(idx1) = im(idx1) + u1_im + u2_im

                re(idx2) = p_re + S1 * v1_im + S2 * v2_im
                im(idx2) = p_im - S1 * v1_re - S2 * v2_re

                re(idx3) = q_re + S2 * v1_im - S1 * v2_im
                im(idx3) = q_im - S2 * v1_re + S1 * v2_re

                re(idx4) = q_re - S2 * v1_im + S1 * v2_im
                im(idx4) = q_im + S2 * v1_re - S1 * v2_re

                re(idx5) = p_re - S1 * v1_im - S2 * v2_im
                im(idx5) = p_im + S1 * v1_re + S2 * v2_re
            end do
        end do

    end subroutine execute_do_classique_dit

    pure module subroutine execute_do_classique_dif(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, idx3, idx4, idx5, num_groups, offset, step
        real(sp) :: u1_re, u1_im, u2_re, u2_im, v1_re, v1_im, v2_re, v2_im
        real(sp) :: p_re, p_im, q_re, q_im
        real(sp) :: tmp1_re, tmp1_im, tmp2_re, tmp2_im, tmp3_re, tmp3_im, tmp4_re, tmp4_im
        real(sp) :: tw1_re, tw1_im, tw2_re, tw2_im, tw3_re, tw3_im, tw4_re, tw4_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride
                idx3 = idx2 + stage%stride
                idx4 = idx3 + stage%stride
                idx5 = idx4 + stage%stride

                tw1_re = stage%twiddles_real(twiddle_idx)
                tw1_im = stage%twiddles_imag(twiddle_idx)
                tw2_re = stage%twiddles_real(twiddle_idx + stage%stride)
                tw2_im = stage%twiddles_imag(twiddle_idx + stage%stride)
                tw3_re = stage%twiddles_real(twiddle_idx + 2 * stage%stride)
                tw3_im = stage%twiddles_imag(twiddle_idx + 2 * stage%stride)
                tw4_re = stage%twiddles_real(twiddle_idx + 3 * stage%stride)
                tw4_im = stage%twiddles_imag(twiddle_idx + 3 * stage%stride)

                u1_re = re(idx2) + re(idx5)
                u1_im = im(idx2) + im(idx5)
                u2_re = re(idx3) + re(idx4)
                u2_im = im(idx3) + im(idx4)
                v1_re = re(idx2) - re(idx5)
                v1_im = im(idx2) - im(idx5)
                v2_re = re(idx3) - re(idx4)
                v2_im = im(idx3) - im(idx4)

                p_re = re(idx1) + C1 * u1_re - C2 * u2_re; p_im = im(idx1) + C1 * u1_im - C2 * u2_im
                q_re = re(idx1) + C1 * u2_re - C2 * u1_re; q_im = im(idx1) + C1 * u2_im - C2 * u1_im

                tmp1_re = p_re + S1 * v1_im + S2 * v2_im
                tmp1_im = p_im - S1 * v1_re - S2 * v2_re

                tmp2_re = q_re + S2 * v1_im - S1 * v2_im
                tmp2_im = q_im - S2 * v1_re + S1 * v2_re

                tmp3_re = q_re - S2 * v1_im + S1 * v2_im
                tmp3_im = q_im + S2 * v1_re - S1 * v2_re

                tmp4_re = p_re - S1 * v1_im - S2 * v2_im
                tmp4_im = p_im + S1 * v1_re + S2 * v2_re

                re(idx1) = re(idx1) + u1_re + u2_re
                im(idx1) = im(idx1) + u1_im + u2_im

                re(idx2) = tmp1_re * tw1_re - tmp1_im * tw1_im
                im(idx2) = tmp1_re * tw1_im + tmp1_im * tw1_re

                re(idx3) = tmp2_re * tw2_re - tmp2_im * tw2_im
                im(idx3) = tmp2_re * tw2_im + tmp2_im * tw2_re

                re(idx4) = tmp3_re * tw3_re - tmp3_im * tw3_im
                im(idx4) = tmp3_re * tw3_im + tmp3_im * tw3_re

                re(idx5) = tmp4_re * tw4_re - tmp4_im * tw4_im
                im(idx5) = tmp4_re * tw4_im + tmp4_im * tw4_re
            end do
        end do

    end subroutine execute_do_classique_dif

end submodule NAFPack_butterfly_radix5
