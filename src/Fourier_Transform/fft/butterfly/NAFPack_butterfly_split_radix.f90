submodule(NAFPack_butterfly) NAFPack_butterfly_split_radix

    use NAFPack_constant, only: pi_sp

    implicit none(type, external)

contains

    pure recursive module subroutine execute_butterfly_split_radix_recursive(stages, N, re, im, stage_idx, decimation_method)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        select case (decimation_method%id)
        case (DIT%id)
            call execute_butterfly_split_radix_recursive_dit(stages, N, re, im, stage_idx)
        case (DIF%id)
            call execute_butterfly_split_radix_recursive_dif(stages, N, re, im, stage_idx)
        case default
            error stop "Error in execute_butterfly_split_radix_recursive: unknown decimation method"
        end select

    end subroutine execute_butterfly_split_radix_recursive

    pure recursive subroutine execute_butterfly_split_radix_recursive_dit(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/2) :: e_re, e_im
        real(sp), dimension(N/4) :: o1_re, o1_im, o3_re, o3_im
        real(sp), dimension(N/4) :: temp1_re, temp1_im, temp2_re, temp2_im
        real(sp), dimension(N/4) :: u1_re, u1_im, td_re, td_im
        integer(isp) :: n2, n4

        if (N == 1) return

        if (N == 2) then
            block
                real(sp) :: tmp_re, tmp_im
                tmp_re = re(1)
                tmp_im = im(1)
                re(1) = tmp_re + re(2)
                im(1) = tmp_im + im(2)
                re(2) = tmp_re - re(2)
                im(2) = tmp_im - im(2)
            end block
            return
        end if

        n2 = N / 2
        n4 = N / 4

        e_re = re(1:N:2)
        e_im = im(1:N:2)
        o1_re = re(2:N:4)
        o1_im = im(2:N:4)
        o3_re = re(4:N:4)
        o3_im = im(4:N:4)

        call execute_butterfly_split_radix_recursive_dit(stages, n2, e_re, e_im, stage_idx - 1)
        call execute_butterfly_split_radix_recursive_dit(stages, n4, o1_re, o1_im, stage_idx - 2)
        call execute_butterfly_split_radix_recursive_dit(stages, n4, o3_re, o3_im, stage_idx - 2)

        associate (tw1_re => stages(stage_idx)%twiddles_real(1:n4), &
                   tw1_im => stages(stage_idx)%twiddles_imag(1:n4), &
                   tw3_re => stages(stage_idx)%twiddles_real(n4 + 1:n2), &
                   tw3_im => stages(stage_idx)%twiddles_imag(n4 + 1:n2))

            temp1_re = o1_re * tw1_re - o1_im * tw1_im
            temp1_im = o1_re * tw1_im + o1_im * tw1_re

            temp2_re = o3_re * tw3_re - o3_im * tw3_im
            temp2_im = o3_re * tw3_im + o3_im * tw3_re
        end associate

        u1_re = temp1_re + temp2_re
        u1_im = temp1_im + temp2_im
        td_re = temp1_re - temp2_re
        td_im = temp1_im - temp2_im

        re(1:n4) = e_re(1:n4) + u1_re(1:n4)
        im(1:n4) = e_im(1:n4) + u1_im(1:n4)

        re(n4 + 1:n2) = e_re(n4 + 1:n2) + td_im(1:n4)
        im(n4 + 1:n2) = e_im(n4 + 1:n2) - td_re(1:n4)

        re(n2 + 1:3 * n4) = e_re(1:n4) - u1_re(1:n4)
        im(n2 + 1:3 * n4) = e_im(1:n4) - u1_im(1:n4)

        re(3 * n4 + 1:N) = e_re(n4 + 1:n2) - td_im(1:n4)
        im(3 * n4 + 1:N) = e_im(n4 + 1:n2) + td_re(1:n4)

    end subroutine execute_butterfly_split_radix_recursive_dit

    pure recursive subroutine execute_butterfly_split_radix_recursive_dif(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/2) :: e_re, e_im
        real(sp), dimension(N/4) :: o1_re, o1_im, o3_re, o3_im
        real(sp), dimension(N/4) :: temp1_re, temp1_im, temp2_re, temp2_im
        real(sp), dimension(N/4) :: u1_re, u1_im, td_re, td_im
        integer(isp) :: n2, n4

        if (N == 1) return

        if (N == 2) then
            block
                real(sp) :: tmp_re, tmp_im
                tmp_re = re(1)
                tmp_im = im(1)
                re(1) = tmp_re + re(2)
                im(1) = tmp_im + im(2)
                re(2) = tmp_re - re(2)
                im(2) = tmp_im - im(2)
            end block
            return
        end if

        n2 = N / 2
        n4 = N / 4

        e_re = re(1:n2)
        e_im = im(1:n2)
        o1_re = re(n2 + 1:n2 + n4)
        o1_im = im(n2 + 1:n2 + n4)
        o3_re = re(n2 + n4 + 1:N)
        o3_im = im(n2 + n4 + 1:N)

        u1_re = e_re(1:n4) - o1_re
        u1_im = e_im(1:n4) - o1_im
        td_re = e_re(n4 + 1:n2) - o3_re
        td_im = e_im(n4 + 1:n2) - o3_im

        e_re(1:n4) = e_re(1:n4) + o1_re
        e_im(1:n4) = e_im(1:n4) + o1_im
        e_re(n4 + 1:n2) = e_re(n4 + 1:n2) + o3_re
        e_im(n4 + 1:n2) = e_im(n4 + 1:n2) + o3_im

        temp1_re = u1_re - td_im
        temp1_im = u1_re + td_im
        temp2_re = u1_im - td_re
        temp2_im = u1_im + td_re

        associate (tw1_re => stages(stage_idx)%twiddles_real(1:n4), &
                   tw1_im => stages(stage_idx)%twiddles_imag(1:n4), &
                   tw3_re => stages(stage_idx)%twiddles_real(n4 + 1:n2), &
                   tw3_im => stages(stage_idx)%twiddles_imag(n4 + 1:n2))
            o1_re = temp1_im * tw1_re - temp2_re * tw1_im
            o1_im = -temp2_re * tw1_re - temp1_im * tw1_im
            o3_re = temp1_re * tw3_re + temp2_im * tw3_im
            o3_im = temp2_im * tw3_re - temp1_re * tw3_im
        end associate

        call execute_butterfly_split_radix_recursive_dif(stages, n2, e_re, e_im, stage_idx + 1)
        call execute_butterfly_split_radix_recursive_dif(stages, n4, o1_re, o1_im, stage_idx + 2)
        call execute_butterfly_split_radix_recursive_dif(stages, n4, o3_re, o3_im, stage_idx + 2)

        re(1:N:2) = e_re
        im(1:N:2) = e_im
        re(2:N:4) = o1_re
        im(2:N:4) = o1_im
        re(4:N:4) = o3_re
        im(4:N:4) = o3_im

    end subroutine execute_butterfly_split_radix_recursive_dif

    module subroutine execute_butterfly_split_radix(stage, N, re, im, decimation_method)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        !TODO : check loop method and decimation method validity
        !TODO : check power of 2 and 4 for split-radix
        select case (decimation_method%id)
        case (DIT%id)
            call execute_do_classique_dit(stage, N, re, im)
        case (DIF%id)
            call execute_do_classique_dif(stage, N, re, im)
        case default
            error stop "Error in execute_butterfly_split_radix: unknown decimation method"
        end select

    end subroutine execute_butterfly_split_radix

    pure module subroutine execute_do_classique_dit(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: depth, group_idx, twiddle_idx, idx1, idx2, idx3, idx4
        integer(isp) :: base_offset, step
        real(sp) :: tw1_re, tw1_im, tw3_re, tw3_im
        real(sp) :: temp1_re, temp1_im, temp2_re, temp2_im
        real(sp) :: u1_re, u1_im, v1_re, v1_im, tmp_re, tmp_im

        if (stage%stride < 1) then
            do depth = 0, size(stage%groups%start_index) - 1
                base_offset = stage%groups%start_index(depth)
                step = stage%groups%group_stride(depth)

                do group_idx = base_offset, N - 1, step
                    idx1 = group_idx + 1
                    idx2 = idx1 + 1
                    if (idx2 > N) cycle

                    tmp_re = re(idx1)
                    tmp_im = im(idx1)

                    re(idx1) = tmp_re + re(idx2)
                    im(idx1) = tmp_im + im(idx2)

                    re(idx2) = tmp_re - re(idx2)
                    im(idx2) = tmp_im - im(idx2)
                end do
            end do
            return
        end if

        do depth = 0, size(stage%groups%start_index) - 1
            base_offset = stage%groups%start_index(depth)
            step = stage%groups%group_stride(depth)

            do group_idx = base_offset, N - 1, step
                do twiddle_idx = 1, stage%stride
                    idx1 = group_idx + twiddle_idx
                    idx2 = idx1 + stage%stride
                    idx3 = idx2 + stage%stride
                    idx4 = idx3 + stage%stride
                    if (idx4 > N) exit

                    tw1_re = stage%twiddles_real(twiddle_idx)
                    tw1_im = stage%twiddles_imag(twiddle_idx)
                    tw3_re = stage%twiddles_real(stage%stride + twiddle_idx)
                    tw3_im = stage%twiddles_imag(stage%stride + twiddle_idx)

                    temp1_re = re(idx3) * tw1_re - im(idx3) * tw1_im
                    temp1_im = re(idx3) * tw1_im + im(idx3) * tw1_re

                    temp2_re = re(idx4) * tw3_re - im(idx4) * tw3_im
                    temp2_im = re(idx4) * tw3_im + im(idx4) * tw3_re

                    u1_re = temp1_re + temp2_re
                    u1_im = temp1_im + temp2_im
                    v1_re = temp1_re - temp2_re
                    v1_im = temp1_im - temp2_im

                    temp1_re = re(idx1); temp1_im = im(idx1)
                    temp2_re = re(idx2); temp2_im = im(idx2)

                    re(idx1) = temp1_re + u1_re
                    im(idx1) = temp1_im + u1_im

                    re(idx2) = temp2_re + v1_im
                    im(idx2) = temp2_im - v1_re

                    re(idx3) = temp1_re - u1_re
                    im(idx3) = temp1_im - u1_im

                    re(idx4) = temp2_re - v1_im
                    im(idx4) = temp2_im + v1_re
                end do
            end do
        end do

    end subroutine execute_do_classique_dit

    pure module subroutine execute_do_classique_dif(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: depth, group_idx, twiddle_idx, idx1, idx2, idx3, idx4
        integer(isp) :: base_offset, step
        real(sp) :: tw1_re, tw1_im, tw3_re, tw3_im
        real(sp) :: tc_re, tc_im, td_re, td_im
        real(sp) :: u1_re, u1_im, v1_re, v1_im, tmp_re, tmp_im

        if (stage%stride < 1) then
            do depth = 0, size(stage%groups%start_index) - 1
                base_offset = stage%groups%start_index(depth)
                step = stage%groups%group_stride(depth)

                do group_idx = base_offset, N - 1, step
                    idx1 = group_idx + 1
                    idx2 = idx1 + 1
                    if (idx2 > N) cycle

                    tmp_re = re(idx1)
                    tmp_im = im(idx1)

                    re(idx1) = tmp_re + re(idx2)
                    im(idx1) = tmp_im + im(idx2)

                    re(idx2) = tmp_re - re(idx2)
                    im(idx2) = tmp_im - im(idx2)
                end do
            end do
            return
        end if

        do depth = 0, size(stage%groups%start_index) - 1
            base_offset = stage%groups%start_index(depth)
            step = stage%groups%group_stride(depth)

            do group_idx = base_offset, N - 1, step
                do twiddle_idx = 1, stage%stride
                    idx1 = group_idx + twiddle_idx
                    idx2 = idx1 + stage%stride
                    idx3 = idx2 + stage%stride
                    idx4 = idx3 + stage%stride
                    if (idx4 > N) exit

                    tw1_re = stage%twiddles_real(twiddle_idx)
                    tw1_im = stage%twiddles_imag(twiddle_idx)
                    tw3_re = stage%twiddles_real(stage%stride + twiddle_idx)
                    tw3_im = stage%twiddles_imag(stage%stride + twiddle_idx)

                    tc_re = re(idx1) - re(idx3)
                    tc_im = im(idx1) - im(idx3)
                    td_re = re(idx2) - re(idx4)
                    td_im = im(idx2) - im(idx4)

                    re(idx1) = re(idx1) + re(idx3)
                    im(idx1) = im(idx1) + im(idx3)
                    re(idx2) = re(idx2) + re(idx4)
                    im(idx2) = im(idx2) + im(idx4)

                    u1_re = tc_re - td_im
                    u1_im = tc_re + td_im
                    v1_re = tc_im - td_re
                    v1_im = tc_im + td_re

                    re(idx3) = u1_im * tw1_re - v1_re * tw1_im
                    im(idx3) = -v1_re * tw1_re - u1_im * tw1_im
                    re(idx4) = u1_re * tw3_re + v1_im * tw3_im
                    im(idx4) = v1_im * tw3_re - u1_re * tw3_im
                end do
            end do
        end do

    end subroutine execute_do_classique_dif

end submodule NAFPack_butterfly_split_radix
