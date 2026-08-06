submodule(NAFPack_butterfly) NAFPack_butterfly_radix2

    implicit none(type, external)

contains

    pure module subroutine execute_butterfly_radix2_recursive(stages, N, re, im, stage_idx, decimation_method)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        type(DecimationMethod), intent(in) :: decimation_method

        if (stage_idx < 1 .or. stage_idx > size(stages)) &
            error stop "stage_idx out of bounds"

        select case (decimation_method%id)
        case (DIT%id)
            call execute_butterfly_radix2_recursive_dit(stages, N, re, im, stage_idx)
        case (DIF%id)
            call execute_butterfly_radix2_recursive_dif(stages, N, re, im, stage_idx)
        case default
            error stop "Error in execute_butterfly_radix2_recursive: unknown decimation method"
        end select

    end subroutine execute_butterfly_radix2_recursive

    pure recursive subroutine execute_butterfly_radix2_recursive_dit(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/2) :: re1, re2, im1, im2
        real(sp), dimension(N/2) :: tmp1_re, tmp1_im
        integer(isp) :: n2

        if (N == 1) return

        n2 = N / 2

        re1 = re(1:N:2)
        im1 = im(1:N:2)
        re2 = re(2:N:2)
        im2 = im(2:N:2)

        call execute_butterfly_radix2_recursive_dit(stages, n2, re1, im1, stage_idx - 1)
        call execute_butterfly_radix2_recursive_dit(stages, n2, re2, im2, stage_idx - 1)

        associate (tw_re => stages(stage_idx)%twiddles_real(1:n2), &
                   tw_im => stages(stage_idx)%twiddles_imag(1:n2))

            tmp1_re = re2 * tw_re - im2 * tw_im
            tmp1_im = re2 * tw_im + im2 * tw_re
        end associate

        re(1:n2) = re1 + tmp1_re
        im(1:n2) = im1 + tmp1_im
        re(n2 + 1:N) = re1 - tmp1_re
        im(n2 + 1:N) = im1 - tmp1_im

    end subroutine execute_butterfly_radix2_recursive_dit

    pure recursive subroutine execute_butterfly_radix2_recursive_dif(stages, N, re, im, stage_idx)
        type(StageView), dimension(:), intent(in) :: stages
        integer(isp), intent(in) :: N, stage_idx
        real(sp), dimension(N), intent(inout) :: re, im
        real(sp), dimension(N/2) :: re1, re2, im1, im2
        real(sp), dimension(N/2) :: tmp1_re, tmp1_im
        integer(isp) :: n2

        if (N == 1) return

        n2 = N / 2

        re1 = re(1:n2)
        im1 = im(1:n2)
        re2 = re(n2 + 1:N)
        im2 = im(n2 + 1:N)

        tmp1_re = re1 - re2
        tmp1_im = im1 - im2

        re1 = re1 + re2
        im1 = im1 + im2

        associate (tw_re => stages(stage_idx)%twiddles_real(1:n2), &
                   tw_im => stages(stage_idx)%twiddles_imag(1:n2))

            re2 = tmp1_re * tw_re - tmp1_im * tw_im
            im2 = tmp1_re * tw_im + tmp1_im * tw_re
        end associate

        call execute_butterfly_radix2_recursive_dif(stages, n2, re1, im1, stage_idx + 1)
        call execute_butterfly_radix2_recursive_dif(stages, n2, re2, im2, stage_idx + 1)

        re(1:N:2) = re1
        im(1:N:2) = im1
        re(2:N:2) = re2
        im(2:N:2) = im2

    end subroutine execute_butterfly_radix2_recursive_dif

    module subroutine execute_butterfly_radix2(stage, N, re, im, loop_method, decimation_method)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        type(LoopMethod), intent(in) :: loop_method
        type(DecimationMethod), intent(in) :: decimation_method

        if (loop_method%use_do_classic) then
            select case (decimation_method%id)
            case (DIT%id)
                call execute_do_classique_dit(stage, N, re, im)
            case (DIF%id)
                call execute_do_classique_dif(stage, N, re, im)
            case default
                error stop "Error in execute_butterfly_radix2: unknown decimation method"
            end select
        else if (loop_method%vectorization%use_simd_openmp) then
            select case (decimation_method%id)
            case (DIT%id)
                call execute_simd_openmp_dit(stage, N, re, im)
            case (DIF%id)
                call execute_simd_openmp_dif(stage, N, re, im)
            case default
                error stop "Error in execute_butterfly_radix2: unknown decimation method"
            end select
        else if (loop_method%vectorization%use_array_syntax) then
            select case (decimation_method%id)
            case (DIT%id)
                call execute_array_syntax_dit(stage, N, re, im)
            case (DIF%id)
                call execute_array_syntax_dif(stage, N, re, im)
            case default
                error stop "Error in execute_butterfly_radix2: unknown decimation method"
            end select
        else if (loop_method%use_do_concurrent) then
            select case (decimation_method%id)
            case (DIT%id)
                call execute_do_concurrent_dit(stage, N, re, im)
            case (DIF%id)
                call execute_do_concurrent_dif(stage, N, re, im)
            case default
                error stop "Error in execute_butterfly_radix2: unknown decimation method"
            end select
        else if (loop_method%parallel%use_openmp) then
            select case (decimation_method%id)
            case (DIT%id)
                call execute_openmp_dit(stage, N, re, im, loop_method%parallel%num_threads)
            case (DIF%id)
                call execute_openmp_dif(stage, N, re, im, loop_method%parallel%num_threads)
            case default
                error stop "Error in execute_butterfly_radix2: unknown decimation method"
            end select
        else
            error stop "Error in execute_butterfly_radix2: no valid loop method selected"
        end if
    end subroutine execute_butterfly_radix2

    pure subroutine kernel_butterfly_radix2_dit(re_idx1, im_idx1, re_idx2, im_idx2, tw_re, tw_im)
        real(sp), intent(inout) :: re_idx1, im_idx1, re_idx2, im_idx2
        real(sp), intent(in) :: tw_re, tw_im
        real(sp) :: tmp_re, tmp_im

        !Complex multiplication
        tmp_re = tw_re * re_idx2 - tw_im * im_idx2
        tmp_im = tw_re * im_idx2 + tw_im * re_idx2

        !Butterfly computations
        re_idx2 = re_idx1 - tmp_re
        im_idx2 = im_idx1 - tmp_im

        re_idx1 = re_idx1 + tmp_re
        im_idx1 = im_idx1 + tmp_im

    end subroutine kernel_butterfly_radix2_dit

    pure subroutine kernel_butterfly_radix2_dif(re_idx1, im_idx1, re_idx2, im_idx2, tw_re, tw_im)
        real(sp), intent(inout) :: re_idx1, im_idx1, re_idx2, im_idx2
        real(sp), intent(in) :: tw_re, tw_im
        real(sp) :: tmp_re, tmp_im

        ! Butterfly computations
        tmp_re = re_idx1 - re_idx2
        tmp_im = im_idx1 - im_idx2

        re_idx1 = re_idx1 + re_idx2
        im_idx1 = im_idx1 + im_idx2

        ! Complex multiplication
        re_idx2 = tmp_re * tw_re - tmp_im * tw_im
        im_idx2 = tmp_re * tw_im + tmp_im * tw_re

    end subroutine kernel_butterfly_radix2_dif

    pure module subroutine execute_do_classique_dit(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, num_groups, offset, step
        real(sp) :: tw_re, tw_im, tmp1_re, tmp1_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride

                tw_re = stage%twiddles_real(twiddle_idx)
                tw_im = stage%twiddles_imag(twiddle_idx)

                tmp1_re = re(idx2) * tw_re - im(idx2) * tw_im
                tmp1_im = re(idx2) * tw_im + im(idx2) * tw_re

                re(idx2) = re(idx1) - tmp1_re
                im(idx2) = im(idx1) - tmp1_im

                re(idx1) = re(idx1) + tmp1_re
                im(idx1) = im(idx1) + tmp1_im
            end do
        end do

    end subroutine execute_do_classique_dit

    pure module subroutine execute_do_classique_dif(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, num_groups, offset, step
        real(sp) :: tw_re, tw_im, tmp1_re, tmp1_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride

                tw_re = stage%twiddles_real(twiddle_idx)
                tw_im = stage%twiddles_imag(twiddle_idx)

                tmp1_re = re(idx1) - re(idx2)
                tmp1_im = im(idx1) - im(idx2)

                re(idx1) = re(idx1) + re(idx2)
                im(idx1) = im(idx1) + im(idx2)

                re(idx2) = tmp1_re * tw_re - tmp1_im * tw_im
                im(idx2) = tmp1_re * tw_im + tmp1_im * tw_re
            end do
        end do

    end subroutine execute_do_classique_dif

    module subroutine execute_simd_openmp_dit(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, num_groups, offset, step
        real(sp) :: current_tw_re, current_tw_im, tmp_re, tmp_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            !$omp simd private(idx1, idx2, current_tw_re, current_tw_im, tmp_re, tmp_im)
            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride

                current_tw_re = stage%twiddles_real(twiddle_idx)
                current_tw_im = stage%twiddles_imag(twiddle_idx)

                !Complex multiplication
                tmp_re = current_tw_re * re(idx2) - current_tw_im * im(idx2)
                tmp_im = current_tw_re * im(idx2) + current_tw_im * re(idx2)

                !Butterfly computations
                re(idx2) = re(idx1) - tmp_re
                im(idx2) = im(idx1) - tmp_im

                re(idx1) = re(idx1) + tmp_re
                im(idx1) = im(idx1) + tmp_im
            end do
            !$omp end simd
        end do

    end subroutine execute_simd_openmp_dit

    module subroutine execute_simd_openmp_dif(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, num_groups, offset, step
        real(sp) :: current_tw_re, current_tw_im, tmp_re, tmp_im

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            !$omp simd private(idx1, idx2, current_tw_re, current_tw_im, tmp_re, tmp_im)
            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride

                current_tw_re = stage%twiddles_real(twiddle_idx)
                current_tw_im = stage%twiddles_imag(twiddle_idx)

                ! Butterfly computations
                tmp_re = re(idx1) - re(idx2)
                tmp_im = im(idx1) - im(idx2)

                re(idx1) = re(idx1) + re(idx2)
                im(idx1) = im(idx1) + im(idx2)

                ! Complex multiplication
                re(idx2) = tmp_re * current_tw_re - tmp_im * current_tw_im
                im(idx2) = tmp_re * current_tw_im + tmp_im * current_tw_re
            end do
            !$omp end simd
        end do

    end subroutine execute_simd_openmp_dif

    pure module subroutine execute_array_syntax_dit(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp), dimension(stage%stride) :: twiddle_idx, idx1, idx2
        integer(isp) :: group_idx, num_groups, offset, step
        real(sp), dimension(stage%stride) :: tmp_re, tmp_im
        integer(isp) :: i

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        twiddle_idx = [(i, i=1, stage%stride)]

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            idx1 = offset + twiddle_idx
            idx2 = idx1 + stage%stride

            !Complex multiplication
            tmp_re = stage%twiddles_real * re(idx2) - stage%twiddles_imag * im(idx2)
            tmp_im = stage%twiddles_real * im(idx2) + stage%twiddles_imag * re(idx2)

            !Butterfly computations
            re(idx2) = re(idx1) - tmp_re
            im(idx2) = im(idx1) - tmp_im

            re(idx1) = re(idx1) + tmp_re
            im(idx1) = im(idx1) + tmp_im

        end do

    end subroutine execute_array_syntax_dit

    pure module subroutine execute_array_syntax_dif(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp), dimension(stage%stride) :: twiddle_idx, idx1, idx2
        integer(isp) :: group_idx, num_groups, offset, step
        real(sp), dimension(stage%stride) :: tmp_re, tmp_im
        integer(isp) :: i

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        twiddle_idx = [(i, i=1, stage%stride)]

        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            idx1 = offset + twiddle_idx
            idx2 = idx1 + stage%stride

            ! Butterfly computations
            tmp_re = re(idx1) - re(idx2)
            tmp_im = im(idx1) - im(idx2)

            re(idx1) = re(idx1) + re(idx2)
            im(idx1) = im(idx1) + im(idx2)

            ! Complex multiplication
            re(idx2) = tmp_re * stage%twiddles_real - tmp_im * stage%twiddles_imag
            im(idx2) = tmp_re * stage%twiddles_imag + tmp_im * stage%twiddles_real
        end do

    end subroutine execute_array_syntax_dif

    pure module subroutine execute_do_concurrent_dit(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, num_groups, offset, step

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        !TODO Fortran 2023 : local(offset, idx1, idx2) — not yet supported by gfortran
        do concurrent(group_idx=0:num_groups - 1, twiddle_idx=1:stage%stride)
            offset = group_idx * step
            idx1 = offset + twiddle_idx
            idx2 = idx1 + stage%stride
            call kernel_butterfly_radix2_dit( &
                re(idx1), im(idx1), &
                re(idx2), im(idx2), &
                stage%twiddles_real(twiddle_idx), stage%twiddles_imag(twiddle_idx))
        end do

    end subroutine execute_do_concurrent_dit

    pure module subroutine execute_do_concurrent_dif(stage, N, re, im)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, num_groups, offset, step

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        !TODO Fortran 2023 : local(offset, idx1, idx2) — not yet supported by gfortran
        do concurrent(group_idx=0:num_groups - 1, twiddle_idx=1:stage%stride)
            offset = group_idx * step
            idx1 = offset + twiddle_idx
            idx2 = idx1 + stage%stride
            call kernel_butterfly_radix2_dif( &
                re(idx1), im(idx1), &
                re(idx2), im(idx2), &
                stage%twiddles_real(twiddle_idx), stage%twiddles_imag(twiddle_idx))
        end do

    end subroutine execute_do_concurrent_dif

    module subroutine execute_openmp_dit(stage, N, re, im, threads)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp), intent(in) :: threads
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, num_groups, offset, step

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        !$omp parallel do default(none) private(group_idx, twiddle_idx, idx1, idx2, offset) &
        !$omp& shared(re, im, stage, step, num_groups) &
        !$omp& num_threads(threads)
        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride
                call kernel_butterfly_radix2_dit( &
                    re(idx1), im(idx1), &
                    re(idx2), im(idx2), &
                    stage%twiddles_real(twiddle_idx), stage%twiddles_imag(twiddle_idx))
            end do
        end do
        !$omp end parallel do

    end subroutine execute_openmp_dit

    module subroutine execute_openmp_dif(stage, N, re, im, threads)
        type(StageView), intent(in) :: stage
        integer(isp), intent(in) :: N
        real(sp), dimension(N), intent(inout) :: re, im
        integer(isp), intent(in) :: threads
        integer(isp) :: group_idx, twiddle_idx, idx1, idx2, num_groups, offset, step

        num_groups = N / stage%butterfly_size
        step = stage%radix * stage%stride

        !$omp parallel do default(none) private(group_idx, twiddle_idx, idx1, idx2, offset) &
        !$omp& shared(re, im, stage, step,num_groups) &
        !$omp& num_threads(threads)
        do group_idx = 0, num_groups - 1
            offset = group_idx * step

            do twiddle_idx = 1, stage%stride
                idx1 = offset + twiddle_idx
                idx2 = idx1 + stage%stride
                call kernel_butterfly_radix2_dif( &
                    re(idx1), im(idx1), &
                    re(idx2), im(idx2), &
                    stage%twiddles_real(twiddle_idx), stage%twiddles_imag(twiddle_idx))
            end do
        end do
        !$omp end parallel do

    end subroutine execute_openmp_dif

end submodule NAFPack_butterfly_radix2
