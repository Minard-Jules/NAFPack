module NAFPack_butterfly

    use NAFPack_kinds, only: isp, sp, dp, qp
    use NAFPack_plan_type, only: StageView
    use NAFPack_loop_method_type, only: LoopMethod
    use NAFPack_loop_method, only: check_loop_method
    use NAFPack_plan_type, only: DecimationMethod, DIT, DIF

    implicit none(type, external)

    private
    public :: execute_butterfly_dft, execute_butterfly_mixed_radix_recursive
    public :: execute_butterfly_radix2_recursive, execute_butterfly_radix2
    public :: execute_butterfly_radix3_recursive, execute_butterfly_radix3
    public :: execute_butterfly_radix4_recursive, execute_butterfly_radix4
    public :: execute_butterfly_radix5_recursive, execute_butterfly_radix5
    public :: execute_butterfly_split_radix_recursive, execute_butterfly_split_radix

    interface
        pure module subroutine execute_butterfly_dft(stage, N, re, im, decimation_method)
            type(StageView), intent(in) :: stage
            integer(isp), intent(in) :: N
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_dft
    end interface

    interface
        pure module subroutine execute_butterfly_mixed_radix_recursive(stages, N, re, im, stage_idx, decimation_method)
            type(StageView), dimension(:), intent(in) :: stages
            integer(isp), intent(in) :: N, stage_idx
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_mixed_radix_recursive
    end interface

    interface
        pure module subroutine execute_butterfly_radix2_recursive(stages, N, re, im, stage_idx, decimation_method)
            type(StageView), dimension(:), intent(in) :: stages
            integer(isp), intent(in) :: N, stage_idx
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_radix2_recursive

        module subroutine execute_butterfly_radix2(stage, N, re, im, loop_method, decimation_method)
            type(StageView), intent(in) :: stage
            integer(isp), intent(in) :: N
            real(sp), dimension(N), intent(inout) :: re, im
            type(LoopMethod), intent(in) :: loop_method
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_radix2
    end interface

    interface
        pure module subroutine execute_butterfly_radix3_recursive(stages, N, re, im, stage_idx, decimation_method)
            type(StageView), dimension(:), intent(in) :: stages
            integer(isp), intent(in) :: N, stage_idx
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_radix3_recursive

        pure module subroutine execute_butterfly_radix3(stage, N, re, im, decimation_method)
            type(StageView), intent(in) :: stage
            integer(isp), intent(in) :: N
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_radix3
    end interface

    interface
        pure module subroutine execute_butterfly_radix4_recursive(stages, N, re, im, stage_idx, decimation_method)
            type(StageView), dimension(:), intent(in) :: stages
            integer(isp), intent(in) :: N, stage_idx
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_radix4_recursive

        pure module subroutine execute_butterfly_radix4(stage, N, re, im, decimation_method)
            type(StageView), intent(in) :: stage
            integer(isp), intent(in) :: N
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_radix4
    end interface

    interface
        pure module subroutine execute_butterfly_radix5_recursive(stages, N, re, im, stage_idx, decimation_method)
            type(StageView), dimension(:), intent(in) :: stages
            integer(isp), intent(in) :: N, stage_idx
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_radix5_recursive

        pure module subroutine execute_butterfly_radix5(stage, N, re, im, decimation_method)
            type(StageView), intent(in) :: stage
            integer(isp), intent(in) :: N
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_radix5
    end interface

    interface
        pure recursive module subroutine execute_butterfly_split_radix_recursive(stages, N, re, im, stage_idx, decimation_method)
            type(StageView), dimension(:), intent(in) :: stages
            integer(isp), intent(in) :: N, stage_idx
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_split_radix_recursive

        module subroutine execute_butterfly_split_radix(stage, N, re, im, decimation_method)
            type(StageView), intent(in) :: stage
            integer(isp), intent(in) :: N
            real(sp), dimension(N), intent(inout) :: re, im
            type(DecimationMethod), intent(in) :: decimation_method
        end subroutine execute_butterfly_split_radix
    end interface

end module NAFPack_butterfly
