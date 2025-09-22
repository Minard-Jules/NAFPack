module NAFPack_memory_management

    use NAFPack_kinds, only: sp, dp, qp, i8, i16, isp, idp

    implicit none(type, external)

    private
    public :: realloc

    interface realloc
        module procedure realloc_vec_integer_i8_1D
        module procedure realloc_vec_integer_i16_1D
        module procedure realloc_vec_integer_isp_1D
        module procedure realloc_vec_integer_idp_1D
        module procedure realloc_vec_integer_i8_2D
        module procedure realloc_vec_integer_i16_2D
        module procedure realloc_vec_integer_isp_2D
        module procedure realloc_vec_integer_idp_2D
        module procedure realloc_vec_integer_i8_3D
        module procedure realloc_vec_integer_i16_3D
        module procedure realloc_vec_integer_isp_3D
        module procedure realloc_vec_integer_idp_3D

        module procedure realloc_vec_real_sp_1D
        module procedure realloc_vec_real_dp_1D
        module procedure realloc_vec_real_qp_1D
        module procedure realloc_vec_real_sp_2D
        module procedure realloc_vec_real_dp_2D
        module procedure realloc_vec_real_qp_2D
        module procedure realloc_vec_real_sp_3D
        module procedure realloc_vec_real_dp_3D
        module procedure realloc_vec_real_qp_3D

        module procedure realloc_vec_complex_sp_1D
        module procedure realloc_vec_complex_dp_1D
        module procedure realloc_vec_complex_qp_1D
        module procedure realloc_vec_complex_sp_2D
        module procedure realloc_vec_complex_dp_2D
        module procedure realloc_vec_complex_qp_2D
        module procedure realloc_vec_complex_sp_3D
        module procedure realloc_vec_complex_dp_3D
        module procedure realloc_vec_complex_qp_3D
        
    end interface realloc

    interface
        module subroutine realloc_vec_integer_i8_1D(vec, new_size)
            integer(i8), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_integer_i8_1D

        module subroutine realloc_vec_integer_i16_1D(vec, new_size)
            integer(i16), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_integer_i16_1D

        module subroutine realloc_vec_integer_isp_1D(vec, new_size)
            integer(isp), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_integer_isp_1D

        module subroutine realloc_vec_integer_idp_1D(vec, new_size)
            integer(idp), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_integer_idp_1D
    end interface
    interface
        module subroutine realloc_vec_integer_i8_2D(vec, new_size)
            integer(i8), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_integer_i8_2D

        module subroutine realloc_vec_integer_i16_2D(vec, new_size)
            integer(i16), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_integer_i16_2D

        module subroutine realloc_vec_integer_isp_2D(vec, new_size)
            integer(isp), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_integer_isp_2D

        module subroutine realloc_vec_integer_idp_2D(vec, new_size)
            integer(idp), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_integer_idp_2D
    end interface
    interface
        module subroutine realloc_vec_integer_i8_3D(vec, new_size)
            integer(i8), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_integer_i8_3D

        module subroutine realloc_vec_integer_i16_3D(vec, new_size)
            integer(i16), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_integer_i16_3D

        module subroutine realloc_vec_integer_isp_3D(vec, new_size)
            integer(isp), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_integer_isp_3D

        module subroutine realloc_vec_integer_idp_3D(vec, new_size)
            integer(idp), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_integer_idp_3D
    end interface

    interface
        module subroutine realloc_vec_real_sp_1D(vec, new_size)
            real(sp), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_real_sp_1D

        module subroutine realloc_vec_real_dp_1D(vec, new_size)
            real(dp), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_real_dp_1D

        module subroutine realloc_vec_real_qp_1D(vec, new_size)
            real(qp), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_real_qp_1D
    end interface
    interface
        module subroutine realloc_vec_real_sp_2D(vec, new_size)
            real(sp), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_real_sp_2D

        module subroutine realloc_vec_real_dp_2D(vec, new_size)
            real(dp), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_real_dp_2D

        module subroutine realloc_vec_real_qp_2D(vec, new_size)
            real(qp), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_real_qp_2D
    end interface
    interface
        module subroutine realloc_vec_real_sp_3D(vec, new_size)
            real(sp), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_real_sp_3D

        module subroutine realloc_vec_real_dp_3D(vec, new_size)
            real(dp), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_real_dp_3D

        module subroutine realloc_vec_real_qp_3D(vec, new_size)
            real(qp), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_real_qp_3D
    end interface

    interface
        module subroutine realloc_vec_complex_sp_1D(vec, new_size)
            complex(sp), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_complex_sp_1D

        module subroutine realloc_vec_complex_dp_1D(vec, new_size)
            complex(dp), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_complex_dp_1D

        module subroutine realloc_vec_complex_qp_1D(vec, new_size)
            complex(qp), dimension(:), allocatable, intent(inout) :: vec
            integer(isp), intent(in) :: new_size
        end subroutine realloc_vec_complex_qp_1D
    end interface
    interface
        module subroutine realloc_vec_complex_sp_2D(vec, new_size)
            complex(sp), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_complex_sp_2D

        module subroutine realloc_vec_complex_dp_2D(vec, new_size)
            complex(dp), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_complex_dp_2D

        module subroutine realloc_vec_complex_qp_2D(vec, new_size)
            complex(qp), dimension(:, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(2), intent(in) :: new_size
        end subroutine realloc_vec_complex_qp_2D
    end interface
    interface
        module subroutine realloc_vec_complex_sp_3D(vec, new_size)
            complex(sp), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_complex_sp_3D

        module subroutine realloc_vec_complex_dp_3D(vec, new_size)
            complex(dp), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_complex_dp_3D

        module subroutine realloc_vec_complex_qp_3D(vec, new_size)
            complex(qp), dimension(:, :, :), allocatable, intent(inout) :: vec
            integer(isp), dimension(3), intent(in) :: new_size
        end subroutine realloc_vec_complex_qp_3D
    end interface

end module NAFPack_memory_management
