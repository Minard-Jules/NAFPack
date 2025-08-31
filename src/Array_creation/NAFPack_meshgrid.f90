!> Module for creating a meshgrid from two vectors
!>
!> This module provides a subroutine to create a meshgrid.
module NAFPack_meshgrid

    use NAFPack_kinds, only: dp, sp, qp, i8, i16, isp, idp
    use NAFPack_loop_method, only: LoopMethod, count_true_methods, default_loop_method
    USE OMP_LIB

    implicit none(type, external)

    private
    public :: meshgrid
    public :: INDEXING_XY, INDEXING_IJ
    public :: check_indexing, check_loop_method

    type :: meshgrid_indexing
        integer :: id
        character(len=2) :: name
    end type meshgrid_indexing

    type(meshgrid_indexing), parameter :: INDEXING_XY = meshgrid_indexing(1, "XY"), &
                                          INDEXING_IJ = meshgrid_indexing(2, "IJ")

    !> Make N-dimensional meshgrid from two vectors **x_vector** and **y_vector**
    interface meshgrid
        module procedure meshgrid_sp_2D
        module procedure meshgrid_dp_2D
        module procedure meshgrid_qp_2D

        module procedure meshgrid_sp_3D
        module procedure meshgrid_dp_3D
        module procedure meshgrid_qp_3D

        module procedure meshgrid_i8_2D
        module procedure meshgrid_i16_2D
        module procedure meshgrid_isp_2D
        module procedure meshgrid_idp_2D

        module procedure meshgrid_i8_3D
        module procedure meshgrid_i16_3D
        module procedure meshgrid_isp_3D
        module procedure meshgrid_idp_3D

        module procedure meshgrid_cmplx_sp_2D
        module procedure meshgrid_cmplx_dp_2D
        module procedure meshgrid_cmplx_qp_2D

        module procedure meshgrid_cmplx_sp_3D
        module procedure meshgrid_cmplx_dp_3D
        module procedure meshgrid_cmplx_qp_3D
    end interface meshgrid

    interface
        module subroutine meshgrid_sp_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            real(sp), dimension(:), intent(in) :: x_vector, y_vector
            real(sp), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_sp_2D

        module subroutine meshgrid_dp_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            real(dp), dimension(:), intent(in) :: x_vector, y_vector
            real(dp), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_dp_2D

        module subroutine meshgrid_qp_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            real(qp), dimension(:), intent(in) :: x_vector, y_vector
            real(qp), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_qp_2D
    end interface

    interface
        module subroutine meshgrid_sp_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            real(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            real(sp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_sp_3D

        module subroutine meshgrid_dp_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            real(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            real(dp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_dp_3D

        module subroutine meshgrid_qp_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            real(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            real(qp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_qp_3D
    end interface

    interface
        module subroutine meshgrid_i8_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            integer(i8), dimension(:), intent(in) :: x_vector, y_vector
            integer(i8), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_i8_2D

        module subroutine meshgrid_i16_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            integer(i16), dimension(:), intent(in) :: x_vector, y_vector
            integer(i16), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_i16_2D

        module subroutine meshgrid_isp_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            integer(isp), dimension(:), intent(in) :: x_vector, y_vector
            integer(isp), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_isp_2D

        module subroutine meshgrid_idp_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            integer(idp), dimension(:), intent(in) :: x_vector, y_vector
            integer(idp), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_idp_2D
    end interface

    interface
        module subroutine meshgrid_i8_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            integer(i8), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            integer(i8), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_i8_3D

        module subroutine meshgrid_i16_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            integer(i16), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            integer(i16), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_i16_3D

        module subroutine meshgrid_isp_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            integer(isp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            integer(isp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_isp_3D

        module subroutine meshgrid_idp_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            integer(idp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            integer(idp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_idp_3D
    end interface

    interface
        module subroutine meshgrid_cmplx_sp_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            complex(sp), dimension(:), intent(in) :: x_vector, y_vector
            complex(sp), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_cmplx_sp_2D

        module subroutine meshgrid_cmplx_dp_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            complex(dp), dimension(:), intent(in) :: x_vector, y_vector
            complex(dp), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_cmplx_dp_2D

        module subroutine meshgrid_cmplx_qp_2D( &
            x_vector, y_vector, X, Y, indexing, strict_mode, loop_method)
            implicit none(type, external)
            complex(qp), dimension(:), intent(in) :: x_vector, y_vector
            complex(qp), dimension(:, :), allocatable, intent(out) :: X, Y
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_cmplx_qp_2D
    end interface

    interface
        module subroutine meshgrid_cmplx_sp_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            complex(sp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            complex(sp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_cmplx_sp_3D

        module subroutine meshgrid_cmplx_dp_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            complex(dp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            complex(dp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_cmplx_dp_3D

        module subroutine meshgrid_cmplx_qp_3D( &
            x_vector, y_vector, z_vector, X, Y, Z, indexing, strict_mode, loop_method)
            implicit none(type, external)
            complex(qp), dimension(:), intent(in) :: x_vector, y_vector, z_vector
            complex(qp), dimension(:, :, :), allocatable, intent(out) :: X, Y, Z
            type(meshgrid_indexing), optional, intent(in) :: indexing
            logical, optional, intent(in) :: strict_mode
            type(LoopMethod), optional, intent(in) :: loop_method
        end subroutine meshgrid_cmplx_qp_3D
    end interface

contains

    subroutine check_indexing(indexing, strict_mode, use_ij_indexing, use_xy_indexing)
        type(meshgrid_indexing), optional, intent(in) :: indexing
        logical, optional, intent(in) :: strict_mode
        logical, intent(out) :: use_ij_indexing, use_xy_indexing
        logical :: is_strict

        is_strict = .false.
        if (present(strict_mode)) is_strict = strict_mode

        use_ij_indexing = .false.
        use_xy_indexing = .false.
        if (present(indexing)) then
            if (indexing%id == INDEXING_IJ%id) then
                use_ij_indexing = .true.
            else if (indexing%id == INDEXING_XY%id) then
                use_xy_indexing = .true.
            else
                if (is_strict) then
                    error stop "Error: Unknown indexing%id in meshgrid"
                else
                    use_ij_indexing = .true.
                end if
            end if
        else
            use_ij_indexing = .true.
        end if
    end subroutine check_indexing

    function check_loop_method(loop_method) result(loop_method_used)
        type(LoopMethod), intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used
        integer :: nb_of_true

        nb_of_true = count_true_methods(loop_method)
        if(nb_of_true == 0)then
            loop_method_used = default_loop_method
        else if(nb_of_true == 1)then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if
    end function check_loop_method

end module NAFPack_meshgrid
