submodule(NAFPack_memory_management) NAFPack_memory_management_real

    implicit none(type, external)

contains

    module subroutine realloc_vec_real_sp_1D(vec, new_size)
        real(sp), dimension(:), allocatable, intent(inout) :: vec
        integer(isp), intent(in) :: new_size
        real(sp), dimension(:), allocatable :: vec_tmp
        integer(isp) :: old_size

        if (.not. allocated(vec)) then
            allocate (vec(new_size))
            return
        end if

        old_size = size(vec)
        call move_alloc(vec, vec_tmp)
        allocate (vec(new_size))
        vec(:) = vec_tmp(1:min(old_size, new_size))
        deallocate (vec_tmp)
    end subroutine realloc_vec_real_sp_1D

    module subroutine realloc_vec_real_dp_1D(vec, new_size)
        real(dp), dimension(:), allocatable, intent(inout) :: vec
        integer(isp), intent(in) :: new_size
        real(dp), dimension(:), allocatable :: vec_tmp
        integer(isp) :: old_size

        if (.not. allocated(vec)) then
            allocate (vec(new_size))
            return
        end if

        old_size = size(vec)
        call move_alloc(vec, vec_tmp)
        allocate (vec(new_size))
        vec(:) = vec_tmp(1:min(old_size, new_size))
        deallocate (vec_tmp)
    end subroutine realloc_vec_real_dp_1D

    module subroutine realloc_vec_real_qp_1D(vec, new_size)
        real(qp), dimension(:), allocatable, intent(inout) :: vec
        integer(isp), intent(in) :: new_size
        real(qp), dimension(:), allocatable :: vec_tmp
        integer(isp) :: old_size

        if (.not. allocated(vec)) then
            allocate (vec(new_size))
            return
        end if

        old_size = size(vec)
        call move_alloc(vec, vec_tmp)
        allocate (vec(new_size))
        vec(:) = vec_tmp(1:min(old_size, new_size))
        deallocate (vec_tmp)
    end subroutine realloc_vec_real_qp_1D

    module subroutine realloc_vec_real_sp_2D(vec, new_size)
        real(sp), dimension(:, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(2), intent(in) :: new_size
        real(sp), dimension(:, :), allocatable :: vec_tmp
        integer(isp), dimension(2) :: old_size

        if (.not. allocated(vec)) then
            allocate (vec(new_size(1), new_size(2)))
            return
        end if

        old_size = [size(vec, 1), size(vec, 2)]
        call move_alloc(vec, vec_tmp)
        allocate (vec(new_size(1), new_size(2)))
        vec(:, :) = vec_tmp(1:min(old_size(1), new_size(1)), 1:min(old_size(2), new_size(2)))
        deallocate (vec_tmp)
    end subroutine realloc_vec_real_sp_2D

    module subroutine realloc_vec_real_dp_2D(vec, new_size)
        real(dp), dimension(:, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(2), intent(in) :: new_size
        real(dp), dimension(:, :), allocatable :: vec_tmp
        integer(isp), dimension(2) :: old_size

        if (.not. allocated(vec)) then
            allocate (vec(new_size(1), new_size(2)))
            return
        end if

        old_size = [size(vec, 1), size(vec, 2)]
        call move_alloc(vec, vec_tmp)
        allocate (vec(new_size(1), new_size(2)))
        vec(:, :) = vec_tmp(1:min(old_size(1), new_size(1)), 1:min(old_size(2), new_size(2)))
        deallocate (vec_tmp)
    end subroutine realloc_vec_real_dp_2D

    module subroutine realloc_vec_real_qp_2D(vec, new_size)
        real(qp), dimension(:, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(2), intent(in) :: new_size
        real(qp), dimension(:, :), allocatable :: vec_tmp
        integer(isp), dimension(2) :: old_size

        if (.not. allocated(vec)) then
            allocate (vec(new_size(1), new_size(2)))
            return
        end if

        old_size = [size(vec, 1), size(vec, 2)]
        call move_alloc(vec, vec_tmp)
        allocate (vec(new_size(1), new_size(2)))
        vec(:, :) = vec_tmp(1:min(old_size(1), new_size(1)), 1:min(old_size(2), new_size(2)))
        deallocate (vec_tmp)
    end subroutine realloc_vec_real_qp_2D

    module subroutine realloc_vec_real_sp_3D(vec, new_size)
        real(sp), dimension(:, :, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(3), intent(in) :: new_size
        real(sp), dimension(:, :, :), allocatable :: vec_tmp
        integer(isp), dimension(3) :: old_size

        if (.not. allocated(vec)) then
            allocate (vec(new_size(1), new_size(2), new_size(3)))
            return
        end if

        old_size = [size(vec, 1), size(vec, 2), size(vec, 3)]
        call move_alloc(vec, vec_tmp)
        allocate (vec(new_size(1), new_size(2), new_size(3)))
        vec(:, :, :) = vec_tmp(1:min(old_size(1), new_size(1)), &
                               1:min(old_size(2), new_size(2)), &
                               1:min(old_size(3), new_size(3)))
        deallocate (vec_tmp)
    end subroutine realloc_vec_real_sp_3D

    module subroutine realloc_vec_real_dp_3D(vec, new_size)
        real(dp), dimension(:, :, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(3), intent(in) :: new_size
        real(dp), dimension(:, :, :), allocatable :: vec_tmp
        integer(isp), dimension(3) :: old_size

        if (.not. allocated(vec)) then
            allocate (vec(new_size(1), new_size(2), new_size(3)))
            return
        end if

        old_size = [size(vec, 1), size(vec, 2), size(vec, 3)]
        call move_alloc(vec, vec_tmp)
        allocate (vec(new_size(1), new_size(2), new_size(3)))
        vec(:, :, :) = vec_tmp(1:min(old_size(1), new_size(1)), &
                               1:min(old_size(2), new_size(2)), &
                               1:min(old_size(3), new_size(3)))
        deallocate (vec_tmp)
    end subroutine realloc_vec_real_dp_3D

    module subroutine realloc_vec_real_qp_3D(vec, new_size)
        real(qp), dimension(:, :, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(3), intent(in) :: new_size
        real(qp), dimension(:, :, :), allocatable :: vec_tmp
        integer(isp), dimension(3) :: old_size

        if (.not. allocated(vec)) then
            allocate (vec(new_size(1), new_size(2), new_size(3)))
            return
        end if

        old_size = [size(vec, 1), size(vec, 2), size(vec, 3)]
        call move_alloc(vec, vec_tmp)
        allocate (vec(new_size(1), new_size(2), new_size(3)))
        vec(:, :, :) = vec_tmp(1:min(old_size(1), new_size(1)), &
                               1:min(old_size(2), new_size(2)), &
                               1:min(old_size(3), new_size(3)))
        deallocate (vec_tmp)
    end subroutine realloc_vec_real_qp_3D

end submodule NAFPack_memory_management_real
