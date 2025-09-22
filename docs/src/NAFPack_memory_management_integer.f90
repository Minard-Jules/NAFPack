submodule(NAFPack_memory_management) NAFPack_memory_management_integer

    implicit none(type, external)

contains

    module subroutine realloc_vec_integer_i8_1D(vec, new_size)
        integer(i8), dimension(:), allocatable, intent(inout) :: vec
        integer(isp), intent(in) :: new_size
        integer(i8), dimension(:), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_i8_1D

    module subroutine realloc_vec_integer_i16_1D(vec, new_size)
        integer(i16), dimension(:), allocatable, intent(inout) :: vec
        integer(isp), intent(in) :: new_size
        integer(i16), dimension(:), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_i16_1D

    module subroutine realloc_vec_integer_isp_1D(vec, new_size)
        integer(isp), dimension(:), allocatable, intent(inout) :: vec
        integer(isp), intent(in) :: new_size
        integer(isp), dimension(:), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_isp_1D

    module subroutine realloc_vec_integer_idp_1D(vec, new_size)
        integer(idp), dimension(:), allocatable, intent(inout) :: vec
        integer(isp), intent(in) :: new_size
        integer(idp), dimension(:), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_idp_1D

    module subroutine realloc_vec_integer_i8_2D(vec, new_size)
        integer(i8), dimension(:, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(2), intent(in) :: new_size
        integer(i8), dimension(:, :), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_i8_2D

    module subroutine realloc_vec_integer_i16_2D(vec, new_size)
        integer(i16), dimension(:, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(2), intent(in) :: new_size
        integer(i16), dimension(:, :), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_i16_2D

    module subroutine realloc_vec_integer_isp_2D(vec, new_size)
        integer(isp), dimension(:, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(2), intent(in) :: new_size
        integer(isp), dimension(:, :), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_isp_2D

    module subroutine realloc_vec_integer_idp_2D(vec, new_size)
        integer(idp), dimension(:, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(2), intent(in) :: new_size
        integer(idp), dimension(:, :), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_idp_2D

    module subroutine realloc_vec_integer_i8_3D(vec, new_size)
        integer(i8), dimension(:, :, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(3), intent(in) :: new_size
        integer(i8), dimension(:, :, :), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_i8_3D

    module subroutine realloc_vec_integer_i16_3D(vec, new_size)
        integer(i16), dimension(:, :, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(3), intent(in) :: new_size
        integer(i16), dimension(:, :, :), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_i16_3D

    module subroutine realloc_vec_integer_isp_3D(vec, new_size)
        integer(isp), dimension(:, :, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(3), intent(in) :: new_size
        integer(isp), dimension(:, :, :), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_isp_3D

    module subroutine realloc_vec_integer_idp_3D(vec, new_size)
        integer(idp), dimension(:, :, :), allocatable, intent(inout) :: vec
        integer(isp), dimension(3), intent(in) :: new_size
        integer(idp), dimension(:, :, :), allocatable :: vec_tmp
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
    end subroutine realloc_vec_integer_idp_3D

end submodule NAFPack_memory_management_integer
