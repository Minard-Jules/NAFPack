submodule(test_realloc) test_realloc_integer_3D

    implicit none(type, external)

contains

    module subroutine test_realloc_integer_i8_3D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        integer(i8), allocatable, dimension(:, :, :) :: array
        integer(isp) :: i, j, k, n1, n2, m1, m2, l1, l2

        n1 = 2
        n2 = 4
        m1 = 3
        m2 = 5
        l1 = 2
        l2 = 3
        allocate (array(n1, m1, l1))
        do k = 1, l1
            do j = 1, m1
                array(:, j, k) = [(int(i + j + k, i8), i=1, n1)]
            end do
        end do
        call realloc(array, [n2, m2, l2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        call check(error, size(array, 3) == l2)
        if (allocated(error)) return

        do k = 1, l1
            do j = 1, m1
                do i = 1, n1
                    call check(error, array(i, j, k) == int(i + j + k, i8))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_integer_i8_3D_grow

    module subroutine test_realloc_integer_i8_3D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        integer(i8), allocatable, dimension(:, :, :) :: array
        integer(isp) :: i, j, k, n1, n2, m1, m2, l1, l2

        n1 = 4
        n2 = 2
        m1 = 5
        m2 = 3
        l1 = 3
        l2 = 2
        allocate (array(n1, m1, l1))
        do k = 1, l1
            do j = 1, m1
                array(:, j, k) = [(int(i + j + k, i8), i=1, n1)]
            end do
        end do
        call realloc(array, [n2, m2, l2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        call check(error, size(array, 3) == l2)
        if (allocated(error)) return

        do k = 1, l2
            do j = 1, m2
                do i = 1, n2
                    call check(error, array(i, j, k) == int(i + j + k, i8))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_integer_i8_3D_shrink

    module subroutine test_realloc_integer_i16_3D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        integer(i16), allocatable, dimension(:, :, :) :: array
        integer(isp) :: i, j, k, n1, n2, m1, m2, l1, l2

        n1 = 2
        n2 = 4
        m1 = 3
        m2 = 5
        l1 = 2
        l2 = 3
        allocate (array(n1, m1, l1))
        do k = 1, l1
            do j = 1, m1
                array(:, j, k) = [(int(i + j + k, i16), i=1, n1)]
            end do
        end do
        call realloc(array, [n2, m2, l2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        call check(error, size(array, 3) == l2)
        if (allocated(error)) return

        do k = 1, l1
            do j = 1, m1
                do i = 1, n1
                    call check(error, array(i, j, k) == int(i + j + k, i16))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_integer_i16_3D_grow

    module subroutine test_realloc_integer_i16_3D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        integer(i16), allocatable, dimension(:, :, :) :: array
        integer(isp) :: i, j, k, n1, n2, m1, m2, l1, l2

        n1 = 4
        n2 = 2
        m1 = 5
        m2 = 3
        l1 = 3
        l2 = 2
        allocate (array(n1, m1, l1))
        do k = 1, l1
            do j = 1, m1
                array(:, j, k) = [(int(i + j + k, i16), i=1, n1)]
            end do
        end do
        call realloc(array, [n2, m2, l2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        call check(error, size(array, 3) == l2)
        if (allocated(error)) return

        do k = 1, l2
            do j = 1, m2
                do i = 1, n2
                    call check(error, array(i, j, k) == int(i + j + k, i16))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_integer_i16_3D_shrink

    module subroutine test_realloc_integer_isp_3D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        integer(isp), allocatable, dimension(:, :, :) :: array
        integer(isp) :: i, j, k, n1, n2, m1, m2, l1, l2

        n1 = 2
        n2 = 4
        m1 = 3
        m2 = 5
        l1 = 2
        l2 = 3
        allocate (array(n1, m1, l1))
        do k = 1, l1
            do j = 1, m1
                array(:, j, k) = [(int(i + j + k, isp), i=1, n1)]
            end do
        end do
        call realloc(array, [n2, m2, l2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        call check(error, size(array, 3) == l2)
        if (allocated(error)) return

        do k = 1, l1
            do j = 1, m1
                do i = 1, n1
                    call check(error, array(i, j, k) == int(i + j + k, isp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_integer_isp_3D_grow

    module subroutine test_realloc_integer_isp_3D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        integer(isp), allocatable, dimension(:, :, :) :: array
        integer(isp) :: i, j, k, n1, n2, m1, m2, l1, l2

        n1 = 4
        n2 = 2
        m1 = 5
        m2 = 3
        l1 = 3
        l2 = 2
        allocate (array(n1, m1, l1))
        do k = 1, l1
            do j = 1, m1
                array(:, j, k) = [(int(i + j + k, isp), i=1, n1)]
            end do
        end do
        call realloc(array, [n2, m2, l2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        call check(error, size(array, 3) == l2)
        if (allocated(error)) return

        do k = 1, l2
            do j = 1, m2
                do i = 1, n2
                    call check(error, array(i, j, k) == int(i + j + k, isp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_integer_isp_3D_shrink

    module subroutine test_realloc_integer_idp_3D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        integer(idp), allocatable, dimension(:, :, :) :: array
        integer(isp) :: i, j, k, n1, n2, m1, m2, l1, l2

        n1 = 2
        n2 = 4
        m1 = 3
        m2 = 5
        l1 = 2
        l2 = 3
        allocate (array(n1, m1, l1))
        do k = 1, l1
            do j = 1, m1
                array(:, j, k) = [(int(i + j + k, idp), i=1, n1)]
            end do
        end do
        call realloc(array, [n2, m2, l2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        call check(error, size(array, 3) == l2)
        if (allocated(error)) return

        do k = 1, l1
            do j = 1, m1
                do i = 1, n1
                    call check(error, array(i, j, k) == int(i + j + k, idp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_integer_idp_3D_grow

    module subroutine test_realloc_integer_idp_3D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        integer(idp), allocatable, dimension(:, :, :) :: array
        integer(isp) :: i, j, k, n1, n2, m1, m2, l1, l2

        n1 = 4
        n2 = 2
        m1 = 5
        m2 = 3
        l1 = 3
        l2 = 2
        allocate (array(n1, m1, l1))
        do k = 1, l1
            do j = 1, m1
                array(:, j, k) = [(int(i + j + k, idp), i=1, n1)]
            end do
        end do
        call realloc(array, [n2, m2, l2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        call check(error, size(array, 3) == l2)
        if (allocated(error)) return

        do k = 1, l2
            do j = 1, m2
                do i = 1, n2
                    call check(error, array(i, j, k) == int(i + j + k, idp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_integer_idp_3D_shrink

end submodule test_realloc_integer_3D
