submodule(test_realloc) test_realloc_real_2D

    implicit none(type, external)

contains

    module subroutine test_realloc_real_sp_2D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        real(sp), allocatable :: array(:, :)
        integer(isp) :: i, j, n1, n2, m1, m2

        n1 = 3
        n2 = 6
        m1 = 4
        m2 = 5
        allocate (array(n1, m1))
        do j = 1, m1
            array(:, j) = [(real(i + j, sp), i=1, n1)]
        end do
        call realloc(array, [n2, m2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        if (allocated(error)) return

        do j = 1, m1
            do i = 1, n1
                call check(error, array(i, j) == real(i + j, sp))
                if (allocated(error)) return
            end do
        end do
    end subroutine test_realloc_real_sp_2D_grow

    module subroutine test_realloc_real_sp_2D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        real(sp), allocatable :: array(:, :)
        integer(isp) :: i, j, n1, n2, m1, m2

        n1 = 6
        n2 = 3
        m1 = 5
        m2 = 4
        allocate (array(n1, m1))
        do j = 1, m1
            array(:, j) = [(real(i + j, sp), i=1, n1)]
        end do
        call realloc(array, [n2, m2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        if (allocated(error)) return

        do j = 1, m2
            do i = 1, n2
                call check(error, array(i, j) == real(i + j, sp))
                if (allocated(error)) return
            end do
        end do
    end subroutine test_realloc_real_sp_2D_shrink

    module subroutine test_realloc_real_dp_2D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp), allocatable :: array(:, :)
        integer(isp) :: i, j, n1, n2, m1, m2

        n1 = 3
        n2 = 6
        m1 = 4
        m2 = 5
        allocate (array(n1, m1))
        do j = 1, m1
            array(:, j) = [(real(i + j, dp), i=1, n1)]
        end do
        call realloc(array, [n2, m2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        if (allocated(error)) return

        do j = 1, m1
            do i = 1, n1
                call check(error, array(i, j) == real(i + j, dp))
                if (allocated(error)) return
            end do
        end do
    end subroutine test_realloc_real_dp_2D_grow

    module subroutine test_realloc_real_dp_2D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp), allocatable :: array(:, :)
        integer(isp) :: i, j, n1, n2, m1, m2

        n1 = 6
        n2 = 3
        m1 = 5
        m2 = 4
        allocate (array(n1, m1))
        do j = 1, m1
            array(:, j) = [(real(i + j, dp), i=1, n1)]
        end do
        call realloc(array, [n2, m2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        if (allocated(error)) return

        do j = 1, m2
            do i = 1, n2
                call check(error, array(i, j) == real(i + j, dp))
                if (allocated(error)) return
            end do
        end do
    end subroutine test_realloc_real_dp_2D_shrink

    module subroutine test_realloc_real_qp_2D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        real(qp), allocatable :: array(:, :)
        integer(isp) :: i, j, n1, n2, m1, m2

        n1 = 3
        n2 = 6
        m1 = 4
        m2 = 5
        allocate (array(n1, m1))
        do j = 1, m1
            array(:, j) = [(real(i + j, qp), i=1, n1)]
        end do
        call realloc(array, [n2, m2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        if (allocated(error)) return

        do j = 1, m1
            do i = 1, n1
                call check(error, array(i, j) == real(i + j, qp))
                if (allocated(error)) return
            end do
        end do
    end subroutine test_realloc_real_qp_2D_grow

    module subroutine test_realloc_real_qp_2D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        real(qp), allocatable :: array(:, :)
        integer(isp) :: i, j, n1, n2, m1, m2

        n1 = 6
        n2 = 3
        m1 = 5
        m2 = 4
        allocate (array(n1, m1))
        do j = 1, m1
            array(:, j) = [(real(i + j, qp), i=1, n1)]
        end do
        call realloc(array, [n2, m2])

        call check(error, size(array, 1) == n2)
        call check(error, size(array, 2) == m2)
        if (allocated(error)) return

        do j = 1, m2
            do i = 1, n2
                call check(error, array(i, j) == real(i + j, qp))
                if (allocated(error)) return
            end do
        end do
    end subroutine test_realloc_real_qp_2D_shrink

end submodule test_realloc_real_2D
