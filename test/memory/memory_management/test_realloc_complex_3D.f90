submodule(test_realloc) test_memory_management_complex_3D

    implicit none(type, external)

contains

    module subroutine test_realloc_complex_sp_3D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        complex(sp), allocatable, dimension(:, :, :) :: array
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
                array(:, j, k) = [(cmplx(i + j + k, i+1, sp), i=1, n1)]
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
                    call check(error, array(i, j, k) == cmplx(i + j + k, i+1, sp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_complex_sp_3D_grow

    module subroutine test_realloc_complex_sp_3D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        complex(sp), allocatable, dimension(:, :, :) :: array
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
                array(:, j, k) = [(cmplx(i + j + k, i+1, sp), i=1, n1)]
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
                    call check(error, array(i, j, k) == cmplx(i + j + k, i+1, sp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_complex_sp_3D_shrink

    module subroutine test_realloc_complex_dp_3D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        complex(dp), allocatable, dimension(:, :, :) :: array
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
                array(:, j, k) = [(cmplx(i + j + k, i+1, dp), i=1, n1)]
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
                    call check(error, array(i, j, k) == cmplx(i + j + k, i+1, dp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_complex_dp_3D_grow

    module subroutine test_realloc_complex_dp_3D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        complex(dp), allocatable, dimension(:, :, :) :: array
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
                array(:, j, k) = [(cmplx(i + j + k, i+1, dp), i=1, n1)]
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
                    call check(error, array(i, j, k) == cmplx(i + j + k, i+1, dp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_complex_dp_3D_shrink

    module subroutine test_realloc_complex_qp_3D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        complex(qp), allocatable, dimension(:, :, :) :: array
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
                array(:, j, k) = [(cmplx(i + j + k, i+1, qp), i=1, n1)]
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
                    call check(error, array(i, j, k) == cmplx(i + j + k, i+1, qp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_complex_qp_3D_grow

    module subroutine test_realloc_complex_qp_3D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        complex(qp), allocatable, dimension(:, :, :) :: array
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
                array(:, j, k) = [(cmplx(i + j + k, i+1, qp), i=1, n1)]
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
                    call check(error, array(i, j, k) == cmplx(i + j + k, i+1, qp))
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_realloc_complex_qp_3D_shrink

end submodule test_memory_management_complex_3D
