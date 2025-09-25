submodule(test_realloc) test_realloc_real_1D

    implicit none(type, external)

contains

    module subroutine test_realloc_real_sp_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        real(sp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(real(i, sp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == real(i, sp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_real_sp_1D_grow

    module subroutine test_realloc_real_sp_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        real(sp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(real(i, sp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == real(i, sp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_real_sp_1D_shrink

    module subroutine test_realloc_real_dp_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(real(i, dp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == real(i, dp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_real_dp_1D_grow

    module subroutine test_realloc_real_dp_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        real(dp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(real(i, dp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == real(i, sp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_real_dp_1D_shrink

    module subroutine test_realloc_real_qp_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        real(qp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(real(i, qp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == real(i, qp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_real_qp_1D_grow

    module subroutine test_realloc_real_qp_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        real(qp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(real(i, qp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == real(i, sp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_real_qp_1D_shrink

end submodule test_realloc_real_1D
