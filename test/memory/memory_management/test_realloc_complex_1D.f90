submodule(test_realloc) test_memory_management_complex_1D

    implicit none(type, external)

contains

    module subroutine test_realloc_complex_sp_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        complex(sp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(cmplx(i, i+1, sp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == cmplx(i, i+1, sp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_complex_sp_1D_grow

    module subroutine test_realloc_complex_sp_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        complex(sp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(cmplx(i, i+1, sp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == cmplx(i, i+1, sp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_complex_sp_1D_shrink

    module subroutine test_realloc_complex_dp_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        complex(dp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(cmplx(i, i+1, dp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == cmplx(i, i+1, dp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_complex_dp_1D_grow

    module subroutine test_realloc_complex_dp_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        complex(dp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(cmplx(i, i+1, dp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == cmplx(i, i+1, sp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_complex_dp_1D_shrink

    module subroutine test_realloc_complex_qp_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        complex(qp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(cmplx(i, i+1, qp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == cmplx(i, i+1, qp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_complex_qp_1D_grow

    module subroutine test_realloc_complex_qp_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        complex(qp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(cmplx(i, i+1, qp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == cmplx(i, i+1, sp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_complex_qp_1D_shrink

end submodule test_memory_management_complex_1D
