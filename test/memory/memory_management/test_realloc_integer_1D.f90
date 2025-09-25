submodule(test_realloc) test_realloc_integer_1D

    implicit none(type, external)

contains

    module subroutine test_realloc_integer_i8_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        integer(i8), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(int(i, i8), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == int(i, i8))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_integer_i8_1D_grow

    module subroutine test_realloc_integer_i8_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        integer(i8), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(int(i, i8), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == int(i, i8))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_integer_i8_1D_shrink

    module subroutine test_realloc_integer_i16_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        integer(i16), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(int(i, i16), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == int(i, i16))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_integer_i16_1D_grow

    module subroutine test_realloc_integer_i16_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        integer(i16), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(int(i, i16), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == int(i, i16))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_integer_i16_1D_shrink

    module subroutine test_realloc_integer_isp_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        integer(isp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(int(i, isp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == int(i, isp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_integer_isp_1D_grow

    module subroutine test_realloc_integer_isp_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        integer(isp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(int(i, isp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == int(i, isp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_integer_isp_1D_shrink

    module subroutine test_realloc_integer_idp_1D_grow(error)
        type(error_type), allocatable, intent(out) :: error
        integer(idp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 5
        n2 = 10
        allocate (array(n1))
        array = [(int(i, idp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n1
            call check(error, array(i) == int(i, idp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_integer_idp_1D_grow

    module subroutine test_realloc_integer_idp_1D_shrink(error)
        type(error_type), allocatable, intent(out) :: error
        integer(idp), allocatable :: array(:)
        integer(isp) :: i, n1, n2

        n1 = 10
        n2 = 5
        allocate (array(n1))
        array = [(int(i, idp), i=1, n1)]
        call realloc(array, n2)

        call check(error, size(array) == n2)
        if (allocated(error)) return

        do i = 1, n2
            call check(error, array(i) == int(i, idp))
            if (allocated(error)) return
        end do
    end subroutine test_realloc_integer_idp_1D_shrink

end submodule test_realloc_integer_1D
