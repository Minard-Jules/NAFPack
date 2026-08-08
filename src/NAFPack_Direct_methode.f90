!> Module for direct methods in NAFPack
module NAFPack_Direct_method

    use NAFPack_kinds, only: dp
    use NAFPack_constant, only: TOL_PIVOT_dp

    use NAFPack_Direct_types, only: MethodTypeDirect, METHOD_DIRECT_NONE, &
                                    METHOD_CHOLESKY, METHOD_LDL_Cholesky, &
                                    METHOD_FADDEEV_LEVERRIER, &
                                    METHOD_Gauss, METHOD_Gauss_JORDAN, &
                                    METHOD_LU, METHOD_LDU, &
                                    METHOD_QR, METHOD_TDMA, &
                                    DirectMethodRequirements, MethodQR, &
                                    QR_HOUSEHOLDER, QR_GIVENS, &
                                    QR_GRAM_SCHMIDT, QR_GRAM_SCHMIDT_Modified

    use NAFPack_matrix_decomposition, only: pivot_partial, pivot_total, &
                                            backward, forward, &
                                            LU_decomposition, LDU_decomposition, &
                                            Cholesky_decomposition, LDL_Cholesky_decomposition, &
                                            QR_Gram_Schmidt_Classical_decomposition, &
                                            QR_Gram_Schmidt_Modified_decomposition, &
                                            QR_Givens_decomposition, QR_Householder_decomposition

    use NAFPack_matrix_properties, only: is_non_zero_diagonal, is_SPD, is_square_matrix, &
                                         is_symmetric, is_tridiagonal

    use NAFPack_matrix_tools, only: Faddeev_Leverrier

    use NAFPack_matricielle, only: Identity_n

    implicit none(type, external)

    private

    public :: DirectMethod
    public :: METHOD_Gauss, METHOD_Gauss_JORDAN
    public :: METHOD_LU, METHOD_LDU
    public :: METHOD_CHOLESKY, METHOD_LDL_Cholesky
    public :: METHOD_QR
    public :: METHOD_TDMA
    public :: METHOD_FADDEEV_LEVERRIER
    public :: QR_HOUSEHOLDER, QR_GIVENS, QR_GRAM_SCHMIDT, QR_GRAM_SCHMIDT_Modified

    type :: DirectMethod
        private
        type(MethodTypeDirect) :: method_type = METHOD_DIRECT_NONE
        type(MethodQR) :: qr_method = QR_GRAM_SCHMIDT
        logical :: use_partial_pivot = .false.
        logical :: use_total_pivot = .false.
        type(DirectMethodRequirements) :: requirements
        procedure(solve_interface_Direct), pass(this), pointer :: solve_method => null()

    contains

        procedure :: set_method => set_method
        procedure :: set_qr_method => set_qr_method
        procedure :: solve => DirectMethod_solve
        procedure :: test_matrix => test_matrix

    end type DirectMethod

    abstract interface
        function solve_interface_Direct(this, A, b) result(x)
            import :: dp
            import :: DirectMethod
            implicit none(type, external)
            class(DirectMethod), intent(in) :: this
            real(dp), dimension(:, :), intent(in) :: A
            real(dp), dimension(:), intent(in) :: b
            real(dp), dimension(size(A, 1)) :: x
        end function solve_interface_Direct
    end interface

contains

    subroutine set_method(this, method, set_pivot_partial, set_pivot_total)
        class(DirectMethod), intent(inout) :: this
        type(MethodTypeDirect), intent(in) :: method
        logical, optional, intent(in) :: set_pivot_partial, set_pivot_total

        this%use_total_pivot = .false.
        this%use_partial_pivot = .false.
        this%requirements = DirectMethodRequirements()

        select case (method%id)
        case (METHOD_Gauss%id)
            this%solve_method => solve_Gauss
            this%method_type = METHOD_Gauss
            this%requirements%needs_square = .true.
        case (METHOD_Gauss_JORDAN%id)
            this%solve_method => solve_GaussJordan
            this%method_type = METHOD_Gauss_JORDAN
            this%requirements%needs_square = .true.
        case (METHOD_LU%id)
            this%solve_method => solve_LU
            this%method_type = METHOD_LU
            this%requirements%needs_square = .true.
        case (METHOD_LDU%id)
            this%solve_method => solve_LDU
            this%method_type = METHOD_LDU
            this%requirements%needs_square = .true.
            this%requirements%needs_non_zero_diag = .true.
        case (METHOD_CHOLESKY%id)
            this%solve_method => solve_Cholesky
            this%method_type = METHOD_CHOLESKY
            this%requirements%needs_square = .true.
            this%requirements%needs_SPD = .true.
        case (METHOD_LDL_Cholesky%id)
            this%solve_method => solve_LDL_Cholesky
            this%method_type = METHOD_LDL_Cholesky
            this%requirements%needs_square = .true.
            this%requirements%needs_symmetric = .true.
        case (METHOD_QR%id)
            this%solve_method => solve_QR
            this%method_type = METHOD_QR
            this%requirements%needs_square = .true.
        case (METHOD_TDMA%id)
            this%solve_method => solve_TDMA
            this%method_type = METHOD_TDMA
            this%requirements%needs_square = .true.
            this%requirements%needs_tridiagonal = .true.
            this%requirements%needs_non_zero_diag = .true.
        case (METHOD_FADDEEV_LEVERRIER%id)
            this%solve_method => solve_Faddeev_Leverrier
            this%method_type = METHOD_FADDEEV_LEVERRIER
            this%requirements%needs_square = .true.
        case DEFAULT
            stop "ERROR :: Unknown method direct"
        end select

        if (present(set_pivot_partial)) then
            if (set_pivot_partial) this%use_partial_pivot = .true.
        else if (present(set_pivot_total)) then
            if (set_pivot_total) this%use_total_pivot = .true.
        end if

    end subroutine set_method

    subroutine set_qr_method(this, qr_method)
        class(DirectMethod), intent(inout) :: this
        type(MethodQR), intent(in) :: qr_method

        this%qr_method = qr_method

    end subroutine set_qr_method

    subroutine test_matrix(this, A, strict_mode)
        class(DirectMethod), intent(inout) :: this
        real(dp), dimension(:, :), intent(in) :: A
        logical, optional, intent(in) :: strict_mode
        logical :: strict

        strict = .false.
        if (present(strict_mode)) strict = strict_mode

        if (this%requirements%needs_square) then
            print*,"Checking if the matrix is square..."
            if (.not. is_square_matrix(A)) then
                if (strict) then
                    print*,"ERROR :: ", trim(this%method_type%name), &
                        " method requires a square matrix."
                    stop
                else
                    print*,"WARNING :: ", trim(this%method_type%name), &
                        " method requires a square matrix."
                end if
            end if
        end if

        if (this%requirements%needs_SPD) then
            print*,"Checking if the matrix is symmetric positive definite (SPD)..."
            if (.not. is_SPD(A)) then
                if (strict) then
                    print*,"ERROR :: ", trim(this%method_type%name), &
                        " method requires a symmetric positive definite matrix."
                    stop
                else
                    print*,"WARNING :: ", trim(this%method_type%name), &
                        " method requires a symmetric positive definite matrix."
                end if
            end if
        end if

        if (this%requirements%needs_non_zero_diag) then
            print*,"Checking if the matrix has a non-zero diagonal..."
            if (.not. is_non_zero_diagonal(A)) then
                if (strict) then
                    print*,"ERROR :: ", trim(this%method_type%name), &
                        " method requires a non-zero diagonal matrix."
                    stop
                else
                    print*,"WARNING :: ", trim(this%method_type%name), &
                        " method requires a non-zero diagonal matrix."
                end if
            end if
        end if

        if (this%requirements%needs_tridiagonal) then
            print*,"Checking if the matrix is tridiagonal..."
            if (.not. is_tridiagonal(A)) then
                if (strict) then
                    print*,"ERROR :: ", trim(this%method_type%name), &
                        " method requires a tridiagonal matrix."
                    stop
                else
                    print*,"WARNING :: ", trim(this%method_type%name), &
                        " method requires a tridiagonal matrix."
                end if
            end if
        end if

        if (this%requirements%needs_symmetric) then
            print*,"Checking if the matrix is symmetric..."
            if (.not. is_symmetric(A)) then
                if (strict) then
                    print*,"ERROR :: ", trim(this%method_type%name), &
                        " method requires a symmetric matrix."
                    stop
                else
                    print*,"WARNING :: ", trim(this%method_type%name), &
                        " method requires a symmetric matrix."
                end if
            end if
        end if

    end subroutine test_matrix

    function DirectMethod_solve(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x

        if (.not. associated(this%solve_method)) then
            stop "ERROR :: No solution method has been set. Call set_method first."
        end if

        x = this%solve_method(A, b)

    end function DirectMethod_solve

    function solve_Gauss(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1), size(A, 2)) :: A_tmp
        real(dp), dimension(:, :), allocatable :: P
        real(dp), dimension(:, :), allocatable :: Q
        real(dp), dimension(size(b)) :: b_tmp
        integer :: i, k, N, M, allocate_status
        real(dp) :: pivot, multiplier

        N = size(A, 1)
        M = size(A, 2)

        if (this%use_partial_pivot) then
            allocate (P(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate P"

            call pivot_partial(A, P)

            A_tmp = matmul(P, A)
            b_tmp = matmul(P, b)
        else if (this%use_total_pivot) then
            allocate (P(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate P"
            P = Identity_n(N)
            allocate (Q(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate Q"
            Q = Identity_n(N)

            call pivot_total(A, P, Q)

            A_tmp = matmul(P, A)
            A_tmp = matmul(A, Q)

            b_tmp = matmul(P, b)
        else
            A_tmp = A
            b_tmp = b
        end if

        do k = 1, N - 1
            pivot = A_tmp(k, k)
            if (abs(pivot) < TOL_PIVOT_dp) stop "ERROR :: Near-zero pivot – matrix may be singular"

            do i = k + 1, N
                multiplier = A_tmp(i, k) / pivot
                A_tmp(i, k) = 0

                ! Vectorized operation
                A_tmp(i, k + 1:N) = A_tmp(i, k + 1:N) - multiplier * A_tmp(k, k + 1:N)
                b_tmp(i) = b_tmp(i) - multiplier * b_tmp(k)
            end do
        end do

        x = backward(A_tmp, b_tmp)
        if (this%use_total_pivot) x = matmul(Q, x)

    end function solve_Gauss

    function solve_GaussJordan(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1), size(A, 2)) :: A_tmp
        real(dp), dimension(:, :), allocatable :: P
        real(dp), dimension(:, :), allocatable :: Q
        real(dp), dimension(size(b)) :: b_tmp
        integer :: i, k, N, M, allocate_status
        real(dp) :: pivot, factor

        N = size(A_tmp, 1)
        M = size(A_tmp, 2)

        if (this%use_partial_pivot) then
            allocate (P(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate P"

            call pivot_partial(A, P)

            A_tmp = matmul(P, A)
            b_tmp = matmul(P, b)
        else if (this%use_total_pivot) then
            allocate (P(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate P"
            P = Identity_n(N)
            allocate (Q(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate Q"
            Q = Identity_n(N)

            call pivot_total(A, P, Q)

            A_tmp = matmul(P, A)
            A_tmp = matmul(A, Q)

            b_tmp = matmul(P, b)
        else
            A_tmp = A
            b_tmp = b
        end if

        do k = 1, N
            pivot = A_tmp(k, k)
            if (abs(pivot) < TOL_PIVOT_dp) stop "ERROR :: Near-zero pivot – matrix may be singular"

            ! Normalisation du pivot
            A_tmp(k, :) = A_tmp(k, :) / pivot
            b_tmp(k) = b_tmp(k) / pivot

            ! Élimination dans toutes les autres lignes
            do i = 1, N
                if (i /= k) then
                    factor = A_tmp(i, k)
                    A_tmp(i, :) = A_tmp(i, :) - factor * A_tmp(k, :)
                    b_tmp(i) = b_tmp(i) - factor * b_tmp(k)
                end if
            end do
        end do

        x = b_tmp
        if (this%use_total_pivot) x = matmul(Q, x)

    end function solve_GaussJordan

    function solve_LU(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1), size(A, 1)) :: L, U
        real(dp), dimension(size(A, 1), size(A, 2)) :: A_tmp
        real(dp), dimension(:, :), allocatable :: P
        real(dp), dimension(:, :), allocatable :: Q
        real(dp), dimension(size(b)) :: b_tmp
        integer :: N, M, allocate_status

        N = size(A, 1)
        M = size(A_tmp, 2)

        if (this%use_partial_pivot) then
            allocate (P(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate P"

            call pivot_partial(A, P)

            A_tmp = matmul(P, A)
            b_tmp = matmul(P, b)
        else if (this%use_total_pivot) then
            allocate (P(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate P"
            P = Identity_n(N)
            allocate (Q(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate Q"
            Q = Identity_n(N)

            call pivot_total(A, P, Q)

            A_tmp = matmul(P, A)
            A_tmp = matmul(A, Q)

            b_tmp = matmul(P, b)
        else
            A_tmp = A
            b_tmp = b
        end if

        call LU_decomposition(A_tmp, L, U)

        x = forward(L, b_tmp)

        x = backward(U, x)

        if (this%use_total_pivot) x = matmul(Q, x)

    end function solve_LU

    function solve_LDU(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1), size(A, 1)) :: L, D, U
        real(dp), dimension(size(A, 1), size(A, 2)) :: A_tmp
        real(dp), dimension(:, :), allocatable :: P
        real(dp), dimension(:, :), allocatable :: Q
        real(dp), dimension(size(b)) :: b_tmp
        integer :: N, M, allocate_status

        N = size(A, 1)
        M = size(A, 2)

        if (this%use_partial_pivot) then
            allocate (P(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate P"

            call pivot_partial(A, P)

            A_tmp = matmul(P, A)
            b_tmp = matmul(P, b)
        else if (this%use_total_pivot) then
            allocate (P(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate P"
            P = Identity_n(N)
            allocate (Q(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate Q"
            Q = Identity_n(N)

            call pivot_total(A, P, Q)

            A_tmp = matmul(P, A)
            A_tmp = matmul(A, Q)

            b_tmp = matmul(P, b)
        else
            A_tmp = A
            b_tmp = b
        end if

        call LDU_decomposition(A_tmp, L, D, U)

        x = forward(L, b_tmp)

        x = forward(D, x)

        x = backward(U, x)
        if (this%use_total_pivot) x = matmul(Q, x)

    end function solve_LDU

    function solve_Cholesky(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1), size(A, 1)) :: L

        call Cholesky_decomposition(A, L)

        x = forward(L, b)

        x = backward(transpose(L), x)

    end function solve_Cholesky

    function solve_LDL_Cholesky(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1), size(A, 1)) :: L, D

        call LDL_Cholesky_decomposition(A, L, D)

        x = forward(L, b)

        x = forward(D, x)

        x = backward(transpose(L), x)

    end function solve_LDL_Cholesky

    function solve_QR(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1), size(A, 2)) :: Q, R

        select case (this%qr_method%id)
        case (QR_HOUSEHOLDER%id)
            call QR_Householder_decomposition(A, Q, R)
        case (QR_GIVENS%id)
            call QR_Givens_decomposition(A, Q, R)
        case (QR_GRAM_SCHMIDT%id)
            call QR_Gram_Schmidt_Classical_decomposition(A, Q, R)
        case (QR_GRAM_SCHMIDT_Modified%id)
            call QR_Gram_Schmidt_Modified_decomposition(A, Q, R)
        case DEFAULT
            stop "ERROR :: Unknown QR method"
        end select

        x = backward(R, matmul(transpose(Q), b))

    end function solve_QR

    function solve_TDMA(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: alpha, beta
        real(dp) :: denom
        integer :: n, i

        N = size(A, 1)

        alpha = 0.0_dp
        beta = 0.0_dp

        alpha(1) = A(1, 2) / A(1, 1)
        beta(1) = b(1) / A(1, 1)
        do i = 2, N
            denom = A(i, i) - A(i, i - 1) * alpha(i - 1)
            if (i < N) alpha(i) = A(i, i + 1) / denom
            beta(i) = (b(i) - A(i, i - 1) * beta(i - 1)) / denom
        end do

        x(n) = beta(n)
        do i = n - 1, 1, -1
            x(i) = beta(i) - alpha(i) * x(i + 1)
        end do

    end function solve_TDMA

    function solve_Faddeev_Leverrier(this, A, b) result(x)
        class(DirectMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1), size(A, 2)) :: Ainv
        real(dp), dimension(size(A, 1) + 1) :: c
        logical :: success

        call Faddeev_Leverrier(A, c, Ainv=Ainv, success=success, check=.false.)
        if (.not. success) then
            print*,"WARNING :: Faddeev-Leverrier method failed, using LU decomposition instead"
            x = solve_LU(this, A, b)
        else
            x = matmul(Ainv, b)
        end if

    end function solve_Faddeev_Leverrier

end module NAFPack_Direct_method
