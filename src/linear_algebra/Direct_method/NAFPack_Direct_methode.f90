!> Module for direct methods in NAFPack
MODULE NAFPack_Direct_method

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    USE NAFPack_matrix_properties
    USE NAFPack_Direct_types
    USE NAFPack_matrix_tools

    IMPLICIT NONE

    PRIVATE
    
    PUBLIC :: DirectMethod
    PUBLIC :: METHOD_Gauss
    PUBLIC :: METHOD_LU, METHOD_LDU
    PUBLIC :: METHOD_CHOLESKY, METHOD_LDL_Cholesky
    PUBLIC :: METHOD_QR
    PUBLIC :: METHOD_TDMA
    PUBLIC :: METHOD_FADDEEV_LEVERRIER
    PUBLIC :: QR_HOUSEHOLDER, QR_GIVENS, QR_GRAM_SCHMIDT,QR_GRAM_SCHMIDT_Modified

    TYPE :: DirectMethod
        PRIVATE
        TYPE(MethodTypeDirect) :: method_type = METHOD_DIRECT_NONE
        TYPE(MethodQR) :: qr_method = QR_GRAM_SCHMIDT
        LOGICAL :: use_partial_pivot = .FALSE.
        LOGICAL :: use_total_pivot = .FALSE.
        PROCEDURE(solve_interface_Direct), PASS(this), POINTER :: solve_method => NULL()

        CONTAINS

        PROCEDURE :: set_method => set_method
        PROCEDURE :: set_qr_method => set_qr_method
        PROCEDURE :: solve => DirectMethod_solve
        PROCEDURE :: test_matrix => test_matrix

    END TYPE DirectMethod

    ABSTRACT INTERFACE
        FUNCTION solve_interface_Direct(this, A, b) RESULT(x)
            IMPORT :: dp
            IMPORT :: DirectMethod
            CLASS(DirectMethod), INTENT(IN) :: this
            REAL(dp), DIMENSION(:,:), INTENT(IN) :: A
            REAL(dp), DIMENSION(:), INTENT(IN)  :: b
            REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        END FUNCTION solve_interface_Direct
    END INTERFACE

    CONTAINS

    SUBROUTINE set_method(this, method, set_pivot_partial, set_pivot_total)
        CLASS(DirectMethod), INTENT(INOUT) :: this
        TYPE(MethodTypeDirect), INTENT(IN) :: method
        LOGICAL, OPTIONAL :: set_pivot_partial, set_pivot_total

        this%use_total_pivot = .FALSE.
        this%use_partial_pivot = .FALSE.

        SELECT CASE (method%value)
        CASE (METHOD_Gauss%value)
            this%solve_method => solve_Gauss
            this%method_type = METHOD_Gauss
        CASE (METHOD_LU%value)
            this%solve_method => solve_LU
            this%method_type = METHOD_LU
        CASE (METHOD_LDU%value)
            this%solve_method => solve_LDU
            this%method_type = METHOD_LDU
        CASE (METHOD_CHOLESKY%value)
            this%solve_method => solve_Cholesky
            this%method_type = METHOD_CHOLESKY
        CASE (METHOD_LDL_Cholesky%value)
            this%solve_method => solve_LDL_Cholesky
            this%method_type = METHOD_LDL_Cholesky
        CASE (METHOD_QR%value)
            ! this%solve_method => solve_QR
            this%method_type = METHOD_QR
        CASE (METHOD_TDMA%value)
            this%solve_method => solve_TDMA
            this%method_type = METHOD_TDMA
        CASE (METHOD_FADDEEV_LEVERRIER%value)
            this%solve_method => solve_Faddeev_Leverrier
            this%method_type = METHOD_FADDEEV_LEVERRIER
        CASE DEFAULT
            STOP "ERROR :: Unknown method direct"
        END SELECT

        IF(PRESENT(set_pivot_partial))THEN
            IF(set_pivot_partial) this%use_partial_pivot = .TRUE.
        ELSE IF(PRESENT(set_pivot_total))THEN
            IF(set_pivot_total) this%use_total_pivot = .TRUE.
        END IF
        
    END SUBROUTINE set_method

    SUBROUTINE set_qr_method(this, qr_method)
        CLASS(DirectMethod), INTENT(INOUT) :: this
        TYPE(MethodQR), INTENT(IN) :: qr_method

        this%qr_method = qr_method

    END SUBROUTINE set_qr_method

    SUBROUTINE test_matrix(this, A)
        CLASS(DirectMethod), INTENT(INOUT) :: this
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: A

        SELECT CASE (this%method_type%value)
        CASE (METHOD_Gauss%value)
            IF(.NOT. is_square_matrix(A)) STOP "ERROR :: Gauss method requires a square matrix."
        CASE (METHOD_LU%value)
            IF(.NOT. is_square_matrix(A)) STOP "ERROR :: LU method requires a square matrix."
        CASE (METHOD_LDU%value)
            IF(.NOT. is_non_zero_diagonal(A)) STOP "ERROR :: LDU method requires a non-zero diagonal matrix."
        CASE (METHOD_CHOLESKY%value)
            IF(.NOT. is_SPD(A)) STOP "ERROR :: Cholesky method requires a symmetric positive definite matrix."
        CASE (METHOD_LDL_Cholesky%value)
            IF(.NOT. is_symmetric(A)) STOP "ERROR :: LDL_Cholesky method requires a symmetric matrix."
        CASE (METHOD_QR%value)
            IF(.NOT. is_square_matrix(A)) STOP "ERROR :: QR method requires a square matrix."
        CASE (METHOD_TDMA%value)
            IF(.NOT. is_tridiagonal(A)) STOP "ERROR :: TDMA method requires a tridiagonal matrix."
            IF(.NOT. is_non_zero_diagonal(A, .TRUE.)) STOP "ERROR :: TDMA method requires a non-zero diagonal matrix."
        CASE (METHOD_FADDEEV_LEVERRIER%value)
            IF(.NOT. is_square_matrix(A)) STOP "ERROR :: Faddeev_Leverrier method requires a square matrix."
        CASE DEFAULT
            STOP "ERROR :: Unknown method direct for testing matrix properties"
        END SELECT

    END SUBROUTINE test_matrix

    FUNCTION DirectMethod_solve(this, A, b) RESULT(x)
        CLASS(DirectMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        
        IF (.NOT. ASSOCIATED(this%solve_method)) THEN
            STOP "ERROR :: No solution method has been set. Call set_method first."
        END IF

        x = this%solve_method(A, b)
        
    END FUNCTION DirectMethod_solve

    FUNCTION solve_Gauss(this, A, b) RESULT(x)
        CLASS(DirectMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: A_tmp
        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Q_final
        REAL(dp), DIMENSION(SIZE(b)) :: b_tmp
        INTEGER :: i, k, n, allocate_status
        REAL(dp) :: pivot, m

        A_tmp = A
        b_tmp = b

        N=SIZE(A_tmp, 1)

        IF(this%use_total_pivot)THEN
            ALLOCATE(Q_final(N,N), STAT=allocate_status)
            IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate Q_final"
            Q_final = Identity_n(N)
        END IF

        IF(this%use_partial_pivot) CALL pivot_partial(A_tmp, b_tmp)
        IF(this%use_total_pivot) CALL pivot_total(A_tmp, b_tmp, Q_final)

        DO k = 1, N-1
            pivot = A_tmp(k, k)
            IF (ABS(pivot) < epsi) STOP "ERROR :: Near-zero pivot – matrix may be singular"

            DO i = k+1, N
                m = A_tmp(i, k) / pivot
                A_tmp(i, k) = 0 

                ! Vectorized operation
                A_tmp(i, k+1:N) = A_tmp(i, k+1:N) - m * A_tmp(k, k+1:N)
                b_tmp(i) = b_tmp(i) - m * b_tmp(k)
            END DO
        END DO

        x = backward(A_tmp, b_tmp)
        IF(this%use_total_pivot) x = MATMUL(Q_final, x)

    END FUNCTION solve_Gauss

    FUNCTION solve_LU(this, A, b) RESULT(x)
        CLASS(DirectMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: L, U
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: A_tmp
        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Q_final
        REAL(dp), DIMENSION(SIZE(b)) :: b_tmp
        INTEGER :: N, allocate_status

        A_tmp = A
        b_tmp = b

        N = SIZE(A, 1)

        IF(this%use_total_pivot)THEN
            ALLOCATE(Q_final(N,N), STAT=allocate_status)
            IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate Q_final"
            Q_final = Identity_n(N)
        END IF

        IF(this%use_partial_pivot) CALL pivot_partial(A_tmp, b_tmp)
        IF(this%use_total_pivot) CALL pivot_total(A_tmp, b_tmp, Q_final)

        CALL LU_decomposition(A_tmp, L, U)

        x = forward(L, b_tmp)

        x = backward(U, x)
        IF(this%use_total_pivot) x = MATMUL(Q_final, x)

    END FUNCTION solve_LU

    FUNCTION solve_LDU(this, A, b) RESULT(x)
        CLASS(DirectMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: L, D, U
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: A_tmp
        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Q_final
        REAL(dp), DIMENSION(SIZE(b)) :: b_tmp
        INTEGER :: N, allocate_status

        A_tmp = A
        b_tmp = b

        N = SIZE(A, 1)

        IF(this%use_total_pivot)THEN
            ALLOCATE(Q_final(N,N), STAT=allocate_status)
            IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate Q_final"
            Q_final = Identity_n(N)
        END IF

        IF(this%use_partial_pivot) CALL pivot_partial(A_tmp, b_tmp)
        IF(this%use_total_pivot) CALL pivot_total(A_tmp, b_tmp, Q_final)

        CALL LDU_decomposition(A_tmp, L, D, U)

        x = forward(L, b_tmp)

        x = forward(D, x)

        x = backward(U, x)
        IF(this%use_total_pivot) x = MATMUL(Q_final, x)

    END FUNCTION solve_LDU

    FUNCTION solve_Cholesky(this, A, b) RESULT(x)
        CLASS(DirectMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: L

        CALL Cholesky_decomposition(A, L)
          
        x = forward(L, b)

        x = backward(TRANSPOSE(L), x)

    END FUNCTION solve_Cholesky

    FUNCTION solve_LDL_Cholesky(this, A, b) RESULT(x)
        CLASS(DirectMethod), INTENT(IN) :: this
        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        REAL(dp),DIMENSION(SIZE(A,1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: L, D

        CALL LDL_Cholesky_decomposition(A, L, D)

        x = forward(L, b)

        x = forward(D, x)

        x = backward(TRANSPOSE(L), x)
    
    END FUNCTION solve_LDL_Cholesky

    FUNCTION solve_QR(this, A, b) RESULT(x)
        CLASS(DirectMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)) :: Q, R

        SELECT CASE (this%qr_method%value)
        CASE (QR_HOUSEHOLDER%value)
            CALL QR_Householder_decomposition(A, Q, R)
        CASE (QR_GIVENS%value)
            CALL QR_Givens_decomposition(A, Q, R)
        CASE (QR_GRAM_SCHMIDT%value)
            CALL QR_Gram_Schmidt_Classical_decomposition(A, Q, R)
        CASE (QR_GRAM_SCHMIDT_Modified%value)
            CALL QR_Gram_Schmidt_Modified_decomposition(A, Q, R)
        CASE DEFAULT
            STOP "ERROR :: Unknown QR method"
        END SELECT

        x = backward(R, MATMUL(TRANSPOSE(Q), b))

    END FUNCTION solve_QR

    FUNCTION solve_TDMA(this, A, b) RESULT(x)
        CLASS(DirectMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A,1)) :: x
        REAL(dp), DIMENSION(SIZE(A,1)) :: alpha, beta
        REAL(dp) :: denom
        INTEGER :: n, i

        N = SIZE(A,1)

        alpha = 0.0_dp
        beta = 0.0_dp

        alpha(1) = A(1,2) / A(1,1)
        beta(1) = b(1) / A(1,1)
        DO i = 2, N
            denom = A(i,i) - A(i,i-1)*alpha(i-1)
            IF (i < N) alpha(i) = A(i,i+1) / denom
            beta(i) = (b(i) - A(i,i-1)*beta(i-1)) / denom
        END DO

        x(n) = beta(n)
        DO i = n-1, 1, -1
            x(i) = beta(i) - alpha(i)*x(i+1)
        END DO

    END FUNCTION solve_TDMA

    FUNCTION solve_Faddeev_Leverrier(this, A, b) RESULT(x)
        CLASS(DirectMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)) :: Ainv
        REAL(dp), DIMENSION(SIZE(A, 1)+1) :: c
        LOGICAL :: success

        CALL Faddeev_Leverrier(A, c, Ainv = Ainv, success = success, check = .FALSE.)
        IF (.NOT. success) THEN
            PRINT*, "WARNING :: Faddeev-Leverrier method failed, using LU decomposition instead"
            x = solve_LU(this, A, b)
        ELSE
            x = MATMUL(Ainv, b)
        END IF

    END FUNCTION solve_Faddeev_Leverrier

    SUBROUTINE pivot_partial(A, b)
        REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: A
        REAL(dp), DIMENSION(:), INTENT(INOUT) :: b
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 1)) :: P
        INTEGER, DIMENSION(1) :: vlmax
        INTEGER :: N, lmax, k

        N=SIZE(A, 1)

        DO k = 1, N-1

            ! Find the maximum absolute value in the column from row k to N
            vlmax = MAXLOC(ABS(A(k:N, k)))
            lmax = vlmax(1) + k - 1
            
            !calculate permutation matrix P
            P = Identity_n(N)
            IF (k /= lmax) THEN
                P = rotation_matrix(P, [k, lmax])
            END IF

            A = MATMUL(P, A)
            b = MATMUL(P, b)

        END DO

    END SUBROUTINE pivot_partial

    SUBROUTINE pivot_total(A, b, Q_final)
        REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: A
        REAL(dp), DIMENSION(:), INTENT(INOUT) :: b
        REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: Q_final
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 1)) :: P, Q
        INTEGER, DIMENSION(2) :: vlmax
        INTEGER :: N, lmax, cmax, k

        N=SIZE(A, 1)

        DO k = 1, N-1
            ! Find max abs element in submatrix
            vlmax = MAXLOC(ABS(A(k:N,k:N)))
            lmax = vlmax(1) + k - 1
            cmax = vlmax(2) + k - 1

            ! permute line if necessary
            P = Identity_n(N)
            IF (lmax /= k) THEN
                P = rotation_matrix(P, [k, lmax])
            END IF

            ! permute column if necessary
            Q = Identity_n(N)
            IF (cmax /= k) THEN
                Q = rotation_matrix(Q, [k, cmax])
            END IF

            Q_final = MATMUL(Q, Q_final)

            ! Apply permutations
            A = MATMUL(P, A)
            A = MATMUL(A, Q)

            b = MATMUL(P, b)
        END DO
        
    END SUBROUTINE pivot_total

END MODULE NAFPack_Direct_method