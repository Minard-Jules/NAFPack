!> Module for matrix decomposition methods
!>
!> This module provides subroutines for various matrix decomposition methods including LU, LDU, Cholesky, and QR decompositions.
MODULE NAFPack_matrix_decomposition

    USE NAFPack_constant
    USE NAFPack_matricielle

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: forward, backward
    PUBLIC :: LU_decomposition, LDU_decomposition, ILU_decomposition
    PUBLIC :: Cholesky_decomposition, LDL_Cholesky_decomposition, Incomplete_Cholesky_decomposition
    PUBLIC :: QR_decomposition
    PUBLIC :: QR_Householder_decomposition, QR_Givens_decomposition, &
              QR_Gram_Schmidt_Classical_decomposition, QR_Gram_Schmidt_Modified_decomposition

    CONTAINS

    
    !> forward algorithm, 
    !> solves the system 
    !> \[ L * y = b \]
    !> where **L** is a lower triangular matrix and **b** is a vector
    FUNCTION forward(L, b) RESULT(y)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: L
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(L, 1)) :: y
        INTEGER :: i, N
    
        N = SIZE(L, 1)
    
        y(1) = b(1) / L(1, 1)

        DO i = 2,N
            y(i) = (b(i) - DOT_PRODUCT(L(i,1:i-1), y(1:i-1))) / L(i,i)
        END DO
    
    END FUNCTION forward

    !> backward algorithm, 
    !> solves the system 
    !> \[ U * x = y \]
    !> where **U** is an upper triangular matrix and **y** is a vector
    FUNCTION backward(U, y) RESULT(x)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: U
        REAL(dp), DIMENSION(:), INTENT(IN) :: y
        REAL(dp), DIMENSION(SIZE(U, 1)) :: x
        INTEGER :: i, N

        N = SIZE(U, 1)

        x(N) = y(N) / U(N, N)

        DO i = N-1, 1, -1
            x(i) = (y(i) - DOT_PRODUCT(U(i, i+1:N), x(i+1:N))) / U(i, i)
        END DO

    END FUNCTION backward

    !> LU decomposition of a matrix A
    !> \[ A = LU \]
    !> This subroutine performs LU decomposition of a given matrix **A**, where **L** is a lower triangular matrix and **U** is an upper triangular matrix.
    SUBROUTINE LU_decomposition(A, L, U)

        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)), INTENT(OUT) :: L, U
        INTEGER :: i, j, N

        N = SIZE(A, 1)

        L = 0.d0
        U = 0.d0

        DO j = 1, N
            L(j, j) = 1.d0

            DO i = 1, j
                U(i, j) = A(i, j) - DOT_PRODUCT(L(i, 1:i-1), U(1:i-1, j))
            END DO
    
            DO i = j+1, N
                L(i, j) = (A(i, j) - DOT_PRODUCT(L(i, 1:j-1), U(1:j-1, j))) / U(j, j)
            END DO
        END DO

    END SUBROUTINE LU_decomposition

    !> LDU decomposition of a matrix A
    !> \[ A = LDU \]
    !> This subroutine performs LDU decomposition of a given matrix **A**, where **L** is a lower triangular matrix, **D** is a diagonal matrix, and **U** is an upper triangular matrix.
    SUBROUTINE LDU_decomposition(A, L, D, U)

        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)), INTENT(OUT) :: L, U, D
        INTEGER :: i, j, k, N

        N = SIZE(A, 1)

        L = 0.d0
        D = 0.d0
        U = 0.d0

        DO j = 1, N
            L(j, j) = 1.d0
            U(j, j) = 1.d0
            
            DO i = 1, j-1
                U(i, j) = (A(i, j) - DOT_PRODUCT(L(i, 1:i-1), U(1:i-1, j) * [ (D(k,k), k = 1, i-1) ])) / D(i, i)
            END DO

            i = j
            D(j, j) = A(j, j) - DOT_PRODUCT(L(j, 1:j-1), U(1:j-1, j) * [ (D(k,k), k = 1, j-1) ])

            DO i = j+1, N
                L(i, j) = (A(i, j) - DOT_PRODUCT(L(i, 1:j-1), U(1:j-1, j) * [ (D(k,k), k = 1, j-1) ])) / D(j, j)
            END DO
        END DO

    END SUBROUTINE LDU_decomposition

    !> Incomplete LU decomposition of a matrix A
    !> \[ A \approx LU \]
    !> This subroutine performs incomplete LU decomposition of a given matrix **A**, where **L** is a lower triangular matrix and **U** is an upper triangular matrix.
    SUBROUTINE ILU_decomposition(A, L, U)

        REAL(dp), DIMENSION(:, :), INTENT(IN)  :: A
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 1)), INTENT(OUT) :: L, U
        LOGICAL, DIMENSION(SIZE(A, 1), SIZE(A, 1)) :: S
        INTEGER :: i, j, N

        N = SIZE(A, 1)

        L = Identity_n(N)
        U = 0.d0

        S = A /= 0

        DO i = 1, N
            DO j = 1, i-1
                IF (S(i,j)) L(i,j) = (A(i,j) - DOT_PRODUCT(L(i, 1:j-1), U(1:j-1, j))) / U(j,j)
            END DO
            DO j = i, N
                IF (S(i,j)) U(i,j) = A(i,j) - DOT_PRODUCT(L(i, 1:i-1), U(1:i-1, j))
            END DO
        END DO

    END SUBROUTINE ILU_decomposition

    !> Cholesky decomposition of a matrix A
    !> \[ A = LL^T \]
    !> This subroutine performs Cholesky decomposition of a given symmetric positive definite matrix **A**, where **L** is a lower triangular matrix.
    SUBROUTINE Cholesky_decomposition(A, L)

        REAL(dp),DIMENSION(: ,:), INTENT(IN) :: A
        REAL(dp),DIMENSION(SIZE(A, 1), SIZE(A, 1)), INTENT(OUT) :: L
        INTEGER :: i, j, N

        N = SIZE(A, 1)

        DO j = 1, N
            L(j, j) = SQRT(A(j, j) - DOT_PRODUCT(L(j, 1:j-1), L(j, 1:j-1)))

            DO i = j+1, N
                L(i, j) = (A(i, j) - DOT_PRODUCT(L(i, 1:j-1), L(j, 1:j-1))) / L(j, j)
            END DO
        END DO

    END SUBROUTINE Cholesky_decomposition

    !> Alternative Cholesky decomposition of a matrix A
    !> \[ A = LDL^T \]
    !> This subroutine performs alternative Cholesky decomposition of a given symmetric positive definite matrix **A**, where **L** is a lower triangular matrix and **D** is a diagonal matrix.
    SUBROUTINE LDL_Cholesky_decomposition(A, L, D)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 1)), INTENT(OUT) :: L, D
        INTEGER :: i, j, N, k

        N = SIZE(A, 1)

        L = Identity_n(N)
        D = 0.d0

        DO j = 1, N
            D(j, j) = A(j, j) - DOT_PRODUCT(L(j, 1:j-1), L(j, 1:j-1) * [ (D(k,k), k = 1, j-1) ])


            DO i = j+1, N
                L(i, j) = (A(i, j) - DOT_PRODUCT(L(i, 1:j-1), L(j, 1:j-1) * [ (D(k,k), k = 1, j-1) ])) / D(j, j)
            END DO
        END DO

    END SUBROUTINE LDL_Cholesky_decomposition
    
    !> Incomplete Cholesky decomposition of a matrix A
    !> \[ A \approx LL^T \]
    !> This subroutine performs incomplete Cholesky decomposition of a given matrix **A**, where **L** is a lower triangular matrix and **U** is an upper triangular matrix.
    SUBROUTINE Incomplete_Cholesky_decomposition(A, L)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 1)), INTENT(OUT) :: L
        LOGICAL, DIMENSION(SIZE(A, 1), SIZE(A, 1)) :: S
        INTEGER :: i, j, N

        N = SIZE(A, 1)

        L = Identity_n(N)
        
        S = A /= 0

        DO i = 1, N
            DO j = 1, i-1
                IF (S(i,j)) L(i,j) = (A(i,j) - DOT_PRODUCT(L(i, 1:j-1), L(j, 1:j-1))) / L(j,j)
            END DO
            IF (S(i,i)) L(i,i) = SQRT(A(i,i) - DOT_PRODUCT(L(i, 1:i-1), L(i, 1:i-1)))
        END DO

    END SUBROUTINE Incomplete_Cholesky_decomposition

    !> QR decomposition of a matrix **A** using various methods
    !> \[ A = QR \]
    !> This subroutine performs QR decomposition of a given matrix **A** using the specified method (Householder, Givens, Classical Gram-Schmidt, or Modified Gram-Schmidt).
    !> The output matrices **Q** is an orthogonal matrix and **R** is an upper triangular matrix.
    SUBROUTINE QR_decomposition(A, method, Q, R)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)), INTENT(OUT) :: Q, R

        IF(method == "QR_Householder")THEN
            CALL QR_Householder_decomposition(A, Q, R)
        ELSE IF(method == "QR_Givens")THEN
            CALL QR_Givens_decomposition(A, Q, R)
        ELSE IF(method == "QR_Gram_Schmidt_Classical")THEN
            CALL QR_Gram_Schmidt_Classical_decomposition(A, Q, R) 
        ELSE IF(method == "QR_Gram_Schmidt_Modified")THEN
            CALL QR_Gram_Schmidt_Modified_decomposition(A, Q, R)
        END IF

    END SUBROUTINE QR_decomposition

    !> QR decomposition using Householder method
    SUBROUTINE QR_Householder_decomposition(A, Q, R)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)), INTENT(OUT) :: Q, R
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2))::Id, H, v_mat_tmp
        REAL(dp), DIMENSION(SIZE(A, 1)) :: v, u, x
        INTEGER :: N, i, j, k
        REAL(dp) :: alpha, w, signe, norm_u

        N= SIZE(A, 1)

        R = A

        Id = Identity_n(N)
        Q = Identity_n(N)

        DO k = 1, N

            x = 0.d0
            u = 0.d0
            v = 0.d0
            v_mat_tmp = 0.d0
            x(k:N) = R(K:N, K)

            alpha = NORM2(R(k:N, k))

            signe = - SIGN(alpha,x(k))
            u(k:N) =  x(k:N) - signe * Id(k:N, k)

            norm_u = NORM2(u)
            IF (norm_u < epsi) CYCLE
            v(k:N) = u(k:N) / norm_u

            w = 1.d0
            DO i = k, N
                DO j = k, N
                    v_mat_tmp(i, j) = v(i) * v(j)
                END DO
            END DO

            H = Id
            H(k:N, k:N) = Id(k:N, k:N) - (1.d0 + w) * v_mat_tmp(k:N, k:N)

            Q = MATMUL(Q, H)
        
            R(k:N, k:N) = MATMUL(H(k:N, k:N), R(k:N, k:N))
    
        END DO
    END SUBROUTINE QR_Householder_decomposition

    !> QR decomposition using Givens rotations
    SUBROUTINE QR_Givens_decomposition(A, Q, R)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)), INTENT(OUT) :: Q, R
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)) :: G
        INTEGER :: N, i, j

        N= SIZE(A, 1)

        R = A

        Q = Identity_n(N)

        DO j = 1, N-1
            DO i = j+1, N

                G = rotation_matrix(R, [i, j])

                R = MATMUL(G, R)

                Q = MATMUL(Q, TRANSPOSE(G))

            END DO
        END DO

    END SUBROUTINE QR_Givens_decomposition

    !> QR decomposition using Classical Gram-Schmidt method
    SUBROUTINE QR_Gram_Schmidt_Classical_decomposition(A, Q, R)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)), INTENT(OUT) :: Q, R
        REAL(dp), DIMENSION(SIZE(A, 1)) :: u
        INTEGER :: N, i, j

        N= SIZE(A, 1)
        Q = 0.d0
        R = 0.d0

        DO j = 1, N
            u = A(:, j)
            DO i = 1, j-1
                R(i, j) = DOT_PRODUCT(Q(:, i),A(:, j))
                u = u - (R(i, j) * Q(:, i))
            END DO
            R(j, j) = NORM2(u)
            Q(:, j) = u / R(j, j)
        END DO

    END SUBROUTINE QR_Gram_Schmidt_Classical_decomposition

    !> QR decomposition using Modified Gram-Schmidt method
    SUBROUTINE QR_Gram_Schmidt_Modified_decomposition(A, Q, R)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)), INTENT(OUT) :: Q, R
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 2)) :: u
        INTEGER :: N, i, j

        N = SIZE(A,1)
        u = A
        Q = 0.d0
        R = 0.d0

        DO i = 1, N
            R(i, i) = NORM2(u(:, i))
            Q(:, i) = u(:, i)/R(i, i)
            DO j = i+1, N
                R(i, j) = DOT_PRODUCT(Q(:, i),u(:, j))
                u(:, j) = u(:,j) - R(i, j)*Q(:, i)
            END DO
        END DO

    END SUBROUTINE QR_Gram_Schmidt_Modified_decomposition
    
END MODULE NAFPack_matrix_decomposition