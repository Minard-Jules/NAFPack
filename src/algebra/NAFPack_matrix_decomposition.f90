MODULE NAFPack_matrix_decomposition

    USE NAFPack_constant
    USE NAFPack_matricielle

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: LU_decomposition, LDU_decomposition, Cholesky_decomposition, QR_decomposition

    CONTAINS

    SUBROUTINE LU_decomposition(A, L, U)

        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)), INTENT(OUT) :: L, U
        INTEGER :: i, j, k, N
        REAL(dp) :: S

        N = SIZE(A, 1)

        L = 0.d0
        U = 0.d0

        DO j = 1, N
            L(j, j) = 1.d0
          
            DO i = 1, j
                S = 0.d0
                DO k = 1, i-1
                    S = S + L(i, k) * U(k, j)
                END DO
                U(i, j) = A(i, j) - S
            END DO
    
            DO i = j+1, N
                S = 0.d0
                DO k = 1, j-1
                    S = S + L(i, k) * U(k, j)
                END DO
                L(i, j) = (A(i, j) - S) / U(j, j)
            END DO
        END DO

    END SUBROUTINE LU_decomposition

    SUBROUTINE LDU_decomposition(A, L, D, U)

        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)), INTENT(OUT) :: L, U, D
        INTEGER :: i, j, k, N
        REAL(dp) :: S

        N = SIZE(A, 1)

        L = 0.d0
        D = 0.d0
        U = 0.d0

        DO j = 1, N
            L(j, j) = 1.d0
            U(j, j) = 1.d0
            
            DO i = 1, j-1
                S = 0.d0
                DO k = 1, i-1
                    S = S + L(i, k) * D(k, k) * U(k, j)
                END DO
                U(i, j) = (A(i, j) - S) / D(i, i)
            END DO

            i = j
            S = 0
            DO k = 1, j-1
                S = S + L(j, k) * D(k, k) * U(k, j)
            END DO
            D(j, j) = A(j, j) - S
    
            DO i = j+1, N
                S = 0.d0
                DO k = 1, j-1
                    S = S + L(i, k) * D(k, k) * U(k, j)
                END DO
                L(i, j) = (A(i, j) - S) / D(j, j)
            END DO
        END DO

    END SUBROUTINE LDU_decomposition

    SUBROUTINE Cholesky_decomposition(A, L)

        REAL(dp),DIMENSION(: ,:), INTENT(IN) :: A
        REAL(dp),DIMENSION(SIZE(A, 1), SIZE(A, 1)), INTENT(OUT) :: L
        INTEGER :: i, j, k, N
        REAL(dp) :: S

        N = SIZE(A, 1)

        DO j = 1, N
            S = 0.d0
            DO k = 1, j-1
                S = S + L(j, k) * L(j, k)
            END DO
            L(j, j) = SQRT(A(j, j) - S)
        
            DO i = j+1, N
                S = 0.d0
                DO k = 1, j-1
                    S = S + L(i, k) * L(j, k)
                END DO
                L(i, j) = (A(i, j) - S) / L(j, j)
            END DO
        END DO

    END SUBROUTINE Cholesky_decomposition

    SUBROUTINE QR_decomposition(A, method, Q, R)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)), INTENT(OUT) :: Q, R

        IF(method == "QR_Householder")THEN
            CALL QR_Householder(A, Q, R)
        ELSE IF(method == "QR_Givens")THEN
            CALL QR_Givens(A, Q, R)
        ELSE IF(method == "QR_Gram_Schmidt_Classical")THEN
            CALL QR_Gram_Schmidt_Classical(A, Q, R) 
        ELSE IF(method == "QR_Gram_Schmidt_Modified")THEN
            CALL QR_Gram_Schmidt_Modified(A, Q, R)
        END IF

    END SUBROUTINE QR_decomposition

    SUBROUTINE QR_Householder(A, Q, R)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)), INTENT(OUT) :: Q, R
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2))::Id, H, v_mat_tmp
        REAL(dp), DIMENSION(SIZE(A, 1)) :: v, u, x
        INTEGER :: N, i, j, k
        REAL(dp) :: alpha, w, signe

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

            v(k:N) = u(k:N) / NORM2(u)

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
    END SUBROUTINE QR_Householder

    FUNCTION rotation_matrix(A,rotation) RESULT(G)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        INTEGER, DIMENSION(2), INTENT(IN) :: rotation
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 2)) :: G
        REAL(dp) :: frac, val_1, val_2
        INTEGER :: k, i, j

        i = rotation(1)
        j = rotation(2)

        G= 0.d0
        DO k = 1, SIZE(A, 1)
            G(k, k) = 1.d0
        END DO

        val_1 = A(j, j)
        val_2 = A(i, j)

        frac = SQRT(val_1**2 + val_2**2)

        G(i,i) = val_1 / frac
        G(j,j) = val_1 / frac
        G(i,j) = - val_2 / frac
        G(j,i) = val_2 / frac

    END FUNCTION rotation_matrix

    SUBROUTINE QR_Givens(A, Q, R)
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

    END SUBROUTINE QR_Givens

    SUBROUTINE QR_Gram_Schmidt_Classical(A, Q, R)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)), INTENT(OUT) :: Q, R
        REAL(dp), DIMENSION(SIZE(A, 1)) :: u
        INTEGER :: N, i, j

        N= SIZE(A, 1)

        Q = Identity_n(N)

        DO j = 1, N
            u = A(:, j)
            DO i = 1, j-1
                R(i, j) = DOT_PRODUCT(Q(:, i),A(:, j))
                u = u - (R(i, j) * Q(:, i))
            END DO
            R(j, j) = SQRT(SUM(u**2))
            Q(:, j) = u / R(j, j)
        END DO

    END SUBROUTINE QR_Gram_Schmidt_Classical

    SUBROUTINE QR_Gram_Schmidt_Modified(A, Q, R)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)), INTENT(OUT) :: Q, R
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 2)) :: u
        INTEGER :: N, i, j

        N = SIZE(A, 1)

        Q = Identity_n(N)
        DO i = 1, N
            u(:, i) = A(:,i)
        END DO

        DO i = 1, N
            R(i, i) = SQRT(SUM(u(:, i)**2))
            Q(:, i) = u(:, i)/R(i, i)
            DO j = i+1, N
                R(i, j) = DOT_PRODUCT(Q(:, i),u(:, j))
                u(:, j) = u(:,j) - R(i, j)*Q(:, i)
            END DO
        END DO

    END SUBROUTINE QR_Gram_Schmidt_Modified
    
END MODULE NAFPack_matrix_decomposition