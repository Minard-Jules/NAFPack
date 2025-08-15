MODULE NAFPack_Krylov_method

    USE NAFPack_constant
    USE NAFPack_matricielle

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE

CONTAINS

    SUBROUTINE lanczos(A, q1, m, Q, T)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: q1
        INTEGER, INTENT(IN) :: m
        REAL(dp), DIMENSION(:, :), INTENT(OUT) :: Q
        REAL(dp), DIMENSION(:, :), INTENT(OUT) :: T
        REAL(dp), DIMENSION(size(A, 1)) :: y, z
        REAL(dp), DIMENSION(m) :: alpha
        REAL(dp), DIMENSION(m - 1) :: beta
        INTEGER :: N, k

        N = size(A, 1)

        alpha = 0.d0
        beta = 0.d0
        Q = 0.d0

        Q(:, 1) = q1 / norm2(q1)

        DO k = 1, m
            IF (k == 1) THEN
                y = matmul(A, Q(:, k))
            ELSE
                y = matmul(A, Q(:, k)) - beta(k - 1) * Q(:, k - 1)
            END IF

            alpha(k) = dot_product(Q(:, k), y)
            z = y - alpha(k) * Q(:, k)
            
            IF (k < m) THEN
                beta(k) = norm2(z)
                IF (beta(k) < 1.0d-12) EXIT
                Q(:, k + 1) = z / beta(k)
            END IF
        END DO

        T = Make_Tridiagonal(beta, alpha, beta)

    END SUBROUTINE lanczos

    SUBROUTINE Arnoldi(A, q1, m, Q, H)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: q1
        INTEGER, INTENT(IN) :: m
        REAL(dp), DIMENSION(:, :), INTENT(OUT) :: Q
        REAL(dp), DIMENSION(:, :), INTENT(OUT) :: H
        INTEGER :: k, j, N

        N = size(A, 1)
        Q(:, 1) = q1 / norm2(q1)
        H = 0.0d0

        DO k = 2, m
            Q(:, k) = matmul(A, Q(:, k - 1))
            DO j = 1, k - 1
                H(j, k - 1) = dot_product(Q(:, j), Q(:, k))
                Q(:, k) = Q(:, k) - H(j, k - 1) * Q(:, j)
            END DO
            H(k, k - 1) = norm2(Q(:, k))
            Q(:, k) = Q(:, k) / H(k, k - 1)
        END DO

    END SUBROUTINE Arnoldi

    SUBROUTINE Arnoldi_MGS(A, q1, m, Q, H)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: q1
        INTEGER, INTENT(IN) :: m
        REAL(dp), DIMENSION(:, :), INTENT(OUT) :: Q ! (n, m+1)
        REAL(dp), DIMENSION(:, :), INTENT(OUT) :: H ! (m+1, m)
        INTEGER :: k, j, N
        REAL(dp), DIMENSION(size(A, 1)) :: w

        N = size(A, 1)
        Q(:, 1) = q1 / norm2(q1)
        H = 0.0d0

        DO k = 1, m
            ! w = A * q_k
            w = matmul(A, Q(:, k))

            ! Modified Gram-Schmidt orthonormalization
            DO j = 1, k
                H(j, k) = dot_product(Q(:, j), w)
                w = w - H(j, k) * Q(:, j)
            END DO

            H(k + 1, k) = norm2(w)
            IF (H(k + 1, k) > 0.0d0) THEN
                Q(:, k + 1) = w / H(k + 1, k)
            ELSE
                Q(:, k + 1) = 0.0d0
            END IF

        END DO

    END SUBROUTINE Arnoldi_MGS

END MODULE NAFPack_Krylov_method
