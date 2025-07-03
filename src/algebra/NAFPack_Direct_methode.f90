MODULE NAFPack_Direct_methode

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    USE NAFPack_matricielle
    USE NAFPack_Eigen

    IMPLICIT NONE
    
    PRIVATE
    
    PUBLIC :: Gauss, A_LU, A_LDU, Cholesky, A_LDL_Cholesky, A_QR, TDMA, Faddeev_Leverrier

    CONTAINS

    !> Gaussian elimination method
    !>
    !> This function implements the Gaussian elimination method for solving linear systems.
    !> It performs partial or total pivoting based on the user's choice.
    FUNCTION Gauss(A, b, pivot_method) RESULT(x)
        
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: pivot_method
        REAL(dp), DIMENSION(SIZE(A,1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: A_tmp, P, Q, Q_final
        REAL(dp), DIMENSION(SIZE(b)) :: b_tmp
        INTEGER, DIMENSION(1) :: vlmax_partial
        INTEGER, DIMENSION(2) :: vlmax_total
        REAL(dp) :: m , pivot
        INTEGER :: k, i, N, lmax, cmax

        A_tmp = A
        b_tmp = b

        N=SIZE(A_tmp, 1)

        Q_final = Identity_n(N)
        
        DO k = 1, N-1
            IF (.NOT. PRESENT(pivot_method)) THEN 
                PRINT*, "WARNING :: No pivot method specified, using normal pivot"
            END IF

            IF(pivot_method == "partial") THEN
                ! Find the maximum absolute value in the column from row k to N
                vlmax_partial = MAXLOC(ABS(A_tmp(k:N, k)))
                lmax = vlmax_partial(1) + k - 1
                
                !calculate permutation matrix P
                P = Identity_n(N)
                IF (k /= lmax) THEN
                    P = rotation_matrix(P, [k, lmax])
                END IF

                A_tmp = MATMUL(P, A_tmp)
                b_tmp = MATMUL(P, b_tmp)
            ELSE IF(pivot_method == "total") THEN
                ! Find max abs element in submatrix
                vlmax_total = MAXLOC(ABS(A_tmp(k:N,k:N)))
                lmax = vlmax_total(1) + k - 1
                cmax = vlmax_total(2) + k - 1

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
                A_tmp = MATMUL(P, A_tmp)
                A_tmp = MATMUL(A_tmp, Q)

                b_tmp = MATMUL(P, b_tmp)
            END IF

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
        IF (pivot_method == "total") x = MATMUL(Q_final, x)

    END FUNCTION Gauss

    !> LU decomposition method
    !>
    !> This function implements the LU decomposition method for solving linear systems.
    FUNCTION A_LU(A, b) RESULT(x)

        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        REAL(dp),DIMENSION(SIZE(A, 1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: L, U
        REAL(dp),DIMENSION(SIZE(A, 1)) :: y

        CALL LU_decomposition(A, L, U)
    
        y = forward(L, b)
        
        x = backward(U, y)

    END FUNCTION A_LU

    !> LDU decomposition method
    !>
    !> This function implements the LDU decomposition method for solving linear systems.
    FUNCTION A_LDU(A, b) RESULT(x)

        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        REAL(dp),DIMENSION(SIZE(A,1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: L, U, D
        REAL(dp),DIMENSION(SIZE(A, 1)) :: y, z

        CALL LDU_decomposition(A, L, D, U)
    
        z = forward(L, b)

        y = forward(D, z)
        
        x = backward(U, y)
    
    END FUNCTION A_LDU

    !> Cholesky method
    !>
    !> This function implements the Cholesky decomposition method for solving linear systems.
    FUNCTION Cholesky(A, b, check) RESULT(x)

        REAL(dp),DIMENSION(: ,:), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        LOGICAL, OPTIONAL, INTENT(IN) :: check
        REAL(dp),DIMENSION(SIZE(A,1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1), SIZE(A, 1)) :: L
        REAL(dp),DIMENSION(SIZE(A, 1)) :: y, lambda
        LOGICAL :: do_check = .TRUE.

        IF(PRESENT(check)) do_check = check

        IF(do_check) THEN
            PRINT*, "Checking if the matrix A is positive definite and symmetric"
            CALL Eigen(A, lambda, method = "Power_iteration")
            IF(MINVAL(lambda) < 0) STOP "ERROR :: A is not a definite matrix (Cholesky)"
            IF(MAXVAL(ABS(A - TRANSPOSE(A))) > epsi) STOP "ERROR :: A is not symmetric (Cholesky)"
        END IF
        
        CALL Cholesky_decomposition(A, L)
          
        y = forward(L, b)
    
        x = backward(TRANSPOSE(L), y)

    END FUNCTION Cholesky

    !> Alternative Cholesky method 
    !>
    !> This function implements the alternative Cholesky decomposition (LDL^T) method for solving linear systems.
    FUNCTION A_LDL_Cholesky(A, b, check) RESULT(x)

        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        LOGICAL, OPTIONAL, INTENT(IN) :: check
        REAL(dp),DIMENSION(SIZE(A,1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: L, D
        REAL(dp),DIMENSION(SIZE(A, 1)) :: y, z, lambda
        LOGICAL :: do_check = .TRUE.

        IF(PRESENT(check)) do_check = check

        IF (do_check) THEN
            PRINT*, "Checking if the matrix A is positive definite and symmetric"
            CALL Eigen(A, lambda, method = "Power_iteration")
            IF(MINVAL(lambda) < 0) STOP "ERROR :: A is not a definite matrix (Cholesky)"
            IF(MAXVAL(ABS(A - TRANSPOSE(A))) > epsi) STOP "ERROR :: A is not symmetric (Cholesky)"
        END IF

        CALL LDL_Cholesky_decomposition(A, L, D)

        z = forward(L, b)

        y = forward(D, z)

        x = backward(TRANSPOSE(L), y)
    
    END FUNCTION A_LDL_Cholesky

    !> QR decomposition method
    !>
    !> This function implements the QR decomposition method for solving linear systems.
    !> It allows the user to choose the method for QR decomposition (e.g., Householder).
    FUNCTION A_QR(A, b, method) RESULT(x)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)) :: Q, R

        CALL QR_decomposition(A, method, Q, R)

        x = backward(R, MATMUL(TRANSPOSE(Q), b))

    END FUNCTION A_QR

    !> TDMA method (or Thomas algorithm)
    !>
    !> This function implements the Thomas algorithm for solving tridiagonal systems.
    FUNCTION TDMA(A, b, check) RESULT(x)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        LOGICAL, OPTIONAL, INTENT(IN) :: check
        REAL(dp), DIMENSION(SIZE(A,1)) :: x
        REAL(dp), DIMENSION(SIZE(A,1)) :: alpha, beta
        REAL(dp) :: denom
        INTEGER :: n, i, j
        LOGICAL :: do_check = .TRUE.

        N = SIZE(A,1)

        IF(PRESENT(check)) do_check = check 

        IF (do_check) THEN
            PRINT*, "Checking if the matrix A is tridiagonal"
            DO i = 1, N
                DO j = 1, N
                    IF (ABS(i-j) > 1) THEN
                        IF (ABS(A(i,j)) > epsi) STOP "ERROR :: Matrix is not tridiagonal"
                    END IF
                END DO
            END DO
        END IF

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

    END FUNCTION TDMA

    !> Faddeev-Leverrier method for computing the inverse of a matrix and its characteristic polynomial coefficients
    !> \[ P(\lambda) = \det(\lambda I - A) = c(1)\lambda^N + c(2)\lambda^{N-1} + \ldots + c(N+1) \]
    SUBROUTINE Faddeev_Leverrier(A, c, Ainv, success, check)
        INTEGER, PARAMETER :: dp = KIND(1.0d0)
        REAL(dp), DIMENSION(:, :), INTENT(IN)  :: A
        LOGICAL, OPTIONAL, INTENT(IN) :: check
        REAL(dp), DIMENSION(:),     INTENT(OUT) :: c
        REAL(dp), DIMENSION(SIZE(A,1), SIZE(A,1)), OPTIONAL, INTENT(OUT) :: Ainv
        LOGICAL, OPTIONAL, INTENT(OUT) :: success
        REAL(dp), DIMENSION(SIZE(A,1), SIZE(A,1)) :: Bk, I, B_Nm1, AB
        LOGICAL :: do_check = .TRUE.
        INTEGER :: N, k

        N = SIZE(A,1)

        IF (PRESENT(check)) do_check = check

        IF (do_check) THEN
            PRINT*, "Checking if the matrix A is square and size of c is correct"
            IF (SIZE(A,2) /= N .OR. SIZE(c) < N+1) THEN
                PRINT *, "Error : Matrix A must be square and size of c must be at least N+1"
                STOP
            END IF
        END IF

        ! Initialization
        I = Identity_n(N)
        c = 0.0_dp
        c(1) = 1.0_dp
        c(2) = -Trace(A)
        Bk = A + c(2)*I

        DO k = 2, N
            AB = MATMUL(A, Bk)
            c(k+1) = -Trace(AB) / REAL(k, dp)
            Bk = AB + c(k+1)*I
            IF (k == N-1 .AND. PRESENT(Ainv)) B_Nm1 = -Bk
        END DO

        IF (PRESENT(Ainv) .AND. PRESENT(success)) THEN
            IF (ABS(c(N+1)) < 1.0e-12_dp) THEN
                success = .FALSE.
                Ainv = 0.0_dp
            ELSE
                success = .TRUE.
                Ainv = B_Nm1 / c(N+1)
            END IF
        ELSE IF (PRESENT(Ainv)) THEN
            IF (ABS(c(N+1)) < 1.0e-12_dp) THEN
                Ainv = 0.0_dp
            ELSE
                Ainv = B_Nm1 / c(N+1)
            END IF
        END IF

    END SUBROUTINE Faddeev_Leverrier

END MODULE NAFPack_Direct_methode