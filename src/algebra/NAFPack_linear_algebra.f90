MODULE NAFPack_linear_algebra

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    USE NAFPack_matricielle

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: Direct_methode, Iterative_methods
    PUBLIC :: Eigen

    CONTAINS


!================== Linear System =======================================================

    !descent algorithm
    FUNCTION descent(L, b) RESULT(y_descent)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: L
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(L, 1)) :: y_descent
        REAL(dp) :: S
        INTEGER :: i, j, N
    
        N = SIZE(L, 1)
    
        y_descent(1) = b(1) / L(1, 1)
    
        DO i = 2, N
            S = 0
            DO j = 1, i-1
                S = S + L(i, j) * y_descent(j)
            END DO
            y_descent(i) = (b(i) - S) / L(i, i)
        END DO
    
    END FUNCTION descent

    !ascent algorithm
    FUNCTION ascent(U, y) RESULT(result)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: U
        REAL(dp), DIMENSION(:), INTENT(IN) :: y
        REAL(dp), DIMENSION(SIZE(U, 1)) :: result
        REAL(dp) :: S
        INTEGER :: i, j, N

        N = SIZE(U, 1)

        result(N) = y(N) / U(N, N)

        DO i = N-1, 1, -1
            S = 0
            DO j = i+1, N
                S = S + U(i, j) * result(j)
            END DO
            result(i) = (y(i) - S) / U(i, i)
        END DO

    END FUNCTION ascent

    SUBROUTINE exchange_vector(Line1, Line2)

        REAL(dp), DIMENSION(:), INTENT(INOUT) :: Line1, Line2
        REAL(dp), DIMENSION(SIZE(Line1)) :: tmp
        
        tmp = Line1
        Line1 = Line2
        Line2 = tmp

    END SUBROUTINE exchange_vector

    SUBROUTINE scalar_exchange(scal1, scal2)

        REAL(dp) :: scal1, scal2
        REAL(dp) :: tmp
        
        tmp = scal1
        scal1 = scal2
        scal2 = tmp

    END SUBROUTINE scalar_exchange
    

!################## direct methode ######################################################

    FUNCTION Direct_methode(A, b, method) RESULT(x)

        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        REAL(dp),DIMENSION(: ,:), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        REAL(dp),DIMENSION(SIZE(A,1)) :: x
        INTEGER :: N

        N = SIZE(A, 1)
        IF(SIZE(A, 2) /= N) STOP "ERROR :: Matrix A not square"

        IF(method == "Gauss")THEN
            x = Gauss(A, b)
        ELSE IF(method == "Gauss_pivot")THEN
            x = Gauss_pivot(A, b)
        ELSE IF(method == "A_LU")THEN
            x = A_LU(A, b)
        ELSE IF(method == "A_LDU")THEN
            x = A_LDU(A, b)
        ELSE IF(method == "Cholesky")THEN
            x = Cholesky(A, b)
        ELSE IF(method == "QR_Householder" .OR. &
                method == "QR_Givens" .OR. &
                method == "QR_Gram_Schmidt_Classical".OR. &
                method == "QR_Gram_Schmidt_Modified")THEN
            x = A_QR(A, b, method = method)
        ELSE IF(method == "TDMA")THEN
            x = TDMA(A, b)
        ELSE
            STOP "ERROR : Wrong method for linear system (direct_methode)"
        END IF

    END FUNCTION Direct_methode

    FUNCTION Gauss(A, b) RESULT(x)
        
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A,1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: A_tmp
        REAL(dp), DIMENSION(SIZE(b)) :: b_tmp
        REAL(dp) :: q, pivot
        INTEGER :: k, i, j, N

        A_tmp = A
        b_tmp = b

        N=SIZE(A_tmp, 1)
        
        DO k = 1, N-1 
            DO i = k+1, N
                pivot = A_tmp(k, k)
                if (pivot == 0.) stop "ERROR :: A is a singular matrix"

                q = A_tmp(i, k) / pivot
                A_tmp(i, k) = 0 
                DO j = k+1, N
                    A_tmp(i, j) = A_tmp(i, j) - q * A_tmp(k, j)
                END DO

                b_tmp(i) = b_tmp(i) - q * b_tmp(k)
            END DO
        END DO

        x = ascent(A_tmp, b_tmp)

    END FUNCTION Gauss

    FUNCTION Gauss_pivot(A, b) RESULT(x)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A,1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: A_tmp
        REAL(dp), DIMENSION(SIZE(b)) :: b_tmp
        INTEGER, DIMENSION(1) :: vlmax
        REAL(dp) :: q, pivot
        INTEGER :: k, i, j, N, lmax

        A_tmp = A
        b_tmp = b

        N = SIZE(A, 1)
        
        DO k = 1, N-1 
            DO i = k+1, N

                vlmax = MAXLOC(ABS(A_tmp(k:n, k)))
                lmax = vlmax(1) + k - 1
                pivot = A_tmp(lmax, k)
                if (pivot == 0.) stop "ERROR :: A is a singular matrix"

                CALL exchange_vector(A_tmp(k, :), A_tmp(lmax, :))
                CALL scalar_exchange(b_tmp(k), b_tmp(lmax))

                q = A_tmp(i, k) / pivot
                A_tmp(i, k) = 0
                DO j = k+1, N
                    A_tmp(i, j) = A_tmp(i, j) - q * A_tmp(k, j)
                END DO

                b_tmp(i) = b_tmp(i) - q * b_tmp(k)
            END DO
        END DO

        x = ascent(A_tmp, b_tmp)

    END FUNCTION Gauss_pivot

    FUNCTION A_LU(A, b) RESULT(x)

        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        REAL(dp),DIMENSION(SIZE(A, 1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: L, U
        REAL(dp),DIMENSION(SIZE(A, 1)) :: y

        CALL LU_decomposition(A, L, U)
    
        y = descent(L, b)
        
        x = ascent(U, y)

    END FUNCTION A_LU

    FUNCTION A_LDU(A, b) RESULT(x)

        REAL(dp),DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        REAL(dp),DIMENSION(SIZE(A,1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: L, U, D
        REAL(dp),DIMENSION(SIZE(A, 1)) :: y, z

        CALL LDU_decomposition(A, L, D, U)
    
        z = descent(L, b)

        y = descent(D, z)
        
        x = ascent(U, y)
    
    END FUNCTION A_LDU

    FUNCTION Cholesky(A, b) RESULT(x)

        REAL(dp),DIMENSION(: ,:), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        REAL(dp),DIMENSION(SIZE(A,1)) :: x
        REAL(dp),DIMENSION(SIZE(A, 1), SIZE(A, 1)) :: L
        REAL(dp),DIMENSION(SIZE(A, 1)) :: y, lambda

        CALL Eigen(A, lambda, method = "Power_iteration")
        IF(MINVAL(lambda) < 0) STOP "ERROR :: A is not a definite matrix (Cholesky)"
        IF(MAXVAL(ABS(A - TRANSPOSE(A))) > epsi) STOP "ERROR :: A is not symmetric (Cholesky)"
        CALL Cholesky_decomposition(A, L)
          
        y = descent(L, b)
    
        x = ascent(TRANSPOSE(L), y)

    END FUNCTION Cholesky

    FUNCTION A_QR(A, b, method) RESULT(x)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)) :: Q, R

        CALL QR_decomposition(A, method, Q, R)

        x = ascent(R, MATMUL(TRANSPOSE(Q), b))

    END FUNCTION A_QR

    FUNCTION TDMA(A, b) RESULT(x)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(SIZE(A,1)) :: x
        REAL(dp), DIMENSION(SIZE(A,1)) :: alpha, beta
        REAL(dp) :: denom
        INTEGER :: n, i

        n = SIZE(A,1)
        IF (SIZE(A,2) /= n) STOP "ERROR :: Matrix A not square"
        IF (SIZE(b,1) /= n) STOP "ERROR :: Dimension mismatch in TDMA"

        alpha(1) = A(1,2) / A(1,1)
        beta(1) = b(1) / A(1,1)
        DO i = 2, n
            denom = A(i,i) - A(i,i-1)*alpha(i-1)
            alpha(i) = 0.0_dp
            IF (i < n) alpha(i) = A(i,i+1) / denom
            beta(i) = (b(i) - A(i,i-1)*beta(i-1)) / denom
        END DO

        ! Back substitution
        x(n) = beta(n)
        DO i = n-1, 1, -1
            x(i) = beta(i) - alpha(i)*x(i+1)
        END DO

    END FUNCTION TDMA


!################## Iterative methods ###################################################

FUNCTION Iterative_methods(A, b, method, x_init, max_iter, omega) RESULT(x)

    REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
    REAL(dp), DIMENSION(:), INTENT(IN) :: b
    REAL(dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: x_init
    CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
    INTEGER, OPTIONAL, INTENT(IN) :: max_iter
    REAL(dp), OPTIONAL, INTENT(IN) :: omega
    REAL(dp), DIMENSION(SIZE(A, 1)) :: x
    REAL(dp), DIMENSION(SIZE(A, 1)) :: x0, x_new, residu
    INTEGER :: k, max_iter_choice, i, N

    N = SIZE(A, 1)
    IF(SIZE(A, 2) /= N) STOP "ERROR :: Matrix A not square"

    IF (PRESENT(max_iter))THEN
        max_iter_choice = max_iter
    ELSE
        max_iter_choice = kmax
    END IF

    IF (PRESENT(x_init))THEN
        IF (SIZE(x_init, 1) /= SIZE(A, 1)) STOP "ERROR : Dimension of x_init different from A"
        x0 = x_init
    ELSE
        x0 = 0.d0
    END IF

    DO k = 1, max_iter_choice

        IF(k == kmax) THEN
            PRINT*, "WARNING :: non-convergence of the iterative method "//method
        END IF
        
        IF(method == "Jacobi")THEN
            CALL Jacobi(A, b , x0, x_new)
        ELSE IF(method == "Gauss_Seidel")THEN
            CALL Gauss_Seidel(A, b, x0, x_new)
        ELSE IF(method == "SOR")THEN
            IF (PRESENT(omega))THEN
                CALL SOR(A, b, x0, x_new, omega)
            ELSE
                CALL SOR(A, b, x0, x_new, 1.d0)
            END IF
        ELSE
            STOP "ERROR : Wrong method for linear system (Iterative_methods)"
        END IF

        DO i = 1, N
            residu(i) = b(i) - DOT_PRODUCT(A(i, :), x_new)
        END DO

        IF (NORM2(residu) < epsi) EXIT

        x0 = x_new

    END DO

    x = x_new

END FUNCTION Iterative_methods

SUBROUTINE Jacobi(A, b, x0, x)
    REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
    REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
    REAL(dp), DIMENSION(SIZE(A, 1)), INTENT(OUT) :: x
    INTEGER :: i, N

    N = SIZE(A, 1)

    DO i = 1, N
        x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x0(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
        x(i) = x(i) / A(i, i)
    END DO

END SUBROUTINE Jacobi

SUBROUTINE Gauss_Seidel(A, b, x0, x)
    REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
    REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
    REAL(dp), DIMENSION(SIZE(A, 1)), INTENT(OUT) :: x
    INTEGER :: i, N

    N = SIZE(A, 1)

    DO i = 1, N
        x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
        x(i) = x(i) / A(i, i)
    END DO

END SUBROUTINE Gauss_Seidel

SUBROUTINE SOR(A, b, x0, x, omega)
    REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
    REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
    REAL(dp), INTENT(IN) :: omega
    REAL(dp), DIMENSION(SIZE(A, 1)), INTENT(OUT) :: x
    INTEGER :: i, N

    N = SIZE(A, 1)

    DO i = 1, N
        x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
        x(i) = omega*(x(i) / A(i, i) - x0(i))
        x(i) = x(i) + x0(i)
    END DO

END SUBROUTINE SOR


!================== Eigen ===============================================================

    SUBROUTINE Eigen(A, lambda, vp, method)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        REAL(dp), DIMENSION(:, :), OPTIONAL, INTENT(OUT) :: vp
        REAL(dp), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: A_tmp
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: vp_tmp
        INTEGER :: N, i
        
        N = SIZE(A, 1)
        IF(SIZE(A, 2) /= N) STOP "ERROR :: Matrix A not square"

        IF(SIZE(lambda, 1) /= N) STOP "ERROR :: dimension lambda"
        IF(PRESENT(vp) .AND. (SIZE(vp, 1) /= N .OR. SIZE(vp, 2) /= N)) STOP "ERROR :: dimension vp"

        IF(method == "Power_iteration")THEN

            A_tmp = A
            DO i=1, N
                CALL Power_iteration(A_tmp, lambda(i), vp_tmp(i, :))
                A_tmp = deflation(A_tmp, lambda(i), vp_tmp(i, :))
            END DO

            IF(PRESENT(vp)) vp = vp_tmp

        ELSE IF(method == "QR_Householder" .OR. &
            method == "QR_Givens" .OR. &
            method == "QR_Gram_Schmidt_Classical" .OR. &
            method == "QR_Gram_Schmidt_Modified")THEN

            IF(PRESENT(vp)) vp = 0
            IF(PRESENT(vp)) PRINT*, "WARNING :: No solution for eigenvectors with the QR method"

            CALL Eigen_QR(A, lambda, method, N)
        ELSE
            STOP "ERROR :: Wrong method for Eigen"
        END IF

    END SUBROUTINE Eigen

    SUBROUTINE Eigen_QR(A,lambda,method, N)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), INTENT(IN) :: method
        INTEGER, INTENT(IN) :: N
        REAL(dp), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(dp), DIMENSION(SIZE(A, 1)) :: lambda_old
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)) :: Q, R
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: A_tmp
        REAL(dp) :: diff
        INTEGER :: i, j

        A_tmp = A

        DO i = 1, kmax

            lambda_old = lambda

            CALL QR_decomposition(A_tmp, method, Q, R)

            A_tmp = MATMUL(R, Q)

            diff = 1
            DO j = 1, N
                lambda(j) = A_tmp(j, j)
            END DO

            IF(i == kmax)THEN
                PRINT*, "WARNING :: non-convergence of the QR method for eigenvalues "//method
                EXIT
            END IF

            diff = ABS(A_tmp(N, N-1))

            IF(diff <= epsi) EXIT
        END DO

    END SUBROUTINE Eigen_QR

    SUBROUTINE Power_iteration(A, lambda, vp)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(OUT) :: vp
        REAL(dp), INTENT(OUT) :: lambda
        REAL(dp), DIMENSION(SIZE(A, 1)) :: u, vp_tmp, r
        INTEGER :: i, N

        N = SIZE(A, 1)
        CALL RANDOM_NUMBER(u)

        u = normalise(u)
        vp_tmp = MATMUL(A, u)
        lambda = DOT_PRODUCT(vp_tmp, u)
        r = vp_tmp - lambda * u

        DO i = 1, kmax
            u = normalise(vp_tmp)
            vp_tmp = MATMUL(A, u)
            lambda = DOT_PRODUCT(vp_tmp, u)
            IF (NORM2(r) <= epsi)EXIT
            r = vp_tmp - lambda * u
            IF(i == kmax)THEN
                PRINT*, "WARNING :: non-convergence of the power iteration method"
            END IF
        END DO

        vp = u

    END SUBROUTINE Power_iteration
  
    FUNCTION deflation(A, lambda, vp) RESULT(result)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: vp
        REAL(dp), INTENT(IN) :: lambda
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: result
        REAL(dp), DIMENSION(SIZE(A, 1)) :: wp
        INTEGER :: i, j, N
        REAL(dp) :: lambda1
        
        N = SIZE(A, 1)
        result = A
        
        CALL Power_iteration(transpose(A), lambda1, wp)
        DO i = 1, N 
            DO j = 1, N
                result(i, j) = result(i, j) - (lambda * vp(i) * wp(j)) / DOT_PRODUCT(vp, wp)
            END DO
        END DO    

    END FUNCTION deflation

END MODULE NAFPack_linear_algebra