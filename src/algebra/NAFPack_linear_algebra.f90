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

    !forward algorithm
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

    !backward algorithm
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

    FUNCTION Direct_methode(A, b, method, pivot_method, check) RESULT(x)

        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: pivot_method
        REAL(dp),DIMENSION(: ,:), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        LOGICAL, OPTIONAL, INTENT(IN) :: check
        REAL(dp),DIMENSION(SIZE(A,1)) :: x
        INTEGER :: N
        LOGICAL :: do_check = .TRUE.

        IF(PRESENT(check)) do_check = check

        N = SIZE(A, 1)
        IF(SIZE(A, 2) /= N) STOP "ERROR :: Matrix A not square"
        IF(SIZE(b, 1) /= N) STOP "ERROR :: Dimension mismatch in linear system"

        IF (.NOT. PRESENT(method)) THEN
            PRINT*, "WARNING :: No method specified for linear system, using LU decomposition"
            x = A_LU(A, b)
        END IF

        IF(method == "Gauss")THEN
            x = Gauss(A, b, pivot_method=pivot_method)
        ELSE IF(method == "A_LU")THEN
            x = A_LU(A, b)
        ELSE IF(method == "A_LDU")THEN
            x = A_LDU(A, b)
        ELSE IF(method == "Cholesky")THEN
            x = Cholesky(A, b, check = do_check)
        ELSE IF (INDEX(method, "QR") == 1) THEN
            x = A_QR(A, b, method = method)
        ELSE IF(method == "TDMA")THEN
            x = TDMA(A, b, check = do_check)
        ELSE
            STOP "ERROR : Wrong method for linear system (direct_methode)"
        END IF

    END FUNCTION Direct_methode

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
                    P(k, k) = 0.0_dp
                    P(lmax, lmax) = 0.0_dp
                    P(k, lmax) = 1.0_dp
                    P(lmax, k) = 1.0_dp
                END IF

                A_tmp = MATMUL(P, A_tmp)
                b_tmp = MATMUL(P, b_tmp)
            ELSE IF(pivot_method == "total") THEN
                ! Find max abs element in submatrix
                vlmax_total = MAXLOC(ABS(A_tmp(k:N,k:N)))
                lmax = vlmax_total(1) + k - 1
                cmax = vlmax_total(2) + k - 1

                P = Identity_n(N)
                Q = Identity_n(N)

                ! permute line if necessary
                IF (lmax /= k) THEN
                    P(k, k) = 0.0_dp
                    P(lmax, lmax) = 0.0_dp
                    P(k, lmax) = 1.0_dp
                    P(lmax, k) = 1.0_dp
                END IF

                ! permute column if necessary
                IF (cmax /= k) THEN
                    Q(k, k) = 0.0_dp
                    Q(cmax, cmax) = 0.0_dp
                    Q(k, cmax) = 1.0_dp
                    Q(cmax, k) = 1.0_dp
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
            IF(MINVAL(lambda) < epsi) STOP "ERROR :: A is not a definite matrix (Cholesky)"
            IF(MAXVAL(ABS(A - TRANSPOSE(A))) > epsi) STOP "ERROR :: A is not symmetric (Cholesky)"
        END IF
        
        CALL Cholesky_decomposition(A, L)
          
        y = forward(L, b)
    
        x = backward(TRANSPOSE(L), y)

    END FUNCTION Cholesky

    FUNCTION A_QR(A, b, method) RESULT(x)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1) ,SIZE(A, 2)) :: Q, R

        CALL QR_decomposition(A, method, Q, R)

        x = backward(R, MATMUL(TRANSPOSE(Q), b))

    END FUNCTION A_QR

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
    INTEGER :: k, max_iter_choice, N

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

        residu = b - MATMUL(A, x_new)

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

    SUBROUTINE Eigen(A, lambda, vp, method, k)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        INTEGER, OPTIONAL, INTENT(IN) :: k
        REAL(dp), DIMENSION(:, :), OPTIONAL, INTENT(OUT) :: vp
        REAL(dp), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: A_tmp
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: vp_tmp
        CHARACTER(LEN = 50) :: base_method
        INTEGER :: N, i, k_max, pos

        IF(PRESENT(k)) THEN
            IF (k <= 0) STOP "ERROR :: k must be a positive integer"
            k_max = k
        ELSE
            k_max = kmax
        END IF
        
        N = SIZE(A, 1)
        IF(SIZE(A, 2) /= N) STOP "ERROR :: Matrix A not square"

        IF(SIZE(lambda, 1) /= N) STOP "ERROR :: dimension lambda"
        IF(PRESENT(vp) .AND. (SIZE(vp, 1) /= N .OR. SIZE(vp, 2) /= N)) STOP "ERROR :: dimension vp"

        IF(method == "Power_iteration")THEN

            A_tmp = A
            DO i=1, N
                CALL Power_iteration(A_tmp, lambda(i), vp_tmp(i, :), k_max)
                A_tmp = deflation(A_tmp, lambda(i), vp_tmp(i, :), k_max)
            END DO

            IF(PRESENT(vp)) vp = vp_tmp

        ELSE IF (INDEX(method, "QR") == 1) THEN

            IF(PRESENT(vp)) vp = 0
            IF(PRESENT(vp)) PRINT*, "WARNING :: No solution for eigenvectors with the QR method"

            pos = INDEX(TRIM(method), "_Shifted")

            IF (pos > 0 .AND. pos + 7 == LEN_TRIM(method)) THEN
                base_method = method(:pos - 1)
                CALL Eigen_QR_Shifted(A, lambda, base_method, N, k_max)
            ELSE
                CALL Eigen_QR(A, lambda, method, N, k_max)
            END IF

        ELSE
            STOP "ERROR :: Wrong method for Eigen"
        END IF

    END SUBROUTINE Eigen

    SUBROUTINE Eigen_QR(A,lambda,method, N, k)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), INTENT(IN) :: method
        INTEGER, INTENT(IN) :: N, k
        REAL(dp), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(dp), DIMENSION(SIZE(A, 1)) :: lambda_old
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: A_tmp, Q, R
        REAL(dp) :: diff
        INTEGER :: i, j

        A_tmp = A

        DO i = 1, k

            lambda_old = lambda

            CALL QR_decomposition(A_tmp, method, Q, R)

            A_tmp = MATMUL(R, Q)

            diff = ABS(A_tmp(2, 1))
            DO j = 3, N
                IF (MAXVAL(ABS(A_tmp(j, 1:j-1))) > diff) THEN
                    diff = MAXVAL(ABS(A_tmp(j, 1:j-1)))
                END IF
            END DO


            IF(i == k)THEN
                PRINT*, " WARNING :: non-convergence of the QR Algorithm for eigenvalues "//method
                PRINT*, "convergence = ", diff
                EXIT
            END IF

            IF(diff <= epsi) EXIT
        END DO

        ! Extract eigenvalues
        lambda = [(A_tmp(i,i), i=1,N)]

    END SUBROUTINE Eigen_QR

    SUBROUTINE Eigen_QR_Shifted(A, lambda, method, N, k)
        INTEGER, INTENT(IN) :: N, k
        CHARACTER(LEN = *), INTENT(IN) :: method
        REAL(dp), DIMENSION(N, N), INTENT(IN) :: A
        REAL(dp), DIMENSION(N), INTENT(OUT) :: lambda
        INTEGER :: i, j
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: A_tmp, Q, R, Id
        REAL(dp) :: shift, diff

        A_tmp = A
        Id = Identity_n(N)

        DO i = 1, k
            !choice of shift: last diagonal element
            shift = A_tmp(N,N)

            ! Gap : A - µI
            A_tmp = A_tmp - shift * Id

            ! QR Decomposition : A - µI = Q * R
            CALL QR_decomposition(A_tmp, method, Q, R)

            ! A = RQ + µI
            A_tmp = MATMUL(R, Q) + shift * Id

            diff = ABS(A_tmp(2, 1))
            DO j = 3, N
                IF (MAXVAL(ABS(A_tmp(j, 1:j-1))) > diff) THEN
                    diff = MAXVAL(ABS(A_tmp(j, 1:j-1)))
                END IF
            END DO
            
            IF(i == k)THEN
                PRINT*, " WARNING :: non-convergence of the Shifted QR Algorithm for eigenvalues "//method
                PRINT*, "convergence = ", diff
                EXIT
            END IF

            IF(diff <= epsi) EXIT
            
        END DO

        ! Extract eigenvalues
        lambda = [(A_tmp(i,i), i=1,N)]

    END SUBROUTINE Eigen_QR_Shifted


    SUBROUTINE Power_iteration(A, lambda, vp, k)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        INTEGER, INTENT(IN) :: k
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

        DO i = 1, k
            u = normalise(vp_tmp)
            vp_tmp = MATMUL(A, u)
            lambda = DOT_PRODUCT(vp_tmp, u)
            IF (NORM2(r) <= epsi)EXIT
            r = vp_tmp - lambda * u
            IF(i == k)THEN
                PRINT*, "WARNING :: non-convergence of the power iteration method"
            END IF
        END DO

        vp = u

    END SUBROUTINE Power_iteration
  
    FUNCTION deflation(A, lambda, vp, k) RESULT(result)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: vp
        REAL(dp), INTENT(IN) :: lambda
        INTEGER, INTENT(IN) :: k
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 1)) :: result
        REAL(dp), DIMENSION(SIZE(A, 1)) :: wp
        INTEGER :: i, j, N
        REAL(dp) :: lambda1
        
        N = SIZE(A, 1)
        result = A
        
        CALL Power_iteration(transpose(A), lambda1, wp, k)
        DO i = 1, N 
            DO j = 1, N
                result(i, j) = result(i, j) - (lambda * vp(i) * wp(j)) / DOT_PRODUCT(vp, wp)
            END DO
        END DO    

    END FUNCTION deflation

END MODULE NAFPack_linear_algebra