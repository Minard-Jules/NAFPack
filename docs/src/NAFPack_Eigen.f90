!> Module for eigenvalue and eigenvector computations in NAFPack
MODULE NAFPack_Eigen

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    USE NAFPack_matricielle

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: Eigen

    CONTAINS

    !================== Eigen ===============================================================

    !> Computes the eigenvalues and eigenvectors of a matrix A
    !> \[ A * \vec{v} = \lambda * \vec{v} \]
    !> with **A** a square matrix, **λ** the eigenvalue, and **v** the eigenvector.
    !> This subroutine allows you to choose the method for computing eigenvalues and eigenvectors:
    !>
    !> - Power iteration
    !> - QR algorithm (with or without shift)
    !> The default method is Power iteration.
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

    !> QR algorithm for computing eigenvalues
    !>
    !> This subroutine implements the QR algorithm for computing the eigenvalues of a matrix.
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

    !> Shifted QR algorithm for computing eigenvalues
    !>
    !> This subroutine implements the shifted QR algorithm for computing the eigenvalues of a matrix.
    !> The shift is chosen as the last diagonal element of the matrix.
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

    !> Power iteration method for computing the dominant eigenvalue and eigenvector
    !>
    !> This subroutine implements the power iteration method for finding the dominant eigenvalue and eigenvector of a matrix.
    !> It iteratively computes the eigenvector and eigenvalue until convergence
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
  
    !> Deflation method for removing the influence of an eigenvalue and eigenvector
    !>
    !> This function performs deflation on a matrix A by removing the influence of an eigenvalue and its corresponding eigenvector.
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
    

END MODULE NAFPack_Eigen