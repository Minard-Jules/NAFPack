!> Module for resolving linear systems
!> 
!> This module provides functions and subroutines for solving linear systems,
MODULE NAFPack_linear_system

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    USE NAFPack_matricielle
    USE NAFPack_Direct_methode
    USE NAFPack_Iterative_methods
    USE NAFPack_Eigen

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: Direct_methode, Iterative_methods

    CONTAINS

!################## direct methode ######################################################

    !> Direct method for solving linear systems
    !> \[ A * x = b \]
    !> This function allows you to choose the direct method for solving the linear system.
    !> It supports various methods:
    !> 
    !> - Gaussian elimination
    !> - LU decomposition
    !> - LDU decomposition
    !> - Cholesky decomposition
    !> - Alternative Cholesky decomposition
    !> - QR decomposition (Householder, Givens, Classical Gram-Schmidt, Modified Gram-Schmidt)
    !> - TDMA (Thomas algorithm)
    FUNCTION Direct_methode(A, b, method, pivot_method, check) RESULT(x)

        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: pivot_method
        REAL(dp),DIMENSION(: ,:), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        LOGICAL, OPTIONAL, INTENT(IN) :: check
        REAL(dp),DIMENSION(SIZE(A,1)) :: x
        REAL(dp),DIMENSION(SIZE(A,1),SIZE(A,2)) :: Ainv
        REAL(dp),DIMENSION(SIZE(A,1)+1) :: c
        INTEGER :: N
        LOGICAL :: do_check = .TRUE., success

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
        ELSE IF(method == "A_LDL_Cholesky")THEN
            x = A_LDL_Cholesky(A, b, check = do_check)
        ELSE IF (INDEX(method, "QR") == 1) THEN
            x = A_QR(A, b, method = method)
        ELSE IF(method == "TDMA")THEN
            x = TDMA(A, b, check = do_check)
        ELSE IF(method == "Faddeev_Leverrier") THEN
            CALL Faddeev_Leverrier(A, c, Ainv = Ainv, success = success, check = do_check)
            IF (.NOT. success) THEN
                PRINT*, "WARNING :: Faddeev-Leverrier method failed, using LU decomposition instead"
                x = A_LU(A, b)
            ELSE
                x = MATMUL(Ainv, b)
            END IF
        ELSE
            STOP "ERROR : Wrong method for linear system (direct_methode)"
        END IF

    END FUNCTION Direct_methode


!################## Iterative methods ###################################################

    !> Iterative method for solving linear systems
    !> \[ A * x = b \]
    !> This function allows you to choose the iterative method for solving the linear system.
    !> It supports various methods:
    !> 
    !> - Jacobi 
    !> - Gauss-Seidel
    !> - Successive Over-Relaxation (SOR)
    !> - strongly implicit procedure (SIP_ILU, SIP_ICF)
    !> - Symmetric Successive Over-Relaxation (SSOR)
    FUNCTION Iterative_methods(A, b, method, x_init, max_iter, omega) RESULT(x)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: x_init
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        INTEGER, OPTIONAL, INTENT(IN) :: max_iter
        REAL(dp), OPTIONAL, INTENT(IN) :: omega
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: L, U
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

        IF (method == "SIP_ILU")THEN
            residu = x0
            CALL ILU_decomposition(A, L, U)
        ELSE IF (method == "SIP_ICF")THEN
            residu = x0
            CALL Incomplete_Cholesky_decomposition(A, L)
        END IF

        DO k = 1, max_iter_choice

            IF(k == kmax) THEN
                PRINT*, "WARNING :: non-convergence of the iterative method "//method
                EXIT
            END IF

            IF(method == "Jacobi")THEN
                CALL Jacobi(A, b , x0, x_new)
            ELSE IF(method == "Gauss_Seidel")THEN
                CALL Gauss_Seidel(A, b, x0, x_new)
            ELSE IF(method == "SOR")THEN
                IF (PRESENT(omega))THEN
                    IF(omega <= 0.d0 .OR. omega >= 2.d0) THEN
                        PRINT*, "WARNING :: omega must be in (0, 2), using default value of 1"
                        CALL SOR(A, b, x0, x_new, 1.d0)
                    ELSE 
                        CALL SOR(A, b, x0, x_new, omega)
                    END IF
                ELSE
                    CALL SOR(A, b, x0, x_new, 1.d0)
                END IF
            ELSE IF(method == "SIP_ILU") THEN
                IF (PRESENT(omega))THEN
                    IF(omega <= 0.d0 .OR. omega >= 2.d0) THEN
                        PRINT*, "WARNING :: omega must be in (0, 2), using default value of 1"
                        CALL SIP_ILU(L, U, x0, x_new, residu, 1.d0)
                    ELSE 
                        CALL SIP_ILU(L, U, x0, x_new, residu, omega)
                    END IF
                ELSE
                    CALL SIP_ILU(L, U, x0, x_new, residu, 1.d0)
                END IF
            ELSE IF(method == "SIP_ICF") THEN
                IF (PRESENT(omega))THEN
                    IF(omega <= 0.d0 .OR. omega >= 2.d0) THEN
                        PRINT*, "WARNING :: omega must be in (0, 2), using default value of 1"
                        CALL SIP_ICF(L, x0, x_new, residu, 1.d0)
                    ELSE 
                        CALL SIP_ICF(L, x0, x_new, residu, omega)
                    END IF
                ELSE
                    CALL SIP_ICF(L, x0, x_new, residu, 1.d0)
                END IF
            ELSE IF(method == "SSOR")THEN
                IF (PRESENT(omega))THEN
                    IF(omega <= 0.d0 .OR. omega >= 2.d0) THEN
                        PRINT*, "WARNING :: omega must be in (0, 2), using default value of 1"
                        CALL SSOR(A, b, x0, x_new, 1.d0)
                    ELSE 
                        CALL SSOR(A, b, x0, x_new, omega)
                    END IF
                ELSE
                    CALL SSOR(A, b, x0, x_new, 1.d0)
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

END MODULE NAFPack_linear_system