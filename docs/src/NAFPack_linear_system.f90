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
    USE NAFPack_Logger_mod
    USE NAFPack_Preconditioners

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

        REAL(dp),DIMENSION(: ,:), INTENT(IN) :: A
        REAL(dp),DIMENSION(:), INTENT(IN) :: b
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: pivot_method
        LOGICAL, OPTIONAL, INTENT(IN) :: check

        REAL(dp), DIMENSION(SIZE(A,1)) :: x
        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Ainv
        REAL(dp), DIMENSION(:), ALLOCATABLE :: c
        INTEGER :: N, alloc_stat
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

            ALLOCATE(Ainv(SIZE(A,1),SIZE(A,2)), c(SIZE(A,1)+1), stat=alloc_stat)
            IF (alloc_stat /= 0) STOP "ERROR :: Memory allocation failed for Ainv or c in Faddeev_Leverrier"

            CALL Faddeev_Leverrier(A, c, Ainv = Ainv, success = success, check = do_check)
            IF (.NOT. success) THEN
                PRINT*, "WARNING :: Faddeev-Leverrier method failed, using LU decomposition instead"
                x = A_LU(A, b)
            ELSE
                x = MATMUL(Ainv, b)
            END IF

            DEALLOCATE(Ainv, c)
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
    FUNCTION Iterative_methods(A, b, method, omega, x_init, max_iter, tol, verbose) RESULT(x)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        REAL(dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: x_init
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: method
        INTEGER, OPTIONAL, INTENT(IN) :: max_iter
        REAL(dp), OPTIONAL, INTENT(IN) :: omega, tol
        TYPE(Logger), OPTIONAL :: verbose
        REAL(dp), DIMENSION(:, :), ALLOCATABLE :: L, U
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x0, x_new, residu
        REAL(dp) :: epsi_tol
        INTEGER :: k, max_iter_choice, N, alloc_stat
        CHARACTER(LEN=64) :: msg

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
        
        IF (PRESENT(tol)) THEN
            epsi_tol = tol
        ELSE
            epsi_tol = epsi
        END IF

        IF (method == "SIP_ILU")THEN
            ALLOCATE(L(SIZE(A,1),SIZE(A,2)), U(SIZE(A,1),SIZE(A,2)), stat=alloc_stat)
            IF (alloc_stat /= 0) STOP "ERROR :: Memory allocation failed for L or U in SIP_ILU"
            residu = b - MATMUL(A, x0)
            CALL ILU_decomposition(A, L, U)
        ELSE IF (method == "SIP_ICF")THEN
            ALLOCATE(L(SIZE(A,1),SIZE(A,2)), stat=alloc_stat)
            IF (alloc_stat /= 0) STOP "ERROR :: Memory allocation failed for L or U in SIP_ICF"
            residu = b - MATMUL(A, x0)
            CALL Incomplete_Cholesky_decomposition(A, L)
        ELSE IF (INDEX(method, "Richardson") == 1)THEN
            residu = b - MATMUL(A, x0)
        END IF

        DO k = 1, max_iter_choice

            IF(k == max_iter_choice) THEN
                PRINT*, "WARNING :: non-convergence of the iterative method "//method
                PRINT*, "Residual norm: ", NORM2(residu)
                EXIT
            END IF

            IF(method == "Jacobi")THEN
                CALL Jacobi(A, b , x0, x_new)
            ELSE IF(method == "Gauss_Seidel")THEN
                CALL Gauss_Seidel(A, b, x0, x_new)
            ELSE IF(method == "SOR")THEN
                CALL SOR(A, b, x0, x_new, Get_Omega(method, omega))
            ELSE IF(method == "JOR")THEN
                CALL JOR(A, b, x0, x_new, Get_Omega(method, omega))
            ELSE IF(method == "SIP_ILU") THEN
                CALL SIP_ILU(L, U, x0, x_new, residu, Get_Omega(method, omega))
            ELSE IF(method == "SIP_ICF") THEN
                CALL SIP_ICF(L, x0, x_new, residu, Get_Omega(method, omega))
            ELSE IF(method == "SSOR")THEN
                CALL SSOR(A, b, x0, x_new, Get_Omega(method, omega))
            ELSE IF(method == "Richardson_Stationary")THEN
                CALL Richardson(x0, x_new, residu, Get_Omega(method, omega), method)
            ELSE IF(method == "Richardson_ILU_Preconditioned")THEN
                CALL Richardson(x0, x_new, residu, Get_Omega(method, omega), method, preconditioner_choice=ILU_preconditioner, A=A)
            ELSE
                STOP "ERROR : Wrong method for linear system (Iterative_methods)"
            END IF

            residu = b - MATMUL(A, x_new)

            IF (PRESENT(verbose)) THEN
                verbose%step = k
                WRITE(msg, '(A,I3,A,ES14.7)') "Iter ", k, " | Norm residu: ", NORM2(residu)
                CALL verbose%log(msg)
            END IF

            IF (NORM2(residu) < epsi_tol) EXIT

            x0 = x_new

        END DO

        IF (method == "SIP_ILU")THEN
            DEALLOCATE(L, U)
        ELSE IF (method == "SIP_ICF")THEN
            DEALLOCATE(L)
        ELSE IF (INDEX(method, "Richardson") == 1) THEN
            preconditioner_use = preconditioner(.FALSE.)
            IF (method == "Richardson_ILU_Preconditioned") CALL preconditioner_use%DEALLOCATE_ILU()
        END IF

        x = x_new

    END FUNCTION Iterative_methods

    FUNCTION Get_Omega(method, omega) RESULT(w)
        CHARACTER(LEN = *), INTENT(IN) :: method
        REAL(dp), OPTIONAL, INTENT(IN) :: omega
        REAL(dp) :: w
        
        IF (PRESENT(omega)) THEN
            IF (INDEX(method, "SOR") > 0 .AND. omega > 0.d0 .AND. omega < 2.d0) THEN
                w = omega
            ELSE IF (method == "JOR" .AND. omega > 0.d0 .AND. omega <= 1.d0) THEN
                w = omega
            ELSE IF (INDEX(method, "SIP") == 1 .AND. omega > 0.d0)THEN
                w = omega
            ELSE IF (INDEX(method, "Richardson") == 1 .AND. omega > 0.d0) THEN
                w = omega
            ELSE
                IF (INDEX(method, "SOR") > 0) THEN
                    PRINT*, "WARNING :: In method "//method//", omega must be in (0, 2), using default value of 1"
                ELSE IF (method == "JOR")THEN
                    PRINT*, "WARNING :: In method "//method//", omega must be in (0, 1), using default value of 1"
                ELSE
                    PRINT*, "WARNING :: In method "//method//", omega must be > 0, using default value of 1"
                END IF
                w = 1.d0
            END IF
        ELSE
            w = 1.d0
        END IF
    END FUNCTION

END MODULE NAFPack_linear_system