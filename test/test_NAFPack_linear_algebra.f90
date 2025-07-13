MODULE test_NAFPack_linear_algebra

    USE NAFPack_linalg
    USE NAFPack_constant
    USE NAFPack_Eigen
    USE NAFPack_Logger_mod

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: test_linear_algebra
    PUBLIC :: test_linear_system

    CONTAINS

!================ Linear System =========================================================== 

    SUBROUTINE test_direct_method(x_true, x, method, stat, pivot_method)
        REAL(dp), DIMENSION(:), INTENT(IN) :: x, x_true
        CHARACTER(LEN = *), INTENT(IN) :: method
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: pivot_method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(SIZE(x)) :: diff_x

        diff_x = x_true - x
        IF (MAXVAL(ABS(diff_x)) < epsi_test) THEN
            WRITE(*,'(A,T50,A,A)') green_color//method//" "//pivot_method, " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T50,A)') red_color//method//" "//pivot_method, " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF
    
    END SUBROUTINE test_direct_method

    SUBROUTINE test_iterative_method(x_true, x, method, stat, preconditioner)
        REAL(dp), DIMENSION(:), INTENT(IN) :: x, x_true
        CHARACTER(LEN = *), INTENT(IN) :: method
        CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: preconditioner
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(SIZE(x)) :: diff_x

        
        diff_x = x_true - x
        IF (MAXVAL(ABS(diff_x)) < epsi_test) THEN
            WRITE(*,'(A,T50,A,A)') green_color//method//" "//preconditioner, " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T50,A)') red_color//method//" "//preconditioner, " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF
    
    END SUBROUTINE test_iterative_method


    SUBROUTINE test_linear_system(stat)

        LOGICAL, INTENT(INOUT) :: stat
        
        INTEGER, PARAMETER :: N = 3

        REAL(dp), DIMENSION(N, N) :: A
        REAL(dp), DIMENSION(N) :: b
        REAL(dp), DIMENSION(N) :: x

        TYPE(linalg) :: solver
        TYPE(IterativeParams) :: params
        TYPE(Logger) :: verbose

        CALL verbose%init()
        verbose%to_terminal = .TRUE.
        verbose%show_final = .TRUE.
        verbose%show_iteration = .FALSE.
        ! verbose%frequency = 1

        A = RESHAPE([4, -1, 0, &
                    -1, 4, -1, &
                     0, -1, 4], [N, N])
        b = [2, 4, 10]
        x = [1, 2, 3]

        WRITE(*,'(A, A)') cyan_color//"Direct methode"// reset_color
        PRINT*, " "

        !================ Direct Methods =====================================================

        !=====================================================================================
        ! Gauss method
        CALL solver%direct%set_method(METHOD_Gauss)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "Gauss", stat = stat)

        !=====================================================================================
        ! Gauss method with pivot partial
        CALL solver%direct%set_method(METHOD_Gauss, set_pivot_partial=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "Gauss with partial pivoting", stat = stat)

        !=====================================================================================
        ! Gauss method with pivot partial
        CALL solver%direct%set_method(METHOD_Gauss, set_pivot_total=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "Gauss with total pivoting", stat = stat)

        !=====================================================================================
        ! LU decomposition method
        CALL solver%direct%set_method(METHOD_LU)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "LU", stat = stat)
        

        !=====================================================================================
        ! LU decomposition method with pivot partial
        CALL solver%direct%set_method(METHOD_LU, set_pivot_partial=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "LU with partial pivoting", stat = stat)

        !=====================================================================================
        ! LU decomposition method with pivot total
        CALL solver%direct%set_method(METHOD_LU, set_pivot_total=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "LU with total pivoting", stat = stat)

        !=====================================================================================
        ! LDU decomposition method
        CALL solver%direct%set_method(METHOD_LDU)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "LDU", stat = stat)

        !=====================================================================================
        ! LDU decomposition method with pivot partial
        CALL solver%direct%set_method(METHOD_LDU, set_pivot_partial=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "LDU with partial pivoting", stat = stat)

        !=====================================================================================
        ! LDU decomposition method with pivot total
        CALL solver%direct%set_method(METHOD_LDU, set_pivot_total=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "LDU with total pivoting", stat = stat)

        !=====================================================================================
        ! Cholesky decomposition method
        CALL solver%direct%set_method(METHOD_Cholesky)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "Cholesky", stat = stat)

        !=====================================================================================
        ! LDL decomposition method
        CALL solver%direct%set_method(METHOD_LDL_Cholesky)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "LDL", stat = stat)

        !=====================================================================================
        ! QR decomposition method Givens
        CALL solver%direct%set_method(METHOD_QR)
        CALL solver%direct%set_qr_method(QR_GIVENS)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "QR Givens", stat = stat)

        !=====================================================================================
        ! QR decomposition method Householder
        CALL solver%direct%set_method(METHOD_QR)
        CALL solver%direct%set_qr_method(QR_HOUSEHOLDER)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "QR Householder", stat = stat)

        !=====================================================================================
        ! QR decomposition method Gram-Schmidt Classical
        CALL solver%direct%set_method(METHOD_QR)
        CALL solver%direct%set_qr_method(QR_GRAM_SCHMIDT)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "QR Gram-Schmidt", stat = stat)

        !=====================================================================================
        ! QR decomposition method Gram-Schmidt Modified
        CALL solver%direct%set_method(METHOD_QR)
        CALL solver%direct%set_qr_method(QR_GRAM_SCHMIDT_Modified)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "QR Gram-Schmidt Modified", stat = stat)

        !=====================================================================================
        ! TDMA method TDMA
        CALL solver%direct%set_method(METHOD_TDMA)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "QR Gram-Schmidt Modified", stat = stat)

        !=====================================================================================
        ! Faddeev leverrier Le method
        CALL solver%direct%set_method(METHOD_FADDEEV_LEVERRIER)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true = x, x = solver%direct%solve(A, b), method = "Faddeev-Le Verrier", stat = stat)


        PRINT*, " "
        PRINT*, " "
        WRITE(*,'(A, A)') cyan_color//"Iterative_methods"// reset_color
        PRINT*, " "


        !================ Iterative Methods ==================================================

        !=====================================================================================
        !Jacobi method
        CALL solver%iterative%set_method(METHOD_Jacobi)
        params = solver%iterative%Init_IterativeParams(N)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "Jacobi", stat = stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !Gauss-Seidel method
        CALL solver%iterative%set_method(METHOD_Gauss_Seidel)
        params = solver%iterative%Init_IterativeParams(N)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
        method = "Gauss-Seidel", stat = stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !SOR method
        CALL solver%iterative%set_method(METHOD_SOR)
        params = solver%iterative%Init_IterativeParams(N, omega=1.1d0)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "SOR", stat = stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !JOR method
        CALL solver%iterative%set_method(METHOD_JOR)
        params = solver%iterative%Init_IterativeParams(N, omega=0.9d0)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "JOR", stat = stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !SIP_ILU method
        CALL solver%iterative%set_method(METHOD_SIP_ILU)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "SIP ILU", stat = stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !SIP_ICF method
        CALL solver%iterative%set_method(METHOD_SIP_ICF)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "SIP ICF", stat = stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !SSOR method
        CALL solver%iterative%set_method(METHOD_SSOR)
        params = solver%iterative%Init_IterativeParams(N, omega=0.9d0)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "SSOR", stat = stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !Richardson method
        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, alpha=0.2d0)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "Richardson", stat = stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !Richardson method with Jacobi preconditioner
        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, alpha=1.d0, method_preconditioner=METHOD_PRECOND_JACOBI)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "Richardson", stat = stat, preconditioner="Jacobi preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !Richardson method with Gauss-Seidel preconditioner
        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, alpha=1.d0, method_preconditioner=METHOD_PRECOND_GS)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "Richardson", stat = stat, preconditioner="Gauss-Seidel preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !Richardson method with SOR preconditioner
        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.1d0, method_preconditioner=METHOD_PRECOND_SOR)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "Richardson", stat = stat, preconditioner="SOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !Richardson method with JOR preconditioner
        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "Richardson", stat = stat, preconditioner="JOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !Richardson method with ILU preconditioner
        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "Richardson", stat = stat, preconditioner="ILU preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !Richardson method with ICF preconditioner
        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "Richardson", stat = stat, preconditioner="ICF preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================
        !Richardson method with SSOR preconditioner
        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SSOR)
        CALL test_iterative_method(x_true = x, x = solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method = "Richardson", stat = stat, preconditioner="SSOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        !=====================================================================================

        CALL verbose%close()

    END SUBROUTINE test_linear_system

!================ Eigen ===================================================================    

    SUBROUTINE test_Eigen_vecteur_propre(A, method, stat)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), INTENT(IN) :: method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION (SIZE(A, 1), SIZE(A, 2)) :: vp
        REAL(dp), DIMENSION (SIZE(A, 1)) :: lambda, diff
        INTEGER :: i
        LOGICAL :: verif_Eigen = .TRUE.

        WRITE(*,'(A, A)') cyan_color//"Eigen value"// reset_color
        PRINT*, " "

        CALL Eigen(A, lambda, vp = vp, method = method)

        DO i = 1, SIZE(A, 1)
            diff = MATMUL(A, vp(i, :)) - lambda(i) * vp(i, :)
            IF(MAXVAL(ABS(diff)) > epsi_test) verif_Eigen = .FALSE.
        END DO

        IF(verif_Eigen)THEN
            WRITE(*,'(A,T50,A,A)') green_color//method, " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T50,A)') red_color//method, " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

    END SUBROUTINE test_Eigen_vecteur_propre

    SUBROUTINE test_Eigen_without_vecteur_propre(A, method, stat)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), INTENT(IN) :: method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(SIZE(A, 1)) :: lambda_test, lambda, diff
        CHARACTER(LEN = 25) :: method_test
        LOGICAL :: verif_Eigen = .TRUE.

        method_test = "Power_iteration"


        IF(stat)THEN
            verif_Eigen = .FALSE.
        ELSE
            lambda_test = [-4, 2, 1]
            CALL Eigen(A, lambda, method = method)

            diff = lambda_test - lambda
            IF(MAXVAL(ABS(diff)) > epsi_test) verif_Eigen = .FALSE.
        END IF

        IF(verif_Eigen)THEN
            WRITE(*,'(A,T50,A,A)') green_color//method, " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T50,A)') red_color//method, " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF


    END SUBROUTINE test_Eigen_without_vecteur_propre

    SUBROUTINE test_linear_algebra(stat)

        LOGICAL, INTENT(INOUT) :: stat
        LOGICAL :: stat_tmp = .FALSE.

        INTEGER, PARAMETER :: N = 3
        REAL(dp), DIMENSION (N,N) :: A

        A=RESHAPE([0, 3, -2, 2, -2, 2, -1, 0, 1], [N, N])
        CALL test_Eigen_vecteur_propre(A, "Power_iteration", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Householder", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Givens", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Classical", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Modified", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Householder_Shifted", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Givens_Shifted", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Classical_Shifted", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Modified_Shifted", stat_tmp)
        IF(stat_tmp) stat = stat_tmp

    END SUBROUTINE test_linear_algebra

END MODULE test_NAFPack_linear_algebra