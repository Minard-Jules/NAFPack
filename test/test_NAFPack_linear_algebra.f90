MODULE test_NAFPack_linear_algebra

    USE NAFPack_linalg
    USE NAFPack_constant
    USE NAFPack_Eigen
    USE NAFPack_Logger_mod

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE
    PUBLIC :: test_linear_algebra
    PUBLIC :: test_linear_system

CONTAINS

!================ Linear System ===========================================================

    SUBROUTINE test_direct_method(x_true, x, method, stat, pivot_method)
        REAL(dp), DIMENSION(:), INTENT(IN) :: x, x_true
        CHARACTER(LEN=*), INTENT(IN) :: method
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: pivot_method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(size(x)) :: diff_x

        diff_x = x_true - x
        IF (maxval(abs(diff_x)) < epsi_test) THEN
            WRITE (*, '(A,T50,A,A)') green_color//method//" "//pivot_method, " :: OK"//reset_color
        ELSE
            WRITE (*, '(A,T50,A)') red_color//method//" "//pivot_method, " :: ECHEC"//reset_color
            stat = .TRUE.
        END IF

    END SUBROUTINE test_direct_method

    SUBROUTINE test_iterative_method(x_true, x, method, stat, preconditioner)
        REAL(dp), DIMENSION(:), INTENT(IN) :: x, x_true
        CHARACTER(LEN=*), INTENT(IN) :: method
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: preconditioner
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(size(x)) :: diff_x

        diff_x = x_true - x
        IF (maxval(abs(diff_x)) < epsi_test) THEN
            WRITE (*, '(A,T50,A,A)') green_color//method//" "//preconditioner, " :: OK"//reset_color
        ELSE
            WRITE (*, '(A,T50,A)') red_color//method//" "//preconditioner, " :: ECHEC"//reset_color
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

        A = reshape([4, -1, 0, &
                     -1, 4, -1, &
                     0, -1, 4], [N, N])
        b = [2, 4, 10]
        x = [1, 2, 3]

        WRITE (*, '(A, A)') cyan_color//"Direct methode"//reset_color
        PRINT*," "

        !================ Direct Methods =====================================================

        !=====================================================================================
        ! Gauss method
        CALL solver%direct%set_method(METHOD_Gauss)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Gauss", stat=stat)

        !=====================================================================================
        ! Gauss method with pivot partial
        CALL solver%direct%set_method(METHOD_Gauss, set_pivot_partial=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Gauss with partial pivoting", stat=stat)

        !=====================================================================================
        ! Gauss method with pivot partial
        CALL solver%direct%set_method(METHOD_Gauss, set_pivot_total=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Gauss with total pivoting", stat=stat)

        !=====================================================================================
        ! Gauss-Jordan method
        CALL solver%direct%set_method(METHOD_Gauss_JORDAN)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Gauss-Jordan", stat=stat)

        !=====================================================================================
        ! Gauss-Jordan method with pivot partial
        CALL solver%direct%set_method(METHOD_Gauss_JORDAN, set_pivot_partial=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), &
                                method="Gauss-Jordan with partial pivoting", stat=stat)

        !=====================================================================================
        ! Gauss-Jordan method with pivot total
        CALL solver%direct%set_method(METHOD_Gauss_JORDAN, set_pivot_total=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), &
                                method="Gauss-Jordan with total pivoting", stat=stat)

        !=====================================================================================
        ! LU decomposition method
        CALL solver%direct%set_method(METHOD_LU)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LU", stat=stat)

        !=====================================================================================
        ! LU decomposition method with pivot partial
        CALL solver%direct%set_method(METHOD_LU, set_pivot_partial=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LU with partial pivoting", stat=stat)

        !=====================================================================================
        ! LU decomposition method with pivot total
        CALL solver%direct%set_method(METHOD_LU, set_pivot_total=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LU with total pivoting", stat=stat)

        !=====================================================================================
        ! LDU decomposition method
        CALL solver%direct%set_method(METHOD_LDU)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LDU", stat=stat)

        !=====================================================================================
        ! LDU decomposition method with pivot partial
        CALL solver%direct%set_method(METHOD_LDU, set_pivot_partial=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LDU with partial pivoting", stat=stat)

        !=====================================================================================
        ! LDU decomposition method with pivot total
        CALL solver%direct%set_method(METHOD_LDU, set_pivot_total=.TRUE.)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LDU with total pivoting", stat=stat)

        !=====================================================================================
        ! Cholesky decomposition method
        CALL solver%direct%set_method(METHOD_Cholesky)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Cholesky", stat=stat)

        !=====================================================================================
        ! LDL decomposition method
        CALL solver%direct%set_method(METHOD_LDL_Cholesky)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LDL", stat=stat)

        !=====================================================================================
        ! QR decomposition method Givens
        CALL solver%direct%set_method(METHOD_QR)
        CALL solver%direct%set_qr_method(QR_GIVENS)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Givens", stat=stat)

        !=====================================================================================
        ! QR decomposition method Householder
        CALL solver%direct%set_method(METHOD_QR)
        CALL solver%direct%set_qr_method(QR_HOUSEHOLDER)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Householder", stat=stat)

        !=====================================================================================
        ! QR decomposition method Gram-Schmidt Classical
        CALL solver%direct%set_method(METHOD_QR)
        CALL solver%direct%set_qr_method(QR_GRAM_SCHMIDT)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Gram-Schmidt", stat=stat)

        !=====================================================================================
        ! QR decomposition method Gram-Schmidt Modified
        CALL solver%direct%set_method(METHOD_QR)
        CALL solver%direct%set_qr_method(QR_GRAM_SCHMIDT_Modified)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Gram-Schmidt Modified", stat=stat)

        !=====================================================================================
        ! TDMA method TDMA
        CALL solver%direct%set_method(METHOD_TDMA)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Gram-Schmidt Modified", stat=stat)

        !=====================================================================================
        ! Faddeev leverrier Le method
        CALL solver%direct%set_method(METHOD_FADDEEV_LEVERRIER)
        CALL solver%direct%test_matrix(A)
        CALL test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Faddeev-Le Verrier", stat=stat)

        PRINT*," "
        PRINT*," "
        WRITE (*, '(A, A)') cyan_color//"Iterative_methods"//reset_color
        PRINT*," "

        !================ Iterative Methods ==================================================

        verbose%to_terminal = .FALSE.
        verbose%to_file = .TRUE.
        verbose%show_final = .TRUE.
        verbose%show_iteration = .TRUE.
        verbose%verbosity_level = 4
        verbose%frequency = 1
        verbose%file_format = FORMAT_FILE_LOG

        !=====================================================================================
        !Jacobi method

        verbose%filename = "Log/Log_Jacobi"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_Jacobi)
        params = solver%iterative%Init_IterativeParams(N)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose), &
                                   method="Jacobi", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Gauss-Seidel method

        verbose%filename = "Log/Log_Gauss_Seidel"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_Gauss_Seidel)
        params = solver%iterative%Init_IterativeParams(N)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Gauss-Seidel", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !SOR method

        verbose%filename = "Log/Log_SOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_SOR)
        params = solver%iterative%Init_IterativeParams(N, omega=1.d0)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="SOR", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !JOR method

        verbose%filename = "Log/Log_JOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_JOR)
        params = solver%iterative%Init_IterativeParams(N, omega=1.d0)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="JOR", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !SIP_ILU method

        verbose%filename = "Log/Log_SIP_ILU"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_SIP_ILU)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="SIP ILU", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !SIP_ICF method

        verbose%filename = "Log/Log_SIP_ICF"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_SIP_ICF)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="SIP ICF", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !SSOR method

        verbose%filename = "Log/Log_SSOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_SSOR)
        params = solver%iterative%Init_IterativeParams(N, omega=1.d0)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="SSOR", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson method

        verbose%filename = "Log/Log_Richardson"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, alpha=0.2d0)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson method with Jacobi preconditioner

        verbose%filename = "Log/Log_Richardson_Jacobi"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, alpha=1.d0, method_preconditioner=METHOD_PRECOND_JACOBI)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="Jacobi preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson method with Gauss-Seidel preconditioner

        verbose%filename = "Log/Log_Richardson_GS"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, alpha=1.d0, method_preconditioner=METHOD_PRECOND_GS)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="Gauss-Seidel preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson method with SOR preconditioner

        verbose%filename = "Log/Log_Richardson_SOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="SOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson method with JOR preconditioner

        verbose%filename = "Log/Log_Richardson_JOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="JOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson method with ILU preconditioner

        verbose%filename = "Log/Log_Richardson_ILU"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="ILU preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson method with ICF preconditioner

        verbose%filename = "Log/Log_Richardson_ICF"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="ICF preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson method with SSOR preconditioner

        verbose%filename = "Log/Log_Richardson_SSOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SSOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="SSOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson unsteady method

        verbose%filename = "Log/Log_Richardson_unsteady"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, is_stationary=.FALSE.)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson unsteady method with Jacobi preconditioner

        verbose%filename = "Log/Log_Richardson_unsteady_Jacobi"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI, &
                                                       is_stationary=.FALSE.)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="Jacobi preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson unsteady method with Gauss-Seidel preconditioner

        verbose%filename = "Log/Log_Richardson_unsteady_GS"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_GS, is_stationary=.FALSE.)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="Gauss-Seidel preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson unsteady method with SOR preconditioner

        verbose%filename = "Log/Log_Richardson_unsteady_SOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_SOR, is_stationary=.FALSE.)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="SOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson unsteady method with JOR preconditioner

        verbose%filename = "Log/Log_Richardson_unsteady_JOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_JOR, is_stationary=.FALSE.)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="JOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson unsteady method with ILU preconditioner

        verbose%filename = "Log/Log_Richardson_unsteady_ILU"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_ILU, is_stationary=.FALSE.)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="ILU preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson unsteady method with ICF preconditioner

        verbose%filename = "Log/Log_Richardson_unsteady_ICF"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_ICF, is_stationary=.FALSE.)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="ICF preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Richardson unsteady method with SSOR preconditioner

        verbose%filename = "Log/Log_Richardson_unsteady_SSOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_SSOR, is_stationary=.FALSE.)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="SSOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Gradient method

        verbose%filename = "Log/Log_Conjugate_Gradient"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Gradient method with Jacobi preconditioner

        verbose%filename = "Log/Log_Conjugate_Gradient_Jacobi"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="Jacobi preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Gradient method with JOR preconditioner

        verbose%filename = "Log/Log_Conjugate_Gradient_JOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="JOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Gradient method with ILU preconditioner

        verbose%filename = "Log/Log_Conjugate_Gradient_ILU"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="ILU preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Gradient method with ICF preconditioner

        verbose%filename = "Log/Log_Conjugate_Gradient_ICF"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="ICF preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Gradient method with SSOR preconditioner

        verbose%filename = "Log/Log_Conjugate_Gradient_SSOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SSOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="SSOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Residual method

        verbose%filename = "Log/Log_Conjugate_Residual"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Residual method with Jacobi preconditioner

        verbose%filename = "Log/Log_Conjugate_Residual_Jacobi"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="Jacobi preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Residual method with JOR preconditioner

        verbose%filename = "Log/Log_Conjugate_Residual_JOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="JOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Residual method with ILU preconditioner

        verbose%filename = "Log/Log_Conjugate_Residual_ILU"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="ILU preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Residual method with ICF preconditioner

        verbose%filename = "Log/Log_Conjugate_Residual_ICF"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="ICF preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !Conjugate Residual method with SSOR preconditioner

        verbose%filename = "Log/Log_Conjugate_Residual_SSOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.0d0, method_preconditioner=METHOD_PRECOND_SSOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="SSOR preconditioner")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNR method

        verbose%filename = "Log/Log_CGNR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNR method with Jacobi preconditioner

        verbose%filename = "Log/Log_CGNR_Jacobi"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="Jacobi")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNR method with JOR preconditioner

        verbose%filename = "Log/Log_CGNR_JOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="JOR")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNR method with ILU preconditioner

        verbose%filename = "Log/Log_CGNR_ILU"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="ILU")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNR method with ICF preconditioner

        verbose%filename = "Log/Log_CGNR_ICF"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="ICF")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNR method with SSOR preconditioner

        verbose%filename = "Log/Log_CGNR_SSOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SSOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="SSOR")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNE method

        verbose%filename = "Log/Log_CGNE"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat)
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNE method with Jacobi preconditioner

        verbose%filename = "Log/Log_CGNE_Jacobi"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="Jacobi")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNE method with JOR preconditioner

        verbose%filename = "Log/Log_CGNE_JOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="JOR")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNE method with ILU preconditioner

        verbose%filename = "Log/Log_CGNE_ILU"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="ILU")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNE method with ICF preconditioner

        verbose%filename = "Log/Log_CGNE_ICF"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="ICF")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================
        !CGNE method with SSOR preconditioner

        verbose%filename = "Log/Log_CGNE_SSOR"
        CALL verbose%init()

        CALL solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SSOR)
        CALL solver%iterative%test_matrix(A, params, verbose)
        CALL test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="SSOR")
        CALL solver%iterative%Dealocate_IterativeParams(params)

        CALL verbose%CLOSE()

        !=====================================================================================

    END SUBROUTINE test_linear_system

!================ Eigen ===================================================================

    SUBROUTINE test_Eigen_vecteur_propre(A, method, stat)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN=*), INTENT(IN) :: method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(size(A, 1), size(A, 2)) :: vp
        REAL(dp), DIMENSION(size(A, 1)) :: lambda, diff
        INTEGER :: i
        LOGICAL :: verif_Eigen = .TRUE.

        WRITE (*, '(A, A)') cyan_color//"Eigen value"//reset_color
        PRINT*," "

        CALL Eigen(A, lambda, vp=vp, method=method)

        DO i = 1, size(A, 1)
            diff = matmul(A, vp(i, :)) - lambda(i) * vp(i, :)
            IF (maxval(abs(diff)) > epsi_test) verif_Eigen = .FALSE.
        END DO

        IF (verif_Eigen) THEN
            WRITE (*, '(A,T50,A,A)') green_color//method, " :: OK"//reset_color
        ELSE
            WRITE (*, '(A,T50,A)') red_color//method, " :: ECHEC"//reset_color
            stat = .TRUE.
        END IF

    END SUBROUTINE test_Eigen_vecteur_propre

    SUBROUTINE test_Eigen_without_vecteur_propre(A, method, stat)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN=*), INTENT(IN) :: method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(size(A, 1)) :: lambda_test, lambda, diff
        CHARACTER(LEN=25) :: method_test
        LOGICAL :: verif_Eigen = .TRUE.

        method_test = "Power_iteration"

        IF (stat) THEN
            verif_Eigen = .FALSE.
        ELSE
            lambda_test = [-4, 2, 1]
            CALL Eigen(A, lambda, method=method)

            diff = lambda_test - lambda
            IF (maxval(abs(diff)) > epsi_test) verif_Eigen = .FALSE.
        END IF

        IF (verif_Eigen) THEN
            WRITE (*, '(A,T50,A,A)') green_color//method, " :: OK"//reset_color
        ELSE
            WRITE (*, '(A,T50,A)') red_color//method, " :: ECHEC"//reset_color
            stat = .TRUE.
        END IF

    END SUBROUTINE test_Eigen_without_vecteur_propre

    SUBROUTINE test_linear_algebra(stat)

        LOGICAL, INTENT(INOUT) :: stat
        LOGICAL :: stat_tmp = .FALSE.

        INTEGER, PARAMETER :: N = 3
        REAL(dp), DIMENSION(N, N) :: A

        A = reshape([0, 3, -2, 2, -2, 2, -1, 0, 1], [N, N])
        CALL test_Eigen_vecteur_propre(A, "Power_iteration", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Householder", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Givens", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Classical", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Modified", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Householder_Shifted", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Givens_Shifted", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Classical_Shifted", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Modified_Shifted", stat_tmp)
        IF (stat_tmp) stat = stat_tmp

    END SUBROUTINE test_linear_algebra

END MODULE test_NAFPack_linear_algebra
