module test_NAFPack_linear_algebra

    use NAFPack_linalg
    use NAFPack_constant
    use NAFPack_Eigen
    use NAFPack_ANSI, only: ColorsAscii
    use NAFPack_Logger_mod

    implicit none(type, external)

    private
    public :: test_linear_algebra
    public :: test_linear_system

contains

!================ Linear System ===========================================================

    subroutine test_direct_method(x_true, x, method, stat, pivot_method)
        real(dp), dimension(:), intent(in) :: x, x_true
        character(LEN=*), intent(in) :: method
        character(LEN=*), optional, intent(in) :: pivot_method
        logical, intent(inout) :: stat
        real(dp), dimension(size(x)) :: diff_x
        type(ColorsAscii) :: colors

        CALL colors%init()

        diff_x = x_true - x
        if (maxval(abs(diff_x)) < TOL_CONVERGENCE_dp) then
            write (*, '(A,T50,A,A)') colors%green//method//" "//pivot_method, " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//method//" "//pivot_method, " :: ECHEC"//colors%reset
            stat = .true.
        end if

    end subroutine test_direct_method

    subroutine test_iterative_method(x_true, x, method, stat, preconditioner)
        real(dp), dimension(:), intent(in) :: x, x_true
        character(LEN=*), intent(in) :: method
        character(LEN=*), optional, intent(in) :: preconditioner
        logical, intent(inout) :: stat
        real(dp), dimension(size(x)) :: diff_x
        type(ColorsAscii) :: colors

        CALL colors%init()

        diff_x = x_true - x
        if (maxval(abs(diff_x)) < TOL_CONVERGENCE_dp) then
            write (*, '(A,T50,A,A)') colors%green//method//" "//preconditioner, " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//method//" "//preconditioner, " :: ECHEC"//colors%reset
            stat = .true.
        end if

    end subroutine test_iterative_method

    subroutine test_linear_system(stat)

        logical, intent(inout) :: stat

        integer, parameter :: N = 3

        real(dp), dimension(N, N) :: A
        real(dp), dimension(N) :: b
        real(dp), dimension(N) :: x

        type(linalg) :: solver
        type(IterativeParams) :: params
        type(Logger) :: verbose

        character(len=100) :: fname
        logical :: exists_dir
        type(ColorsAscii) :: colors

        CALL colors%init()

        A = reshape([4, -1, 0, &
                     -1, 4, -1, &
                     0, -1, 4], [N, N])
        b = [2, 4, 10]
        x = [1, 2, 3]

        write (*, '(A, A)') colors%cyan//"Direct methode"//colors%reset
        print*," "

        !================ Direct Methods =====================================================

        !=====================================================================================
        ! Gauss method
        call solver%direct%set_method(METHOD_Gauss)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Gauss", stat=stat)

        !=====================================================================================
        ! Gauss method with pivot partial
        call solver%direct%set_method(METHOD_Gauss, set_pivot_partial=.true.)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Gauss with partial pivoting", stat=stat)

        !=====================================================================================
        ! Gauss method with pivot partial
        call solver%direct%set_method(METHOD_Gauss, set_pivot_total=.true.)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Gauss with total pivoting", stat=stat)

        !=====================================================================================
        ! Gauss-Jordan method
        call solver%direct%set_method(METHOD_Gauss_JORDAN)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Gauss-Jordan", stat=stat)

        !=====================================================================================
        ! Gauss-Jordan method with pivot partial
        call solver%direct%set_method(METHOD_Gauss_JORDAN, set_pivot_partial=.true.)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), &
                                method="Gauss-Jordan with partial pivoting", stat=stat)

        !=====================================================================================
        ! Gauss-Jordan method with pivot total
        call solver%direct%set_method(METHOD_Gauss_JORDAN, set_pivot_total=.true.)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), &
                                method="Gauss-Jordan with total pivoting", stat=stat)

        !=====================================================================================
        ! LU decomposition method
        call solver%direct%set_method(METHOD_LU)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LU", stat=stat)

        !=====================================================================================
        ! LU decomposition method with pivot partial
        call solver%direct%set_method(METHOD_LU, set_pivot_partial=.true.)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LU with partial pivoting", stat=stat)

        !=====================================================================================
        ! LU decomposition method with pivot total
        call solver%direct%set_method(METHOD_LU, set_pivot_total=.true.)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LU with total pivoting", stat=stat)

        !=====================================================================================
        ! LDU decomposition method
        call solver%direct%set_method(METHOD_LDU)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LDU", stat=stat)

        !=====================================================================================
        ! LDU decomposition method with pivot partial
        call solver%direct%set_method(METHOD_LDU, set_pivot_partial=.true.)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LDU with partial pivoting", stat=stat)

        !=====================================================================================
        ! LDU decomposition method with pivot total
        call solver%direct%set_method(METHOD_LDU, set_pivot_total=.true.)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LDU with total pivoting", stat=stat)

        !=====================================================================================
        ! Cholesky decomposition method
        call solver%direct%set_method(METHOD_Cholesky)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Cholesky", stat=stat)

        !=====================================================================================
        ! LDL decomposition method
        call solver%direct%set_method(METHOD_LDL_Cholesky)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="LDL", stat=stat)

        !=====================================================================================
        ! QR decomposition method Givens
        call solver%direct%set_method(METHOD_QR)
        call solver%direct%set_qr_method(QR_GIVENS)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Givens", stat=stat)

        !=====================================================================================
        ! QR decomposition method Householder
        call solver%direct%set_method(METHOD_QR)
        call solver%direct%set_qr_method(QR_HOUSEHOLDER)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Householder", stat=stat)

        !=====================================================================================
        ! QR decomposition method Gram-Schmidt Classical
        call solver%direct%set_method(METHOD_QR)
        call solver%direct%set_qr_method(QR_GRAM_SCHMIDT)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Gram-Schmidt", stat=stat)

        !=====================================================================================
        ! QR decomposition method Gram-Schmidt Modified
        call solver%direct%set_method(METHOD_QR)
        call solver%direct%set_qr_method(QR_GRAM_SCHMIDT_Modified)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Gram-Schmidt Modified", stat=stat)

        !=====================================================================================
        ! TDMA method TDMA
        call solver%direct%set_method(METHOD_TDMA)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="QR Gram-Schmidt Modified", stat=stat)

        !=====================================================================================
        ! Faddeev leverrier Le method
        call solver%direct%set_method(METHOD_FADDEEV_LEVERRIER)
        call solver%direct%test_matrix(A)
        call test_direct_method(x_true=x, x=solver%direct%solve(A, b), method="Faddeev-Le Verrier", stat=stat)

        print*," "
        print*," "
        write (*, '(A, A)') colors%cyan//"Iterative_methods"//colors%reset
        print*," "

        !================ Iterative Methods ==================================================

        verbose%to_terminal = .false.
        verbose%to_file = .true.
        verbose%show_final = .true.
        verbose%show_iteration = .true.
        verbose%verbosity_level = 4
        verbose%frequency = 1
        verbose%file_format = FORMAT_FILE_LOG

        fname = "Log"

        inquire(FILE = fname, EXIST = exists_dir)

        if (.not. exists_dir) then
            call execute_command_line("mkdir "//trim(fname))
        end if

        !=====================================================================================
        !Jacobi method

        verbose%filename = trim(fname)//"/Log_Jacobi"
        call verbose%init()

        call solver%iterative%set_method(METHOD_Jacobi)
        params = solver%iterative%Init_IterativeParams(N)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose), &
                                   method="Jacobi", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Gauss-Seidel method

        verbose%filename = trim(fname)//"/Log_Gauss_Seidel"
        call verbose%init()

        call solver%iterative%set_method(METHOD_Gauss_Seidel)
        params = solver%iterative%Init_IterativeParams(N)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Gauss-Seidel", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !SOR method

        verbose%filename = trim(fname)//"/Log_SOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_SOR)
        params = solver%iterative%Init_IterativeParams(N, omega=1.d0)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="SOR", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !JOR method

        verbose%filename = trim(fname)//"/Log_JOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_JOR)
        params = solver%iterative%Init_IterativeParams(N, omega=1.d0)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="JOR", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !SIP_ILU method

        verbose%filename = trim(fname)//"/Log_SIP_ILU"
        call verbose%init()

        call solver%iterative%set_method(METHOD_SIP_ILU)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="SIP ILU", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !SIP_ICF method

        verbose%filename = trim(fname)//"/Log_SIP_ICF"
        call verbose%init()

        call solver%iterative%set_method(METHOD_SIP_ICF)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="SIP ICF", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !SSOR method

        verbose%filename = trim(fname)//"/Log_SSOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_SSOR)
        params = solver%iterative%Init_IterativeParams(N, omega=1.d0)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="SSOR", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson method

        verbose%filename = trim(fname)//"/Log_Richardson"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, alpha=0.2d0)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson method with Jacobi preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_Jacobi"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, alpha=1.d0, method_preconditioner=METHOD_PRECOND_JACOBI)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="Jacobi preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson method with Gauss-Seidel preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_GS"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, alpha=1.d0, method_preconditioner=METHOD_PRECOND_GS)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="Gauss-Seidel preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson method with SOR preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_SOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="SOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson method with JOR preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_JOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="JOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson method with ILU preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_ILU"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="ILU preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson method with ICF preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_ICF"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="ICF preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson method with SSOR preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_SSOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SSOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson", stat=stat, preconditioner="SSOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson unsteady method

        verbose%filename = trim(fname)//"/Log_Richardson_unsteady"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, is_stationary=.false.)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson unsteady method with Jacobi preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_unsteady_Jacobi"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI, &
                                                       is_stationary=.false.)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="Jacobi preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson unsteady method with Gauss-Seidel preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_unsteady_GS"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_GS, is_stationary=.false.)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="Gauss-Seidel preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson unsteady method with SOR preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_unsteady_SOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_SOR, is_stationary=.false.)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="SOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson unsteady method with JOR preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_unsteady_JOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_JOR, is_stationary=.false.)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="JOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson unsteady method with ILU preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_unsteady_ILU"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_ILU, is_stationary=.false.)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="ILU preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson unsteady method with ICF preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_unsteady_ICF"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_ICF, is_stationary=.false.)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="ICF preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Richardson unsteady method with SSOR preconditioner

        verbose%filename = trim(fname)//"/Log_Richardson_unsteady_SSOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_RICHARDSON)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, &
                                                       method_preconditioner=METHOD_PRECOND_SSOR, is_stationary=.false.)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Richardson unsteady", stat=stat, preconditioner="SSOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Gradient method

        verbose%filename = trim(fname)//"/Log_Conjugate_Gradient"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Gradient method with Jacobi preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Gradient_Jacobi"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="Jacobi preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Gradient method with JOR preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Gradient_JOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="JOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Gradient method with ILU preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Gradient_ILU"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="ILU preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Gradient method with ICF preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Gradient_ICF"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="ICF preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Gradient method with SSOR preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Gradient_SSOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_GRADIENT)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SSOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Gradient", stat=stat, preconditioner="SSOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Residual method

        verbose%filename = trim(fname)//"/Log_Conjugate_Residual"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Residual method with Jacobi preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Residual_Jacobi"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="Jacobi preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Residual method with JOR preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Residual_JOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="JOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Residual method with ILU preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Residual_ILU"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="ILU preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Residual method with ICF preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Residual_ICF"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="ICF preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !Conjugate Residual method with SSOR preconditioner

        verbose%filename = trim(fname)//"/Log_Conjugate_Residual_SSOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CONJUGATE_RESIDUAL)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.0d0, method_preconditioner=METHOD_PRECOND_SSOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="Conjugate Residual", stat=stat, preconditioner="SSOR preconditioner")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNR method

        verbose%filename = trim(fname)//"/Log_CGNR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNR method with Jacobi preconditioner

        verbose%filename = trim(fname)//"/Log_CGNR_Jacobi"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="Jacobi")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNR method with JOR preconditioner

        verbose%filename = trim(fname)//"/Log_CGNR_JOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="JOR")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNR method with ILU preconditioner

        verbose%filename = trim(fname)//"/Log_CGNR_ILU"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="ILU")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNR method with ICF preconditioner

        verbose%filename = trim(fname)//"/Log_CGNR_ICF"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="ICF")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNR method with SSOR preconditioner

        verbose%filename = trim(fname)//"/Log_CGNR_SSOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNR)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SSOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNR", stat=stat, preconditioner="SSOR")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNE method

        verbose%filename = trim(fname)//"/Log_CGNE"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat)
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNE method with Jacobi preconditioner

        verbose%filename = trim(fname)//"/Log_CGNE_Jacobi"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, method_preconditioner=METHOD_PRECOND_JACOBI)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="Jacobi")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNE method with JOR preconditioner

        verbose%filename = trim(fname)//"/Log_CGNE_JOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_JOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="JOR")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNE method with ILU preconditioner

        verbose%filename = trim(fname)//"/Log_CGNE_ILU"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ILU)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="ILU")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNE method with ICF preconditioner

        verbose%filename = trim(fname)//"/Log_CGNE_ICF"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_ICF)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="ICF")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================
        !CGNE method with SSOR preconditioner

        verbose%filename = trim(fname)//"/Log_CGNE_SSOR"
        call verbose%init()

        call solver%iterative%set_method(METHOD_CGNE)
        params = solver%iterative%Init_IterativeParams(N, A=A, omega=1.d0, method_preconditioner=METHOD_PRECOND_SSOR)
        call solver%iterative%test_matrix(A, params, verbose)
        call test_iterative_method(x_true=x, x=solver%iterative%solve(A, b, params, verbose=verbose), &
                                   method="CGNE", stat=stat, preconditioner="SSOR")
        call solver%iterative%Dealocate_IterativeParams(params)

        call verbose%close()

        !=====================================================================================

    end subroutine test_linear_system

!================ Eigen ===================================================================

    subroutine test_Eigen_vecteur_propre(A, method, stat)

        real(dp), dimension(:, :), intent(in) :: A
        character(LEN=*), intent(in) :: method
        logical, intent(inout) :: stat
        real(dp), dimension(size(A, 1), size(A, 2)) :: vp
        real(dp), dimension(size(A, 1)) :: lambda, diff
        integer :: i
        logical :: verif_Eigen = .true.
        type(ColorsAscii) :: colors

        CALL colors%init()

        write (*, '(A, A)') colors%cyan//"Eigen value"//colors%reset
        print*," "

        call Eigen(A, lambda, vp=vp, method=method)

        do i = 1, size(A, 1)
            diff = matmul(A, vp(i, :)) - lambda(i) * vp(i, :)
            if (maxval(abs(diff)) > TOL_CONVERGENCE_dp) verif_Eigen = .false.
        end do

        if (verif_Eigen) then
            write (*, '(A,T50,A,A)') colors%green//method, " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//method, " :: ECHEC"//colors%reset
            stat = .true.
        end if

    end subroutine test_Eigen_vecteur_propre

    subroutine test_Eigen_without_vecteur_propre(A, method, stat)

        real(dp), dimension(:, :), intent(in) :: A
        character(LEN=*), intent(in) :: method
        logical, intent(inout) :: stat
        real(dp), dimension(size(A, 1)) :: lambda_test, lambda, diff
        character(LEN=25) :: method_test
        logical :: verif_Eigen = .true.
        type(ColorsAscii) :: colors

        CALL colors%init()

        method_test = "Power_iteration"

        if (stat) then
            verif_Eigen = .false.
        else
            lambda_test = [-4, 2, 1]
            call Eigen(A, lambda, method=method)

            diff = lambda_test - lambda
            if (maxval(abs(diff)) > TOL_CONVERGENCE_dp) verif_Eigen = .false.
        end if

        if (verif_Eigen) then
            write (*, '(A,T50,A,A)') colors%green//method, " :: OK"//colors%reset
        else
            write (*, '(A,T50,A)') colors%red//method, " :: ECHEC"//colors%reset
            stat = .true.
        end if

    end subroutine test_Eigen_without_vecteur_propre

    subroutine test_linear_algebra(stat)

        logical, intent(inout) :: stat
        logical :: stat_tmp = .false.

        integer, parameter :: N = 3
        real(dp), dimension(N, N) :: A

        A = reshape([0, 3, -2, 2, -2, 2, -1, 0, 1], [N, N])
        call test_Eigen_vecteur_propre(A, "Power_iteration", stat_tmp)
        call test_Eigen_without_vecteur_propre(A, "QR_Householder", stat_tmp)
        call test_Eigen_without_vecteur_propre(A, "QR_Givens", stat_tmp)
        call test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Classical", stat_tmp)
        call test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Modified", stat_tmp)
        call test_Eigen_without_vecteur_propre(A, "QR_Householder_Shifted", stat_tmp)
        call test_Eigen_without_vecteur_propre(A, "QR_Givens_Shifted", stat_tmp)
        call test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Classical_Shifted", stat_tmp)
        call test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Modified_Shifted", stat_tmp)
        if (stat_tmp) stat = stat_tmp

    end subroutine test_linear_algebra

end module test_NAFPack_linear_algebra
