!> Module for iterative methods in NAFPack
MODULE NAFPack_Iterative_methods

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    USE NAFPack_matricielle
    USE NAFPack_Iterative_types
    USE NAFPack_Logger_mod
    USE NAFPack_Preconditioners
    USE NAFPack_Iterative_Params

    IMPLICIT NONE
    
    PRIVATE

    PUBLIC :: IterativeMethod
    PUBLIC :: IterativeParams

    PUBLIC :: MethodTypeIterative
    PUBLIC :: METHOD_ITERATIVE_NONE
    PUBLIC :: METHOD_Jacobi, METHOD_JOR
    PUBLIC :: METHOD_GAUSS_SEIDEL, METHOD_SOR, METHOD_SSOR
    PUBLIC :: METHOD_SIP_ILU, METHOD_SIP_ICF
    PUBLIC :: METHOD_RICHARDSON

    PUBLIC :: MethodPreconditioner
    PUBLIC :: METHOD_PRECOND_NONE
    PUBLIC :: METHOD_PRECOND_JACOBI, METHOD_PRECOND_JOR
    PUBLIC :: METHOD_PRECOND_GS, METHOD_PRECOND_SOR, METHOD_PRECOND_SSOR
    PUBLIC :: METHOD_PRECOND_ILU, METHOD_PRECOND_ICF

    TYPE :: IterativeMethod
        PRIVATE
        TYPE(MethodTypeIterative) :: method_type = METHOD_ITERATIVE_NONE
        TYPE(MethodPreconditioner) :: preconditioner_type = METHOD_PRECOND_NONE
        PROCEDURE(solve_interface_Iterative), PASS(this), POINTER :: solve_method => NULL()

        CONTAINS

        PROCEDURE :: set_method => set_method
        PROCEDURE :: solve => IterativeMethod_solve
        PROCEDURE :: Init_IterativeParams => Init_IterativeParams
        PROCEDURE :: Dealocate_IterativeParams => Dealocate_IterativeParams

    END TYPE IterativeMethod

    ABSTRACT INTERFACE
        FUNCTION solve_interface_Iterative(this, A, b, x0, params) RESULT(x)
            IMPORT :: dp
            IMPORT :: IterativeParams
            IMPORT :: IterativeMethod
            CLASS(IterativeMethod), INTENT(IN) :: this
            REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
            REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
            TYPE(IterativeParams), INTENT(IN) :: params
            REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        END FUNCTION solve_interface_Iterative
    END INTERFACE

    CONTAINS

   SUBROUTINE set_method(this, method)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
        TYPE(MethodTypeIterative), INTENT(IN) :: method

        SELECT CASE (method%value)
        CASE (METHOD_Jacobi%value)
            this%solve_method => solve_Jacobi
            this%method_type = METHOD_Jacobi
        CASE (METHOD_GAUSS_SEIDEL%value)
            this%solve_method => solve_Gauss_Seidel
            this%method_type = METHOD_GAUSS_SEIDEL
        CASE (METHOD_SOR%value)
            this%solve_method => solve_SOR
            this%method_type = METHOD_SOR
        CASE (METHOD_JOR%value)
            this%solve_method => solve_JOR
            this%method_type = METHOD_JOR
        CASE (METHOD_SIP_ILU%value)
            this%solve_method => solve_SIP_ILU
            this%method_type = METHOD_SIP_ILU
        CASE (METHOD_SIP_ICF%value)
            this%solve_method => solve_SIP_ICF
            this%method_type = METHOD_SIP_ICF
        CASE (METHOD_SSOR%value)
            this%solve_method => solve_SSOR
            this%method_type = METHOD_SSOR
        CASE (METHOD_RICHARDSON%value)
            this%solve_method => solve_Richardson
            this%method_type = METHOD_RICHARDSON
        CASE DEFAULT
            STOP "ERROR :: Unknown method iterative"
        END SELECT
        
    END SUBROUTINE set_method

    FUNCTION Init_IterativeParams(this, N, A, x0, max_iter_choice, epsi_tol, omega, method_preconditioner, alpha) RESULT(params)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: N
        REAL(dp), DIMENSION(:, :), OPTIONAL, INTENT(IN) :: A
        REAL(dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: x0
        INTEGER, OPTIONAL, INTENT(IN) :: max_iter_choice
        REAL(dp), OPTIONAL, INTENT(IN) :: epsi_tol
        REAL(dp), OPTIONAL, INTENT(IN) :: omega
        REAL(dp), OPTIONAL, INTENT(IN) :: alpha
        TYPE(MethodPreconditioner), OPTIONAL, INTENT(IN) :: method_preconditioner
        TYPE(IterativeParams) :: params
        INTEGER :: allocate_status

        ALLOCATE(params%x_init(N), STAT=allocate_status)
        IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate x_init"
        IF (PRESENT(x0)) THEN
            params%x_init = x0
        ELSE
            params%x_init = 0.0_dp
        END IF
        IF(PRESENT(max_iter_choice)) params%max_iter = 1000
        IF(PRESENT(epsi_tol)) params%tol = 1.0e-6_dp
        IF(PRESENT(omega)) params%omega = omega
        IF(PRESENT(alpha)) params%alpha = alpha

        SELECT CASE (this%method_type%value)
        CASE (METHOD_SIP_ILU%value)
            ALLOCATE(params%L(N,N), STAT=allocate_status)
            IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
            ALLOCATE(params%U(N,N), STAT=allocate_status)
            IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate U"
            CALL ILU_decomposition(A, params%L, params%U)
        CASE (METHOD_SIP_ICF%value)
            ALLOCATE(params%L(N,N), STAT=allocate_status)
            IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
            CALL Incomplete_Cholesky_decomposition(A, params%L)
        END SELECT

        IF(PRESENT(method_preconditioner))THEN
            params%precond => ApplyPreconditioner
            SELECT CASE (method_preconditioner%value)
            CASE (METHOD_PRECOND_JACOBI%value)
                ALLOCATE(params%D(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate D"
                params%D = Calculate_Jacobi_preconditioner(A)
                this%preconditioner_type = METHOD_PRECOND_JACOBI
            CASE (METHOD_PRECOND_GS%value)
                ALLOCATE(params%L(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                params%L = Calculate_Gauss_Seidel_preconditioner(A)
                this%preconditioner_type = METHOD_PRECOND_GS
            CASE (METHOD_PRECOND_SOR%value)
                ALLOCATE(params%L(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                params%L = Calculate_SOR_preconditioner(A, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_SOR
            CASE (METHOD_PRECOND_JOR%value)
                ALLOCATE(params%D(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate D"
                params%D = Calculate_JOR_preconditioner(A, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_JOR
            CASE (METHOD_PRECOND_ILU%value)
                ALLOCATE(params%L(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                ALLOCATE(params%U(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate U"
                CALL Calculate_ILU_preconditioner(A, params%L, params%U, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_ILU
            CASE (METHOD_PRECOND_ICF%value)
                ALLOCATE(params%L(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                params%L = Calculate_ICF_preconditioner(A, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_ICF
            CASE (METHOD_PRECOND_SSOR%value)
                ALLOCATE(params%L(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                ALLOCATE(params%D(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate D"
                ALLOCATE(params%U(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate U"
                CALL Calculate_SSOR_preconditioner(A, params%L, params%D, params%U, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_SSOR
            CASE DEFAULT
                STOP "ERROR :: Unknown method "
            END SELECT
        END IF

    END FUNCTION Init_IterativeParams

    SUBROUTINE Dealocate_IterativeParams(this, params, success)
        CLASS(IterativeMethod), INTENT(IN) :: this
        TYPE(IterativeParams), INTENT(INOUT) :: params
        LOGICAL, OPTIONAL, INTENT(OUT) :: success
        INTEGER :: deallocate_status

        IF(PRESENT(success)) success = .TRUE.

        IF (ALLOCATED(params%x_init)) DEALLOCATE(params%x_init, STAT=deallocate_status)
        IF (ALLOCATED(params%L)) DEALLOCATE(params%L, STAT=deallocate_status)
        IF (ALLOCATED(params%U)) DEALLOCATE(params%U, STAT=deallocate_status)
        IF (ALLOCATED(params%D)) DEALLOCATE(params%D, STAT=deallocate_status)
        IF (ALLOCATED(params%residual)) DEALLOCATE(params%residual, STAT=deallocate_status)

        IF(deallocate_status /= 0 .AND. PRESENT(success)) success = .FALSE.

    END SUBROUTINE Dealocate_IterativeParams

    FUNCTION IterativeMethod_solve(this, A, b, params, verbose) RESULT(x)

        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        TYPE(IterativeParams), INTENT(INOUT) :: params
        TYPE(Logger), OPTIONAL :: verbose
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x0, x_new, residu
        INTEGER :: k, N
        CHARACTER(LEN=64) :: msg

        N = SIZE(A, 1)

        x0 = params%x_init

        params%residual = b - MATMUL(A, x0)

        DO k = 1, params%max_iter

            IF(k == params%max_iter) THEN
                PRINT*, "WARNING :: non-convergence of the iterative method "//this%method_type%name
                PRINT*, "Residual norm: ", NORM2(residu)
                EXIT
            END IF

            x_new = this%solve_method(A, b, x0, params)

            residu = b - MATMUL(A, x_new)
            params%residual = residu

            IF (PRESENT(verbose)) THEN
                verbose%step = k
                WRITE(msg, '(A,I5,A,ES14.7)') "Iter ", k, " | Norm residu: ", NORM2(residu)
                CALL verbose%log(msg)
            END IF

            IF (NORM2(residu) < params%tol) EXIT

            x0 = x_new

        END DO

        x = x_new

    END FUNCTION IterativeMethod_solve

    !> Jacobi iterative method
    !>
    !> This subroutine implements the Jacobi method for solving linear systems.
    FUNCTION solve_Jacobi(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(IN) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x0(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x(i) = x(i) / A(i, i)
        END DO

    END FUNCTION solve_Jacobi

    !> Gauss-Seidel iterative method
    !>
    !> This subroutine implements the Gauss-Seidel method for solving linear systems.
    FUNCTION solve_Gauss_Seidel(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(IN) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x(i) = x(i) / A(i, i)
        END DO

    END FUNCTION solve_Gauss_Seidel

    !> Successive Over-Relaxation (SOR) iterative method
    !>
    !> This subroutine implements the SOR method for solving linear systems.
    FUNCTION solve_SOR(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(IN) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x(i) = params%omega*(x(i) / A(i, i) - x0(i)) + x0(i)
        END DO

    END FUNCTION solve_SOR

    !> Jacobi over-relaxation (JOR) iterative method
    !>
    !> This subroutine implements the Jacobi over-relaxation method for solving linear systems.
    FUNCTION solve_JOR(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(IN) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x0(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x(i) = params%omega*(x(i) / A(i, i) - x0(i)) + x0(i)
        END DO

    END FUNCTION solve_JOR

    !> strongly implicit procedure (SIP) method (or stone's method)
    !>
    !> This subroutine implements the SIP method for solving linear systems.
    !> It uses the incomplete LU decomposition of the matrix A.
    FUNCTION solve_SIP_ILU(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(IN) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1)) :: y, z

        IF(.NOT. ALLOCATED(params%L) .OR. .NOT. ALLOCATED(params%U)) STOP "ERROR :: Incomplete LU decomposition not initialized"

        y = forward(params%L, params%residual)
        
        z = backward(params%U, y)
        
        x = x0 + params%omega * z

    END FUNCTION solve_SIP_ILU

    !> strongly implicit procedure (SIP) method (or stone's method)
    !>
    !> This subroutine implements the SIP method for solving linear systems.
    !> It uses the incomplete Cholesky decomposition of the matrix A.
    FUNCTION solve_SIP_ICF(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(IN) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1)) :: y, z

        IF(.NOT. ALLOCATED(params%L)) STOP "ERROR :: Incomplete LU decomposition not initialized"

        y = forward(params%L, params%residual)

        z = backward(TRANSPOSE(params%L), y)

        x = x0 + params%omega * z

    END FUNCTION solve_SIP_ICF

    !> Symmetric successive Over-Relaxation (SSOR) iterative method
    !>
    !> This subroutine implements the SSOR method for solving linear systems.
    FUNCTION solve_SSOR(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(IN) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x_tmp
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x_tmp(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x_tmp(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x_tmp(i) = params%omega*(x_tmp(i) / A(i, i) - x0(i)) + x0(i)
        END DO

        ! backward
        DO i = N, 1, -1
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x_tmp(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x(i+1:N))
            x(i) = params%omega*(x(i) / A(i, i) - x_tmp(i)) + x_tmp(i)
        END DO

    END FUNCTION solve_SSOR

    FUNCTION solve_Richardson(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(IN) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x

        IF(this%preconditioner_type%value == METHOD_PRECOND_NONE%value)THEN
            x = x0 + params%alpha * params%residual
        ELSE IF(this%preconditioner_type%value /= METHOD_PRECOND_NONE%value)THEN
            x = x0 + params%precond(this%preconditioner_type)
        END IF

    END FUNCTION solve_Richardson

END MODULE NAFPack_Iterative_methods