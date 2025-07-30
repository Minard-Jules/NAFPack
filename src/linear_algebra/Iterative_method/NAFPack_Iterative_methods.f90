!> Module for iterative methods in NAFPack
MODULE NAFPack_Iterative_methods

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    USE NAFPack_matricielle
    USE NAFPack_Iterative_types
    USE NAFPack_Logger_mod
    USE NAFPack_Preconditioners
    USE NAFPack_Iterative_Params
    USE NAFPack_matrix_properties

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
    PUBLIC :: METHOD_CONJUGATE_GRADIENT
    PUBLIC :: METHOD_CONJUGATE_RESIDUAL

    PUBLIC :: MethodPreconditioner
    PUBLIC :: METHOD_PRECOND_NONE
    PUBLIC :: METHOD_PRECOND_JACOBI, METHOD_PRECOND_JOR
    PUBLIC :: METHOD_PRECOND_GS, METHOD_PRECOND_SOR, METHOD_PRECOND_SSOR
    PUBLIC :: METHOD_PRECOND_ILU, METHOD_PRECOND_ICF

    TYPE :: IterativeMethod
        PRIVATE
        TYPE(MethodTypeIterative) :: method_type = METHOD_ITERATIVE_NONE
        TYPE(MethodPreconditioner) :: preconditioner_type = METHOD_PRECOND_NONE
        TYPE(IterativeMethodRequirements) :: requirements
        PROCEDURE(solve_interface_Iterative), PASS(this), POINTER :: solve_method => NULL()

        CONTAINS

        PROCEDURE :: set_method => set_method
        PROCEDURE :: solve => IterativeMethod_solve
        PROCEDURE :: Init_IterativeParams => Init_IterativeParams
        PROCEDURE :: Dealocate_IterativeParams => Dealocate_IterativeParams
        PROCEDURE :: test_matrix => test_matrix

    END TYPE IterativeMethod

    ABSTRACT INTERFACE
        FUNCTION solve_interface_Iterative(this, A, b, x0, params) RESULT(x)
            IMPORT :: dp
            IMPORT :: IterativeParams
            IMPORT :: IterativeMethod
            CLASS(IterativeMethod), INTENT(IN) :: this
            REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
            REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
            TYPE(IterativeParams), INTENT(INOUT) :: params
            REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        END FUNCTION solve_interface_Iterative
    END INTERFACE

    CONTAINS

   SUBROUTINE set_method(this, method)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
        TYPE(MethodTypeIterative), INTENT(IN) :: method

        this%requirements = IterativeMethodRequirements()

        SELECT CASE (method%value)
        CASE (METHOD_Jacobi%value)
            this%solve_method => solve_Jacobi
            this%method_type = METHOD_Jacobi
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_diag_dom = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_GAUSS_SEIDEL%value)
            this%solve_method => solve_Gauss_Seidel
            this%method_type = METHOD_GAUSS_SEIDEL
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_diag_dom = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_SOR%value)
            this%solve_method => solve_SOR
            this%method_type = METHOD_SOR
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_diag_dom = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_JOR%value)
            this%solve_method => solve_JOR
            this%method_type = METHOD_JOR
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_diag_dom = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_SIP_ILU%value)
            this%solve_method => solve_SIP_ILU
            this%method_type = METHOD_SIP_ILU
            this%requirements%needs_square = .TRUE.
        CASE (METHOD_SIP_ICF%value)
            this%solve_method => solve_SIP_ICF
            this%method_type = METHOD_SIP_ICF
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_SSOR%value)
            this%solve_method => solve_SSOR
            this%method_type = METHOD_SSOR
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_RICHARDSON%value)
            this%solve_method => solve_Richardson
            this%method_type = METHOD_RICHARDSON
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_CONJUGATE_GRADIENT%value)
            this%solve_method => solve_ConjugateGradient
            this%method_type = METHOD_CONJUGATE_GRADIENT
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_CONJUGATE_RESIDUAL%value)
            this%solve_method => solve_ConjugateResidual
            this%method_type = METHOD_CONJUGATE_RESIDUAL
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_symetric = .TRUE.
        CASE DEFAULT
            STOP "ERROR :: Unknown method iterative"
        END SELECT
        
    END SUBROUTINE set_method

    FUNCTION Init_IterativeParams(this, N, A, x0, max_iter_choice, epsi_tol, omega, &
                                  method_preconditioner, alpha, is_stationary, is_strict_mode) RESULT(params)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: N
        REAL(dp), DIMENSION(:, :), OPTIONAL, INTENT(IN) :: A
        REAL(dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: x0
        INTEGER, OPTIONAL, INTENT(IN) :: max_iter_choice
        REAL(dp), OPTIONAL, INTENT(IN) :: epsi_tol
        REAL(dp), OPTIONAL, INTENT(IN) :: omega
        REAL(dp), OPTIONAL, INTENT(IN) :: alpha
        TYPE(MethodPreconditioner), OPTIONAL, INTENT(IN) :: method_preconditioner
        LOGICAL, OPTIONAL, INTENT(IN) :: is_stationary
        LOGICAL, OPTIONAL, INTENT(IN) :: is_strict_mode
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

        ALLOCATE(params%residual(N), STAT=allocate_status)
        IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate residual"

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
        CASE (METHOD_CONJUGATE_GRADIENT%value)
            ALLOCATE(params%p(N), STAT=allocate_status)
            IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate optimal descent direction p"
        CASE (METHOD_CONJUGATE_RESIDUAL%value)
            ALLOCATE(params%p(N), STAT=allocate_status)
            IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate optimal descent direction p"
        END SELECT

        IF(PRESENT(method_preconditioner))THEN
            params%precond => ApplyPreconditioner
            SELECT CASE (method_preconditioner%value)
            CASE (METHOD_PRECOND_JACOBI%value)
                ALLOCATE(params%D(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate D"
                params%D = Calculate_Jacobi_preconditioner(A)
                this%preconditioner_type = METHOD_PRECOND_JACOBI
                this%requirements%needs_diag_dom = .TRUE.
                this%requirements%needs_SPD = .TRUE.
            CASE (METHOD_PRECOND_GS%value)
                ALLOCATE(params%L(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                params%L = Calculate_Gauss_Seidel_preconditioner(A)
                this%preconditioner_type = METHOD_PRECOND_GS
                this%requirements%needs_diag_dom = .TRUE.
                this%requirements%needs_SPD = .TRUE.
            CASE (METHOD_PRECOND_SOR%value)
                ALLOCATE(params%L(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                params%L = Calculate_SOR_preconditioner(A, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_SOR
                this%requirements%needs_diag_dom = .TRUE.
                this%requirements%needs_SPD = .TRUE.
            CASE (METHOD_PRECOND_JOR%value)
                ALLOCATE(params%D(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate D"
                params%D = Calculate_JOR_preconditioner(A, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_JOR
                this%requirements%needs_diag_dom = .TRUE.
                this%requirements%needs_SPD = .TRUE.
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
                this%requirements%needs_SPD = .TRUE.
            CASE (METHOD_PRECOND_SSOR%value)
                ALLOCATE(params%L(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                ALLOCATE(params%D(N,N), STAT=allocate_status)
                IF(allocate_status /= 0) STOP "ERROR :: Unable to allocate D"
                CALL Calculate_SSOR_preconditioner(A, params%L, params%D, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_SSOR
                this%requirements%needs_SPD = .TRUE.
            CASE DEFAULT
                STOP "ERROR :: Unknown method "
            END SELECT
        END IF

        IF (PRESENT(is_stationary)) THEN
            IF(is_stationary) THEN
                params%is_stationary = .TRUE.
                this%requirements%needs_SPD = .FALSE.
            ELSE
                params%is_stationary = .FALSE.
            END IF
        END IF

        IF (PRESENT(is_strict_mode)) THEN
            IF(is_strict_mode) THEN
                params%strict_mode = .TRUE.
            ELSE
                params%strict_mode = .FALSE.
            END IF
        END IF

    END FUNCTION Init_IterativeParams

    SUBROUTINE test_matrix(this, A, params)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: A
        TYPE(IterativeParams), INTENT(IN) :: params

        
        IF (this%requirements%needs_square) THEN
            PRINT*, "Checking if the matrix is square..."
            IF (.NOT. is_square_matrix(A)) THEN
                IF(params%strict_mode) THEN
                    PRINT*, "ERROR :: "//this%method_type%name//" method requires a square matrix."
                    STOP
                ELSE
                    PRINT*, "WARNING :: "//this%method_type%name//" method requires a square matrix."
                END IF
            END IF
        END IF

        IF (this%requirements%needs_SPD) THEN
            PRINT*, "Checking if the matrix is symmetric positive definite (SPD)..."
            IF (.NOT. is_SPD(A)) THEN
                IF(params%strict_mode) THEN
                    PRINT*, "ERROR :: "//this%method_type%name//" method requires a symmetric positive definite matrix."
                    STOP
                ELSE
                    PRINT*, "WARNING :: "//this%method_type%name//" method requires a symmetric positive definite matrix."
                END IF
            END IF
        END IF

        IF (this%requirements%needs_diag_dom) THEN
            PRINT*, "Checking if the matrix is diagonally dominant..."
            IF (.NOT. is_diagonally_dominant(A)) THEN
                IF(params%strict_mode) THEN
                    PRINT*, "ERROR :: "//this%method_type%name//" method requires a diagonally dominant matrix."
                    STOP
                ELSE
                    PRINT*, "WARNING :: "//this%method_type%name//" method requires a diagonally dominant matrix."
                END IF
            END IF
        END IF

        IF (this%requirements%needs_symetric) THEN
            PRINT*, "Checking if the matrix is symmetric..."
            IF (.NOT. is_symmetric(A)) THEN
                IF(params%strict_mode) THEN
                    PRINT*, "ERROR :: "//this%method_type%name//" method requires a symmetric matrix."
                    STOP
                ELSE
                    PRINT*, "WARNING :: "//this%method_type%name//" method requires a symmetric matrix."
                END IF
            END IF
        END IF

    END SUBROUTINE test_matrix

    SUBROUTINE Dealocate_IterativeParams(this, params, success)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
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
        IF(this%preconditioner_type%value /= METHOD_PRECOND_NONE%value) this%preconditioner_type%value = METHOD_PRECOND_NONE%value
        this%requirements = IterativeMethodRequirements()

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

        IF(this%preconditioner_type%value == METHOD_PRECOND_NONE%value .AND. &
          (this%method_type%value == METHOD_CONJUGATE_GRADIENT%value .OR. &
           this%method_type%value == METHOD_CONJUGATE_RESIDUAL%value))THEN
            params%p = params%residual
        ELSE IF(this%preconditioner_type%value /= METHOD_PRECOND_NONE%value .AND. &
               (this%method_type%value == METHOD_CONJUGATE_GRADIENT%value .OR. &
                this%method_type%value == METHOD_CONJUGATE_RESIDUAL%value))THEN
            params%p = params%precond(this%preconditioner_type)
        END IF

        DO k = 1, params%max_iter
            params%k = k
            IF(k == params%max_iter) THEN
                PRINT*, "WARNING :: non-convergence of the iterative method "//this%method_type%name
                PRINT*, "Residual norm: ", NORM2(residu)
                EXIT
            END IF

            x_new = this%solve_method(A, b, x0, params)

            residu = b - MATMUL(A, x_new)
            params%residual = residu

            IF (PRESENT(verbose)) THEN
                IF (verbose%show_iteration) THEN
                    verbose%step = k
                    WRITE(msg, '(A,I5,A,ES14.7)') "Iter ", k, " | Norm residu: ", NORM2(residu)
                    CALL verbose%log(msg)
                END IF
            END IF

            IF (NORM2(residu) < params%tol) EXIT

            x0 = x_new

        END DO
        
        IF (PRESENT(verbose)) THEN
            IF (verbose%show_final) THEN
                verbose%step = k
                WRITE(msg, '(A,I5,A,ES14.7)') "Iter ", k, " | Norm residu: ", NORM2(residu)
                CALL verbose%log(msg)
            END IF
        END IF

        x = x_new

    END FUNCTION IterativeMethod_solve

    !> Jacobi iterative method
    !>
    !> This subroutine implements the Jacobi method for solving linear systems.
    FUNCTION solve_Jacobi(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
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
        TYPE(IterativeParams), INTENT(INOUT) :: params
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
        TYPE(IterativeParams), INTENT(INOUT) :: params
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
        TYPE(IterativeParams), INTENT(INOUT) :: params
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
        TYPE(IterativeParams), INTENT(INOUT) :: params
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
        TYPE(IterativeParams), INTENT(INOUT) :: params
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
        TYPE(IterativeParams), INTENT(INOUT) :: params
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
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1)) :: z
    

        IF(this%preconditioner_type%value == METHOD_PRECOND_NONE%value)THEN
            IF(.NOT. params%is_stationary) THEN
                params%alpha = DOT_PRODUCT(params%residual, params%residual) / &
                               DOT_PRODUCT(params%residual, MATMUL(A, params%residual))
            END IF
            x = x0 + params%alpha * params%residual
        ELSE IF(this%preconditioner_type%value /= METHOD_PRECOND_NONE%value)THEN
            z = params%precond(this%preconditioner_type)
            IF(.NOT. params%is_stationary) THEN
                params%alpha = DOT_PRODUCT(params%residual, z) / &
                               DOT_PRODUCT(z, MATMUL(A, z))
            END IF
            x = x0 + params%alpha * z
        END IF

    END FUNCTION solve_Richardson

    FUNCTION solve_ConjugateGradient(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1)) :: z

        IF(this%preconditioner_type%value == METHOD_PRECOND_NONE%value)THEN
            IF(params%k /= 1)THEN
                params%beta = DOT_PRODUCT(params%residual, params%residual) / params%old_dot_product
                params%p = params%residual + params%beta * params%p
            END IF

            params%alpha = DOT_PRODUCT(params%residual, params%residual) / DOT_PRODUCT(params%p, MATMUL(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = DOT_PRODUCT(params%residual, params%residual)
        ELSE IF(this%preconditioner_type%value /= METHOD_PRECOND_NONE%value)THEN
            IF(this%preconditioner_type%value == METHOD_PRECOND_GS%value .AND. &
               this%preconditioner_type%value == METHOD_PRECOND_SOR%value) THEN
                STOP "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            END IF

            IF(params%k == 1) THEN
                z = params%p
            ELSE IF(params%k /= 1)THEN
                z = params%precond(this%preconditioner_type)
                params%beta = DOT_PRODUCT(params%residual, z) / params%old_dot_product
                params%p = z + params%beta * params%p
            END IF

            params%alpha = DOT_PRODUCT(params%residual, z) / DOT_PRODUCT(params%p, MATMUL(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = DOT_PRODUCT(params%residual, z)
        END IF
        
    END FUNCTION solve_ConjugateGradient

    FUNCTION solve_ConjugateResidual(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x
        REAL(dp), DIMENSION(SIZE(A, 1)) :: z

        IF(this%preconditioner_type%value == METHOD_PRECOND_NONE%value)THEN
            IF(params%k /= 1)THEN
                params%beta = DOT_PRODUCT(params%residual, MATMUL(A,params%residual)) / params%old_dot_product
                params%p = params%residual + params%beta * params%p
            END IF

            params%alpha = DOT_PRODUCT(params%residual, MATMUL(A,params%residual)) / &
                            DOT_PRODUCT(MATMUL(A, params%p), MATMUL(A, params%p))
            
            x = x0 + params%alpha * params%p

            params%old_dot_product = DOT_PRODUCT(params%residual, MATMUL(A,params%residual))
        ELSE IF(this%preconditioner_type%value /= METHOD_PRECOND_NONE%value)THEN
            IF(this%preconditioner_type%value == METHOD_PRECOND_GS%value .AND. &
               this%preconditioner_type%value == METHOD_PRECOND_SOR%value) THEN
                STOP "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            END IF

            IF(params%k == 1) THEN
                z = params%p
            ELSE IF(params%k /= 1)THEN
                z = params%precond(this%preconditioner_type)
                params%beta = DOT_PRODUCT(z, MATMUL(A,params%residual)) / params%old_dot_product
                params%p = z + params%beta * params%p
            END IF

            params%alpha = DOT_PRODUCT(z, MATMUL(A,params%residual)) / &
                            DOT_PRODUCT(MATMUL(A, params%p), MATMUL(A, params%p))
            
            x = x0 + params%alpha * params%p

            params%old_dot_product = DOT_PRODUCT(z, MATMUL(A,params%residual))
        END IF
        
    END FUNCTION solve_ConjugateResidual 

END MODULE NAFPack_Iterative_methods