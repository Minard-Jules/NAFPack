MODULE NAFPack_Iterative_Params

    USE NAFPack_constant
    USE NAFPack_Iterative_types
    USE NAFPack_Preconditioners
    USE NAFPack_matrix_decomposition

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE

    PUBLIC :: IterativeParams
    PUBLIC :: ApplyPreconditioner

    TYPE :: IterativeParams
        ! Solution and initial guess
        REAL(dp), DIMENSION(:), ALLOCATABLE :: x_init
        ! Preconditioner matrices
        REAL(dp), DIMENSION(:, :), ALLOCATABLE :: L, U, D
        ! Algorithm vectors
        REAL(dp), DIMENSION(:), ALLOCATABLE :: p
        REAL(dp), DIMENSION(:), ALLOCATABLE :: residual
        ! Norms and tolerances
        REAL(dp) :: norm_residual
        REAL(dp) :: norm_initial_residual = 1.d0
        REAL(dp) :: tol = 1.0d-12
        ! Iteration control
        INTEGER :: k = 0
        INTEGER :: max_iter = 1000
        ! Method parameters
        REAL(dp) :: omega = 1.d0
        REAL(dp) :: alpha = 1.d0
        REAL(dp) :: beta = 1.d0
        ! ILU/IC fill level
        TYPE(FILL_LEVEL_USED) :: fill_level = FILL_LEVEL_NONE
        ! Flags
        LOGICAL :: is_stationary = .TRUE.
        LOGICAL :: strict_mode = .FALSE.
        ! Miscellaneous
        REAL(dp) :: old_dot_product = 0.d0
        TYPE(Norm_used) :: norm = NORM_2
        ! Preconditioner procedure pointer
        PROCEDURE(ApplyPreconditioner), PASS(params), POINTER :: precond

    CONTAINS

        PROCEDURE :: norm_function
    END TYPE IterativeParams

CONTAINS

    FUNCTION ApplyPreconditioner(params, method, x) RESULT(y)
        CLASS(IterativeParams), INTENT(IN) :: params
        CLASS(MethodPreconditioner), INTENT(IN) :: method
        REAL(dp), DIMENSION(:), INTENT(IN) :: x
        REAL(dp), DIMENSION(size(params%x_init)) :: y

        SELECT CASE (method%id)
        CASE (METHOD_PRECOND_JACOBI%id)
            IF (.NOT. allocated(params%D)) STOP "ERROR :: Jacobi preconditioner requires &
                                                &preconditioner matrix D to be allocated"
            y = matmul(params%D, x)
        CASE (METHOD_PRECOND_GS%id)
            IF (.NOT. allocated(params%L)) STOP "ERROR :: Gauss-Seidel preconditioner requires &
                                                &preconditioner matrix L to be allocated"
            y = forward(params%L, x)
        CASE (METHOD_PRECOND_SOR%id)
            IF (.NOT. allocated(params%L)) STOP "ERROR :: SOR preconditioner requires &
                                                &preconditioner matrix L to be allocated"
            y = forward(params%L, x)
        CASE (METHOD_PRECOND_JOR%id)
            IF (.NOT. allocated(params%D)) STOP "ERROR :: JOR preconditioner requires &
                                                &preconditioner matrix D to be allocated"
            y = matmul(params%D, x)
        CASE (METHOD_PRECOND_ILU%id)
            IF (.NOT. allocated(params%L) .OR. &
               .NOT. allocated(params%U)) STOP "ERROR :: ILU preconditioner requires &
                                                &preconditioner matrices L and U to be allocated"
            y = forward(params%L, x)
            y = backward(params%U, y)
        CASE (METHOD_PRECOND_ICF%id)
            IF (.NOT. allocated(params%L)) STOP "ERROR :: ICF preconditioner requires &
                                                &preconditioner matrix L to be allocated"
            y = forward(params%L, x)
            y = backward(transpose(params%L), y)
        CASE (METHOD_PRECOND_SSOR%id)
            IF (.NOT. allocated(params%L) .OR. &
               .NOT. allocated(params%D)) STOP "ERROR :: SSOR preconditioner requires &
                                                &preconditioner matrices L and D to be allocated"
            y = forward(params%L, x)
            y = matmul(params%D, y)
            y = backward(transpose(params%L), y)
        CASE DEFAULT
            STOP "ERROR :: Unknown preconditioner method"
        END SELECT

    END FUNCTION ApplyPreconditioner

    FUNCTION norm_function(this, vector) RESULT(RESULT)
        CLASS(IterativeParams), INTENT(IN) :: this
        REAL(dp), DIMENSION(:), INTENT(IN) :: vector
        REAL(dp) :: RESULT

        SELECT CASE (this%norm%id)
        CASE (NORM_1%id)
            RESULT = sum(abs(vector))
        CASE (NORM_2%id)
            RESULT = norm2(vector)
        CASE (NORM_INF%id)
            RESULT = maxval(abs(vector))
        CASE DEFAULT
            STOP "ERROR :: Unknown norm type"
        END SELECT

    END FUNCTION norm_function

END MODULE NAFPack_Iterative_Params
