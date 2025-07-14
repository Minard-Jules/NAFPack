MODULE NAFPack_Iterative_Params

    USE NAFPack_constant
    USE NAFPack_Iterative_types
    USE NAFPack_Preconditioners
    USE NAFPack_matrix_decomposition

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: IterativeParams
    PUBLIC :: ApplyPreconditioner

    TYPE :: IterativeParams
        REAL(dp), DIMENSION(:), ALLOCATABLE :: x_init
        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: L, U, D
        REAL(dp), DIMENSION(:), ALLOCATABLE :: residual
        REAL(dp), DIMENSION(:), ALLOCATABLE :: p
        REAL(dp) :: omega = 1.d0
        INTEGER :: k = 0
        INTEGER :: max_iter = 1000
        REAL(dp) :: tol = 1.0d-12
        REAL(dp) :: alpha = 1.d0
        REAL(dp) :: beta = 1.d0
        LOGICAL :: is_stationary = .TRUE.
        REAL(dp) :: old_dot_product = 0.d0
        PROCEDURE(ApplyPreconditioner), PASS(params), POINTER :: precond
    END TYPE IterativeParams

    CONTAINS

    FUNCTION ApplyPreconditioner(params, method) RESULT(y)
        CLASS(IterativeParams), INTENT(IN) :: params
        CLASS(MethodPreconditioner), INTENT(IN) :: method
        REAL(dp), DIMENSION(SIZE(params%x_init)) :: y

        SELECT CASE(method%value)
        CASE (METHOD_PRECOND_JACOBI%value)
            IF(.NOT. ALLOCATED(params%D)) STOP "ERROR :: Jacobi preconditioner requires &
                                                &preconditioner matrix D to be allocated"
            y = MATMUL(params%D, params%residual)
        CASE (METHOD_PRECOND_GS%value)
            IF(.NOT. ALLOCATED(params%L)) STOP "ERROR :: Gauss-Seidel preconditioner requires &
                                                &preconditioner matrix L to be allocated"
            y = forward(params%L, params%residual)
        CASE (METHOD_PRECOND_SOR%value)
            IF(.NOT. ALLOCATED(params%L)) STOP "ERROR :: SOR preconditioner requires &
                                                &preconditioner matrix L to be allocated"
            y = forward(params%L, params%residual)
        CASE (METHOD_PRECOND_JOR%value)
            IF(.NOT. ALLOCATED(params%D)) STOP "ERROR :: JOR preconditioner requires &
                                                &preconditioner matrix D to be allocated"
            y = MATMUL(params%D, params%residual)
        CASE (METHOD_PRECOND_ILU%value)
            IF(.NOT. ALLOCATED(params%L) .OR. &
               .NOT. ALLOCATED(params%U)) STOP "ERROR :: ILU preconditioner requires &
                                                &preconditioner matrices L and U to be allocated"
            y = forward(params%L, params%residual)
            y = backward(params%U, y)
        CASE (METHOD_PRECOND_ICF%value)
            IF(.NOT. ALLOCATED(params%L)) STOP "ERROR :: ICF preconditioner requires &
                                                &preconditioner matrix L to be allocated"
            y = forward(params%L, params%residual)
            y = backward(TRANSPOSE(params%L), y)
        CASE (METHOD_PRECOND_SSOR%value)
            IF(.NOT. ALLOCATED(params%L) .OR. &
               .NOT. ALLOCATED(params%D)) STOP "ERROR :: SSOR preconditioner requires &
                                                &preconditioner matrices L and D to be allocated"
            y = forward(params%L, params%residual)
            y = MATMUL(params%D, y)
            y = backward(TRANSPOSE(params%L), y)
        CASE DEFAULT
            STOP "ERROR :: Unknown preconditioner method"
        END SELECT

    END FUNCTION ApplyPreconditioner

END MODULE NAFPack_Iterative_Params