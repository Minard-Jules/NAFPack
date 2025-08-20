module NAFPack_Iterative_Params

    use NAFPack_constant
    use NAFPack_Iterative_types
    use NAFPack_Preconditioners
    use NAFPack_matrix_decomposition

    implicit none(type, external)

    private

    public :: IterativeParams
    public :: ApplyPreconditioner

    type :: IterativeParams
        ! Solution and initial guess
        real(dp), dimension(:), allocatable :: x_init
        ! Preconditioner matrices
        real(dp), dimension(:, :), allocatable :: L, U, D
        ! Algorithm vectors
        real(dp), dimension(:), allocatable :: p
        real(dp), dimension(:), allocatable :: residual
        ! Norms and tolerances
        real(dp) :: norm_residual
        real(dp) :: norm_initial_residual = 1.d0
        real(dp) :: tol = 1.0d-12
        ! Iteration control
        integer :: k = 0
        integer :: max_iter = 1000
        ! Method parameters
        real(dp) :: omega = 1.d0
        real(dp) :: alpha = 1.d0
        real(dp) :: beta = 1.d0
        ! ILU/IC fill level
        type(FILL_LEVEL_USED) :: fill_level = FILL_LEVEL_NONE
        ! Flags
        logical :: is_stationary = .true.
        logical :: strict_mode = .false.
        ! Miscellaneous
        real(dp) :: old_dot_product = 0.d0
        type(Norm_used) :: norm = NORM_2
        ! Preconditioner procedure pointer
        procedure(ApplyPreconditioner), pass(params), pointer :: precond

    contains

        procedure :: norm_function
    end type IterativeParams

contains

    function ApplyPreconditioner(params, method, x) result(y)
        class(IterativeParams), intent(IN) :: params
        class(MethodPreconditioner), intent(IN) :: method
        real(dp), dimension(:), intent(IN) :: x
        real(dp), dimension(size(params%x_init)) :: y

        select case (method%id)
        case (METHOD_PRECOND_JACOBI%id)
            if (.not. allocated(params%D)) stop "ERROR :: Jacobi preconditioner requires &
                                                &preconditioner matrix D to be allocated"
            y = matmul(params%D, x)
        case (METHOD_PRECOND_GS%id)
            if (.not. allocated(params%L)) stop "ERROR :: Gauss-Seidel preconditioner requires &
                                                &preconditioner matrix L to be allocated"
            y = forward(params%L, x)
        case (METHOD_PRECOND_SOR%id)
            if (.not. allocated(params%L)) stop "ERROR :: SOR preconditioner requires &
                                                &preconditioner matrix L to be allocated"
            y = forward(params%L, x)
        case (METHOD_PRECOND_JOR%id)
            if (.not. allocated(params%D)) stop "ERROR :: JOR preconditioner requires &
                                                &preconditioner matrix D to be allocated"
            y = matmul(params%D, x)
        case (METHOD_PRECOND_ILU%id)
            if (.not. allocated(params%L) .or. &
               .not. allocated(params%U)) stop "ERROR :: ILU preconditioner requires &
                                                &preconditioner matrices L and U to be allocated"
            y = forward(params%L, x)
            y = backward(params%U, y)
        case (METHOD_PRECOND_ICF%id)
            if (.not. allocated(params%L)) stop "ERROR :: ICF preconditioner requires &
                                                &preconditioner matrix L to be allocated"
            y = forward(params%L, x)
            y = backward(transpose(params%L), y)
        case (METHOD_PRECOND_SSOR%id)
            if (.not. allocated(params%L) .or. &
               .not. allocated(params%D)) stop "ERROR :: SSOR preconditioner requires &
                                                &preconditioner matrices L and D to be allocated"
            y = forward(params%L, x)
            y = matmul(params%D, y)
            y = backward(transpose(params%L), y)
        case DEFAULT
            stop "ERROR :: Unknown preconditioner method"
        end select

    end function ApplyPreconditioner

    function norm_function(this, vector) result(result)
        class(IterativeParams), intent(IN) :: this
        real(dp), dimension(:), intent(IN) :: vector
        real(dp) :: result

        select case (this%norm%id)
        case (NORM_1%id)
            result = sum(abs(vector))
        case (NORM_2%id)
            result = norm2(vector)
        case (NORM_INF%id)
            result = maxval(abs(vector))
        case DEFAULT
            stop "ERROR :: Unknown norm type"
        end select

    end function norm_function

end module NAFPack_Iterative_Params
