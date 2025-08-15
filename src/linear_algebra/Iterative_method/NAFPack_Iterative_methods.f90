!==========================================================
! NAFPack_Iterative_methods.f90
! Module for iterative methods in NAFPack
!==========================================================

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
    USE NAFPack_memory_monitor

    IMPLICIT NONE(TYPE, EXTERNAL)

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
    PUBLIC :: METHOD_CGNR, METHOD_CGNE
    PUBLIC :: METHOD_GMRES

    PUBLIC :: MethodPreconditioner
    PUBLIC :: METHOD_PRECOND_NONE
    PUBLIC :: METHOD_PRECOND_JACOBI, METHOD_PRECOND_JOR
    PUBLIC :: METHOD_PRECOND_GS, METHOD_PRECOND_SOR, METHOD_PRECOND_SSOR
    PUBLIC :: METHOD_PRECOND_ILU, METHOD_PRECOND_ICF

    PUBLIC :: Norm_used
    PUBLIC :: NORM_2, NORM_1, NORM_INF

    PUBLIC :: FILL_LEVEL_USED
    PUBLIC :: FILL_LEVEL_NONE
    PUBLIC :: FILL_LEVEL_0, FILL_LEVEL_1, FILL_LEVEL_2, FILL_LEVEL_3
    PUBLIC :: FILL_LEVEL_N

    !======================
    ! Type definitions
    !======================

    TYPE :: IterativeMethod
        PRIVATE
        TYPE(MethodTypeIterative) :: method_type = METHOD_ITERATIVE_NONE
        TYPE(MethodPreconditioner) :: preconditioner_type = METHOD_PRECOND_NONE
        TYPE(IterativeMethodRequirements) :: requirements
        TYPE(relaxation_factor_used) :: relaxation_factor = RELAXATION_FACTOR_NONE
        TYPE(relaxation_factor_used) :: relaxation_factor_preconditioner = RELAXATION_FACTOR_NONE
        PROCEDURE(solve_interface_Iterative), PASS(this), POINTER :: solve_method => null()

    CONTAINS

        PROCEDURE :: set_method => set_method
        PROCEDURE :: solve => IterativeMethod_solve
        PROCEDURE :: Init_IterativeParams => Init_IterativeParams
        PROCEDURE :: Dealocate_IterativeParams => Dealocate_IterativeParams
        PROCEDURE :: test_matrix => test_matrix

    END TYPE IterativeMethod

    !======================
    ! Interface
    !======================

    ABSTRACT INTERFACE
        FUNCTION solve_interface_Iterative(this, A, b, x0, params) RESULT(x)
            IMPORT :: dp
            IMPORT :: IterativeParams
            IMPORT :: IterativeMethod
            CLASS(IterativeMethod), INTENT(IN) :: this
            REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
            REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
            TYPE(IterativeParams), INTENT(INOUT) :: params
            REAL(dp), DIMENSION(size(A, 1)) :: x
        END FUNCTION solve_interface_Iterative
    END INTERFACE

CONTAINS

    !======================
    ! Management of iterative methods
    !======================

    SUBROUTINE set_method(this, method)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
        TYPE(MethodTypeIterative), INTENT(IN) :: method

        this%requirements = IterativeMethodRequirements()

        SELECT CASE (method%id)
        CASE (METHOD_Jacobi%id)
            this%solve_method => solve_Jacobi
            this%method_type = METHOD_Jacobi
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_diag_dom = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_GAUSS_SEIDEL%id)
            this%solve_method => solve_Gauss_Seidel
            this%method_type = METHOD_GAUSS_SEIDEL
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_diag_dom = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_SOR%id)
            this%solve_method => solve_SOR
            this%method_type = METHOD_SOR
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_diag_dom = .TRUE.
            this%requirements%needs_SPD = .TRUE.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        CASE (METHOD_JOR%id)
            this%solve_method => solve_JOR
            this%method_type = METHOD_JOR
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_diag_dom = .TRUE.
            this%requirements%needs_SPD = .TRUE.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        CASE (METHOD_SIP_ILU%id)
            this%solve_method => solve_SIP_ILU
            this%method_type = METHOD_SIP_ILU
            this%requirements%needs_square = .TRUE.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        CASE (METHOD_SIP_ICF%id)
            this%solve_method => solve_SIP_ICF
            this%method_type = METHOD_SIP_ICF
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_SPD = .TRUE.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        CASE (METHOD_SSOR%id)
            this%solve_method => solve_SSOR
            this%method_type = METHOD_SSOR
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_SPD = .TRUE.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        CASE (METHOD_RICHARDSON%id)
            this%solve_method => solve_Richardson
            this%method_type = METHOD_RICHARDSON
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_SPD = .TRUE.
            this%relaxation_factor = RELAXATION_FACTOR_ALPHA
        CASE (METHOD_CONJUGATE_GRADIENT%id)
            this%solve_method => solve_ConjugateGradient
            this%method_type = METHOD_CONJUGATE_GRADIENT
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_SPD = .TRUE.
        CASE (METHOD_CONJUGATE_RESIDUAL%id)
            this%solve_method => solve_ConjugateResidual
            this%method_type = METHOD_CONJUGATE_RESIDUAL
            this%requirements%needs_square = .TRUE.
            this%requirements%needs_symetric = .TRUE.
        CASE (METHOD_CGNR%id)
            this%solve_method => solve_CGNR
            this%method_type = METHOD_CGNR
        CASE (METHOD_CGNE%id)
            this%solve_method => solve_CGNE
            this%method_type = METHOD_CGNE
        CASE (METHOD_GMRES%id)
            this%solve_method => solve_GMRES
            this%method_type = METHOD_GMRES
            this%requirements%needs_square = .TRUE.
        CASE DEFAULT
            STOP "ERROR :: Unknown method iterative"
        END SELECT

    END SUBROUTINE set_method

    FUNCTION Init_IterativeParams(this, N, A, x0, max_iter_choice, epsi_tol, omega, Norm_choice, fill_level, &
                                  method_preconditioner, alpha, is_stationary, is_strict_mode) RESULT(params)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: N
        REAL(dp), DIMENSION(:, :), OPTIONAL, INTENT(IN) :: A
        REAL(dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: x0
        INTEGER, OPTIONAL, INTENT(IN) :: max_iter_choice
        REAL(dp), OPTIONAL, INTENT(IN) :: epsi_tol
        REAL(dp), OPTIONAL, INTENT(IN) :: omega
        REAL(dp), OPTIONAL, INTENT(IN) :: alpha
        TYPE(Norm_used), OPTIONAL, INTENT(IN) :: Norm_choice
        TYPE(MethodPreconditioner), OPTIONAL, INTENT(IN) :: method_preconditioner
        LOGICAL, OPTIONAL, INTENT(IN) :: is_stationary
        LOGICAL, OPTIONAL, INTENT(IN) :: is_strict_mode
        TYPE(FILL_LEVEL_USED), OPTIONAL, INTENT(IN) :: fill_level
        TYPE(IterativeParams) :: params
        INTEGER :: allocate_status

        ALLOCATE (params%x_init(N), STAT=allocate_status)
        IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate x_init"
        IF (present(x0)) THEN
            params%x_init = x0
        ELSE
            params%x_init = 0.0_dp
        END IF

        IF (present(max_iter_choice)) params%max_iter = 1000
        IF (present(epsi_tol)) params%tol = 1.0e-6_dp
        IF (present(omega)) params%omega = omega
        IF (present(alpha)) params%alpha = alpha
        IF (present(Norm_choice)) params%norm = Norm_choice
        IF (present(fill_level)) params%fill_level = fill_level

        ALLOCATE (params%residual(N), STAT=allocate_status)
        IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate residual"

        SELECT CASE (this%method_type%id)
        CASE (METHOD_SIP_ILU%id)
            ALLOCATE (params%L(N, N), STAT=allocate_status)
            IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
            ALLOCATE (params%U(N, N), STAT=allocate_status)
            IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate U"
            CALL ILU_decomposition(A, params%L, params%U)
        CASE (METHOD_SIP_ICF%id)
            ALLOCATE (params%L(N, N), STAT=allocate_status)
            IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
            CALL Incomplete_Cholesky_decomposition(A, params%L)
        CASE (METHOD_CONJUGATE_GRADIENT%id)
            ALLOCATE (params%p(N), STAT=allocate_status)
            IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate optimal descent direction p"
        CASE (METHOD_CONJUGATE_RESIDUAL%id)
            ALLOCATE (params%p(N), STAT=allocate_status)
            IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate optimal descent direction p"
        CASE (METHOD_CGNR%id)
            ALLOCATE (params%p(N), STAT=allocate_status)
            IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate optimal descent direction p"
        CASE (METHOD_CGNE%id)
            ALLOCATE (params%p(N), STAT=allocate_status)
            IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate optimal descent direction p"
        END SELECT

        IF (present(method_preconditioner)) THEN
            params%precond => ApplyPreconditioner
            SELECT CASE (method_preconditioner%id)
            CASE (METHOD_PRECOND_JACOBI%id)
                ALLOCATE (params%D(N, N), STAT=allocate_status)
                IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate D"
                params%D = Calculate_Jacobi_preconditioner(A)
                this%preconditioner_type = METHOD_PRECOND_JACOBI
                this%requirements%needs_diag_dom = .TRUE.
                this%requirements%needs_SPD = .TRUE.
            CASE (METHOD_PRECOND_GS%id)
                ALLOCATE (params%L(N, N), STAT=allocate_status)
                IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                params%L = Calculate_Gauss_Seidel_preconditioner(A)
                this%preconditioner_type = METHOD_PRECOND_GS
                this%requirements%needs_diag_dom = .TRUE.
                this%requirements%needs_SPD = .TRUE.
            CASE (METHOD_PRECOND_SOR%id)
                ALLOCATE (params%L(N, N), STAT=allocate_status)
                IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                params%L = Calculate_SOR_preconditioner(A, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_SOR
                this%requirements%needs_diag_dom = .TRUE.
                this%requirements%needs_SPD = .TRUE.
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            CASE (METHOD_PRECOND_JOR%id)
                ALLOCATE (params%D(N, N), STAT=allocate_status)
                IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate D"
                params%D = Calculate_JOR_preconditioner(A, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_JOR
                this%requirements%needs_diag_dom = .TRUE.
                this%requirements%needs_SPD = .TRUE.
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            CASE (METHOD_PRECOND_ILU%id)
                ALLOCATE (params%L(N, N), STAT=allocate_status)
                IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                ALLOCATE (params%U(N, N), STAT=allocate_status)
                IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate U"
                IF (params%fill_level%id /= FILL_LEVEL_NONE%id) THEN
                    CALL Calculate_ILU_preconditioner(A, params%L, params%U, params%omega, params%alpha, params%fill_level%id)
                ELSE
                    CALL Calculate_ILU_preconditioner(A, params%L, params%U, params%omega, params%alpha)
                END IF
                this%preconditioner_type = METHOD_PRECOND_ILU
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            CASE (METHOD_PRECOND_ICF%id)
                ALLOCATE (params%L(N, N), STAT=allocate_status)
                IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                IF (params%fill_level%id /= FILL_LEVEL_NONE%id) THEN
                    params%L = Calculate_ICF_preconditioner(A, params%omega, params%alpha, params%fill_level%id)
                ELSE
                    params%L = Calculate_ICF_preconditioner(A, params%omega, params%alpha)
                END IF
                this%preconditioner_type = METHOD_PRECOND_ICF
                this%requirements%needs_SPD = .TRUE.
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            CASE (METHOD_PRECOND_SSOR%id)
                ALLOCATE (params%L(N, N), STAT=allocate_status)
                IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate L"
                ALLOCATE (params%D(N, N), STAT=allocate_status)
                IF (allocate_status /= 0) STOP "ERROR :: Unable to allocate D"
                CALL Calculate_SSOR_preconditioner(A, params%L, params%D, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_SSOR
                this%requirements%needs_SPD = .TRUE.
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            CASE DEFAULT
                STOP "ERROR :: Unknown method "
            END SELECT
        END IF

        IF (present(is_stationary)) THEN
            IF (is_stationary) THEN
                params%is_stationary = .TRUE.
                this%requirements%needs_SPD = .FALSE.
            ELSE
                params%is_stationary = .FALSE.
            END IF
        END IF

        IF (present(is_strict_mode)) THEN
            IF (is_strict_mode) THEN
                params%strict_mode = .TRUE.
            ELSE
                params%strict_mode = .FALSE.
            END IF
        END IF

    END FUNCTION Init_IterativeParams

    SUBROUTINE test_matrix(this, A, params, verbose)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        TYPE(IterativeParams), INTENT(IN) :: params
        TYPE(Logger), OPTIONAL, INTENT(INOUT) :: verbose
        CHARACTER(KIND=ucs4, LEN=100) :: msg
        LOGICAL :: show_matrix_test

        show_matrix_test = .FALSE.
        IF (present(verbose)) THEN
            IF (verbose%show_matrix_test) show_matrix_test = .TRUE.
        END IF

        IF (show_matrix_test) THEN
            CALL verbose%WRITE(center_with_fill("Testing matrix properties for method:"//trim(this%method_type%name), &
                                                100, fill_char="="), box_style="top")
            CALL verbose%WRITE(ucs4_"", box_style="middle")
        END IF

        IF (this%requirements%needs_square) THEN
            IF (show_matrix_test) CALL verbose%log_info(ucs4_"Checking if the matrix is square...")
            IF (.NOT. is_square_matrix(A)) THEN
                WRITE (msg, '(2A)') trim(this%method_type%name), " requires a square matrix"
                IF (params%strict_mode) THEN
                    IF (show_matrix_test) CALL verbose%log_error(msg)
                    STOP
                ELSE
                    IF (show_matrix_test) CALL verbose%log_warning(msg)
                END IF
            END IF
        END IF

        IF (this%requirements%needs_SPD) THEN
            IF (show_matrix_test) CALL verbose%log_info(ucs4_"Checking if the matrix is symmetric positive definite (SPD)...")
            IF (.NOT. is_SPD(A)) THEN
                WRITE (msg, '(2A)') trim(this%method_type%name), " method requires a symmetric positive definite matrix."
                IF (params%strict_mode) THEN
                    IF (show_matrix_test) CALL verbose%log_error(msg)
                    STOP
                ELSE
                    IF (show_matrix_test) CALL verbose%log_warning(msg)
                END IF
            END IF
        END IF

        IF (this%requirements%needs_diag_dom) THEN
            IF (show_matrix_test) CALL verbose%log_info(ucs4_"Checking if the matrix is diagonally dominant...")
            IF (.NOT. is_diagonally_dominant(A)) THEN
                WRITE (msg, '(2A)') trim(this%method_type%name), " method requires a diagonally dominant matrix."
                IF (params%strict_mode) THEN
                    IF (show_matrix_test) CALL verbose%log_error(msg)
                    STOP
                ELSE
                    IF (show_matrix_test) CALL verbose%log_warning(msg)
                END IF
            END IF
        END IF

        IF (this%requirements%needs_symetric) THEN
            IF (show_matrix_test) CALL verbose%log_info(ucs4_"Checking if the matrix is symmetric...")
            IF (.NOT. is_symmetric(A)) THEN
                WRITE (msg, '(2A)') trim(this%method_type%name), " method requires a symmetric matrix."
                IF (params%strict_mode) THEN
                    IF (show_matrix_test) CALL verbose%log_error(msg)
                    STOP
                ELSE
                    IF (show_matrix_test) CALL verbose%log_warning(msg)
                END IF
            END IF
        END IF

        IF (show_matrix_test) THEN
            CALL verbose%WRITE(center_with_fill("", width=100, fill_char="="), box_style="bottom")
            CALL verbose%WRITE(ucs4_"")
        END IF

    END SUBROUTINE test_matrix

    SUBROUTINE Dealocate_IterativeParams(this, params, success)
        CLASS(IterativeMethod), INTENT(INOUT) :: this
        TYPE(IterativeParams), INTENT(INOUT) :: params
        LOGICAL, OPTIONAL, INTENT(OUT) :: success
        INTEGER :: deallocate_status

        IF (present(success)) success = .TRUE.

        IF (allocated(params%x_init)) DEALLOCATE (params%x_init, STAT=deallocate_status)
        IF (allocated(params%L)) DEALLOCATE (params%L, STAT=deallocate_status)
        IF (allocated(params%U)) DEALLOCATE (params%U, STAT=deallocate_status)
        IF (allocated(params%D)) DEALLOCATE (params%D, STAT=deallocate_status)
        IF (allocated(params%residual)) DEALLOCATE (params%residual, STAT=deallocate_status)
        IF (allocated(params%p)) DEALLOCATE (params%p, STAT=deallocate_status)
        params%norm = NORM_2
        params%fill_level = FILL_LEVEL_NONE

        IF (deallocate_status /= 0 .AND. present(success)) success = .FALSE.

        this%preconditioner_type = METHOD_PRECOND_NONE
        this%method_type = METHOD_ITERATIVE_NONE
        this%relaxation_factor = RELAXATION_FACTOR_NONE
        this%relaxation_factor_preconditioner = RELAXATION_FACTOR_NONE

        this%requirements = IterativeMethodRequirements()

    END SUBROUTINE Dealocate_IterativeParams

    !======================
    ! Solve the system
    !======================

    FUNCTION IterativeMethod_solve(this, A, b, params, verbose) RESULT(x)

        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b
        TYPE(IterativeParams), INTENT(INOUT) :: params
        TYPE(Logger), OPTIONAL, INTENT(INOUT) :: verbose
        REAL(dp), DIMENSION(size(A, 1)) :: x
        REAL(dp), DIMENSION(size(A, 1)) :: x0, x_new
        INTEGER :: k, N, frequency
        INTEGER :: start_system_clock, end_system_clock, rate
        REAL(dp) :: elapsed_time
        LOGICAL :: show_info_solver
        LOGICAL :: show_iteration
        LOGICAL :: show_final

        N = size(A, 1)
        x0 = params%x_init

        params%residual = b - matmul(A, x0)
        params%norm_initial_residual = params%norm_function(params%residual)

        show_info_solver = .FALSE.
        show_iteration = .FALSE.
        show_final = .FALSE.
        frequency = 1
        IF (present(verbose)) THEN
            IF (verbose%show_info_solver) show_info_solver = .TRUE.
            IF (verbose%show_iteration) show_iteration = .TRUE.
            IF (verbose%show_final) show_final = .TRUE.
            frequency = verbose%frequency
        END IF

        IF (show_info_solver) CALL log_solver_info(this, params, verbose, N)
        IF (show_iteration) CALL log_iteration_header(verbose)

        CALL system_clock(start_system_clock, count_rate=rate)

        DO k = 1, params%max_iter

            params%k = k
            IF (k == params%max_iter) THEN
                EXIT
            END IF

            x_new = this%solve_method(A, b, x0, params)

            params%residual = b - matmul(A, x_new)
            params%norm_residual = params%norm_function(params%residual)
            IF (show_iteration .AND. mod(k, frequency) == 0) THEN
                CALL system_clock(end_system_clock)
                elapsed_time = REAL(end_system_clock - start_system_clock, dp) / REAL(rate, dp)
                CALL log_iteration_step(verbose, k, params, elapsed_time)
            END IF

            IF (params%norm_residual < params%tol) EXIT

            x0 = x_new

        END DO

        IF (show_iteration) CALL verbose%WRITE(center_with_fill("", width=100, fill_char="="), box_style="bottom")

        IF (show_final) THEN
            CALL system_clock(end_system_clock)
            elapsed_time = REAL(end_system_clock - start_system_clock, dp) / REAL(rate, dp)
            CALL log_final_result(verbose, k, params, x_new, elapsed_time, N)
        END IF

        x = x_new

    END FUNCTION IterativeMethod_solve

    !======================
    ! Log solver information
    !======================

    SUBROUTINE log_solver_info(solver, params, verbose, N)
        TYPE(IterativeMethod), INTENT(IN) :: solver
        TYPE(IterativeParams), INTENT(IN) :: params
        TYPE(Logger), INTENT(INOUT) :: verbose
        INTEGER, INTENT(IN) :: N
        CHARACTER(10) :: date, time
        CHARACTER(KIND=ucs4, LEN=100) :: msg

        CALL verbose%WRITE(ucs4_"")
        CALL verbose%WRITE(center_with_fill("Starting system solver", width=100, fill_char="="), box_style="top")
        CALL verbose%WRITE(ucs4_"", box_style="middle")

        ! call date_and_time(date,time,zone,values)
        CALL date_and_time(DATE=date, TIME=time)
        WRITE (date, '(A)') date(:4)//"-"//date(5:6)//"-"//date(7:8)
        WRITE (time, '(A)') time(:2)//":"//time(3:4)//":"//time(5:6)
        CALL log_field(verbose, "Date and time", trim(date)//" "//trim(time))

        IF (solver%relaxation_factor%id == RELAXATION_FACTOR_OMEGA%id) THEN
            CALL log_field(verbose, "Method used", trim(solver%method_type%name)//" "//trim(solver%method_type%name2))
            WRITE (msg, '(A, T40, 3A, ES0.4)') "Relaxation factor used", ": ", trim(solver%relaxation_factor%name), &
                " = ", params%omega
            CALL verbose%log_info(msg)
        ELSE IF (solver%relaxation_factor%id == RELAXATION_FACTOR_ALPHA%id) THEN
            CALL log_field(verbose, "Method used", trim(solver%method_type%name))
            WRITE (msg, '(A, T40, 3A, ES0.4)') "Relaxation factor used", ": ", trim(solver%relaxation_factor%name), &
                " = ", params%alpha
            CALL verbose%log_info(msg)
        ELSE
            CALL log_field(verbose, "Method used", trim(solver%method_type%name))
        END IF

        IF (solver%preconditioner_type%id /= METHOD_PRECOND_NONE%id) THEN
            IF (solver%relaxation_factor_preconditioner%id == RELAXATION_FACTOR_OMEGA%id) THEN
                CALL log_field(verbose, "Preconditioner used", trim(solver%preconditioner_type%name))
                WRITE (msg, '(A, T40, 3A, ES0.4)') "Relaxation factor used", ": ", &
                    trim(solver%relaxation_factor_preconditioner%name), " = ", params%omega
                CALL verbose%log_info(msg)
            ELSE
                CALL log_field(verbose, "Preconditioner used", trim(solver%preconditioner_type%name))
            END IF
        END IF

        IF (solver%preconditioner_type%id == METHOD_PRECOND_ILU%id .OR. &
            solver%preconditioner_type%id == METHOD_PRECOND_ICF%id) THEN
            IF (params%fill_level%id /= FILL_LEVEL_NONE%id) THEN
                CALL log_field(verbose, "Fill level used", trim(params%fill_level%name))
            ELSE
                CALL log_field(verbose, "Fill level used", "basic "//trim(solver%preconditioner_type%name))
            END IF
        ELSE IF (solver%method_type%id == METHOD_SIP_ILU%id .OR. &
                 solver%method_type%id == METHOD_SIP_ICF%id) THEN
            IF (params%fill_level%id /= FILL_LEVEL_NONE%id) THEN
                CALL log_field(verbose, "Fill level used", trim(params%fill_level%name))
            ELSE
                CALL log_field(verbose, "Fill level used", "basic "//trim(solver%method_type%name2))
            END IF
        END IF

        WRITE (msg, '(A,T40,A,I0,A,I0)') "System size", ": ", N, " x ", N
        CALL verbose%log_info(msg)

        WRITE (msg, '(A,T40,A,I0,A)') "Memory used", ": ", get_memory_kb(), " KB"
        CALL verbose%log_info(msg)

        CALL log_field(verbose, "Norme used", trim(params%norm%name))

        WRITE (msg, '(A,T40,A, ES0.2)') "Convergence criterion (Tolerance)", ": ||r|| < ", params%tol
        CALL verbose%log_info(msg)

        CALL log_field(verbose, "Max iterations", params%max_iter)

        WRITE (msg, '(A,T36,A,ES0.7)') "Initial residual norm", ": ||r0|| = ", params%norm_initial_residual
        CALL verbose%log_detail(msg)

        CALL verbose%WRITE(center_with_fill("", width=100, fill_char="="), box_style="bottom")
        CALL verbose%WRITE(ucs4_"")

    END SUBROUTINE log_solver_info

    SUBROUTINE log_iteration_header(verbose)
        CLASS(Logger), INTENT(INOUT) :: verbose
        CHARACTER(KIND=ucs4, LEN=100) :: msg

        CALL verbose%WRITE(center_with_fill("Iterations", width=100, fill_char="="), box_style="top")
        CALL verbose%WRITE(ucs4_"", box_style="middle")

        WRITE (msg, '(A,5X,A,5X,A,5X,A)') "Iter", &
            "Residual norm (||r||)", &
            "Relative residual norm (||r||/||r0||)", &
            "Time (s)"
        CALL verbose%log_info(msg)
        CALL verbose%WRITE(repeat(ucs4_"-", 100), box_style="middle")
    END SUBROUTINE log_iteration_header

    SUBROUTINE log_iteration_step(verbose, k, params, elapsed_time)
        CLASS(Logger), INTENT(INOUT) :: verbose
        INTEGER, INTENT(IN) :: k
        REAL(dp), INTENT(IN) :: elapsed_time
        TYPE(IterativeParams), INTENT(IN) :: params
        CHARACTER(KIND=ucs4, LEN=100) :: msg
        INTEGER :: end_system_clock

        CALL system_clock(end_system_clock)
        WRITE (msg, '(T2,I0,T15,ES0.7,T48,ES0.7,T79,ES0.7)') k, &
            params%norm_residual, &
            params%norm_residual / params%norm_initial_residual, &
            elapsed_time
        CALL verbose%log_time(msg)
    END SUBROUTINE log_iteration_step

    SUBROUTINE log_final_result(verbose, k, params, x_new, elapsed_time, N)
        CLASS(Logger), INTENT(INOUT) :: verbose
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: k
        REAL(dp), INTENT(IN) :: elapsed_time
        TYPE(IterativeParams), INTENT(IN) :: params
        REAL(dp), DIMENSION(:), INTENT(IN) :: x_new
        CHARACTER(KIND=ucs4, LEN=100) :: msg
        INTEGER :: end_system_clock
        INTEGER :: i

        CALL verbose%WRITE(ucs4_"")
        CALL verbose%WRITE(center_with_fill("Results", width=100, fill_char="="), box_style="top")
        CALL verbose%WRITE(ucs4_" ", box_style="middle")

        IF (k < params%max_iter) THEN
            CALL log_field(verbose, "Status", "CONVERGED")

            WRITE (msg, '(A,T40,A,ES0.7,A,ES0.1,A)') "Final residual", ": ||r|| = ", params%norm_residual, &
                " < ", params%tol, " (convergence achieved)"
            CALL verbose%log_info(msg)

            WRITE (msg, '(A,T36,A,ES0.7)') "Relative residual norm", ": ||r||/||r0|| = ", &
                params%norm_residual / params%norm_initial_residual
            CALL verbose%log_detail(msg)

            CALL log_field(verbose, "Total iterations", k)

            WRITE (msg, '(A,T40,A,ES0.7)') "Solution", ": x = ["
            IF (N < 6) THEN
                DO i = 1, size(x_new)
                    WRITE (msg, '(2A,ES0.7)') trim(msg), " ", x_new(i)
                END DO
            ELSE
                WRITE (msg, '(2A,ES0.7,3X,ES0.7,3X,ES0.7,3X,A,ES0.7,3X,ES0.7,3X,ES0.7,A)') trim(msg), " ", &
                    x_new(1), x_new(2), x_new(3), x_new(N - 3), x_new(N - 2), x_new(N - 1)
            END IF
            WRITE (msg, '(2A,ES0.7)') trim(msg), "]"
            CALL verbose%log_info(msg)

            CALL system_clock(end_system_clock)
            WRITE (msg, '(A,T40,A,ES0.7,A,I0)') "Solver completed in ", ": ", &
                elapsed_time, " seconds"
            CALL verbose%log_info(msg)

            WRITE (msg, '(A,T40,A,I0,A)') "Memory used", ": ", get_memory_kb(), " KB"
            CALL verbose%log_info(msg)
        ELSE
            CALL log_field(verbose, "Status", "NOT CONVERGED")

            WRITE (msg, '(A,T40,A,ES0.7,A,ES0.1,A)') "Final residual", ": ||r|| = ", params%norm_residual, &
                " >= ", params%tol, " (convergence not achieved)"
            CALL verbose%log_info(msg)
        END IF
        CALL verbose%WRITE(center_with_fill("", width=100, fill_char="="), box_style="bottom")
    END SUBROUTINE

    !======================
    ! Solve methods
    !======================

    !> Jacobi iterative method
    !>
    !> This subroutine implements the Jacobi method for solving linear systems.
    FUNCTION solve_Jacobi(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(size(A, 1)) :: x
        INTEGER :: i, N

        N = size(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - dot_product(A(i, 1:i - 1), x0(1:i - 1)) - dot_product(A(i, i + 1:N), x0(i + 1:N))
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
        REAL(dp), DIMENSION(size(A, 1)) :: x
        INTEGER :: i, N

        N = size(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - dot_product(A(i, 1:i - 1), x(1:i - 1)) - dot_product(A(i, i + 1:N), x0(i + 1:N))
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
        REAL(dp), DIMENSION(size(A, 1)) :: x
        INTEGER :: i, N

        N = size(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - dot_product(A(i, 1:i - 1), x(1:i - 1)) - dot_product(A(i, i + 1:N), x0(i + 1:N))
            x(i) = params%omega * (x(i) / A(i, i) - x0(i)) + x0(i)
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
        REAL(dp), DIMENSION(size(A, 1)) :: x
        INTEGER :: i, N

        N = size(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - dot_product(A(i, 1:i - 1), x0(1:i - 1)) - dot_product(A(i, i + 1:N), x0(i + 1:N))
            x(i) = params%omega * (x(i) / A(i, i) - x0(i)) + x0(i)
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
        REAL(dp), DIMENSION(size(A, 1)) :: x
        REAL(dp), DIMENSION(size(A, 1)) :: y, z

        IF (.NOT. allocated(params%L) .OR. .NOT. allocated(params%U)) STOP "ERROR :: Incomplete LU decomposition not initialized"

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
        REAL(dp), DIMENSION(size(A, 1)) :: x
        REAL(dp), DIMENSION(size(A, 1)) :: y, z

        IF (.NOT. allocated(params%L)) STOP "ERROR :: Incomplete LU decomposition not initialized"

        y = forward(params%L, params%residual)

        z = backward(transpose(params%L), y)

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
        REAL(dp), DIMENSION(size(A, 1)) :: x
        REAL(dp), DIMENSION(size(A, 1)) :: x_tmp
        INTEGER :: i, N

        N = size(A, 1)

        ! forward
        DO i = 1, N
            x_tmp(i) = b(i) - dot_product(A(i, 1:i - 1), x_tmp(1:i - 1)) - dot_product(A(i, i + 1:N), x0(i + 1:N))
            x_tmp(i) = params%omega * (x_tmp(i) / A(i, i) - x0(i)) + x0(i)
        END DO

        ! backward
        DO i = N, 1, -1
            x(i) = b(i) - dot_product(A(i, 1:i - 1), x_tmp(1:i - 1)) - dot_product(A(i, i + 1:N), x(i + 1:N))
            x(i) = params%omega * (x(i) / A(i, i) - x_tmp(i)) + x_tmp(i)
        END DO

    END FUNCTION solve_SSOR

    !> Richardson iterative method
    !>
    !> This subroutine implements the Richardson method for solving linear systems.
    FUNCTION solve_Richardson(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(size(A, 1)) :: x
        REAL(dp), DIMENSION(size(A, 1)) :: z_prec

        IF (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) THEN
            IF (.NOT. params%is_stationary) THEN
                params%alpha = dot_product(params%residual, params%residual) / &
                               dot_product(params%residual, matmul(A, params%residual))
            END IF
            x = x0 + params%alpha * params%residual
        ELSE IF (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) THEN
            z_prec = params%precond(this%preconditioner_type, params%residual)
            IF (.NOT. params%is_stationary) THEN
                params%alpha = dot_product(params%residual, z_prec) / &
                               dot_product(z_prec, matmul(A, z_prec))
            END IF
            x = x0 + params%alpha * z_prec
        END IF

    END FUNCTION solve_Richardson

    !> Conjugate Gradient iterative method
    !>
    !> This subroutine implements the Conjugate Gradient method for solving linear systems.
    FUNCTION solve_ConjugateGradient(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(size(A, 1)) :: x
        REAL(dp), DIMENSION(size(A, 1)) :: z_prec

        IF (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) THEN
            IF (params%k == 1) THEN
                params%p = params%residual
            ELSE IF (params%k /= 1) THEN
                params%beta = dot_product(params%residual, params%residual) / params%old_dot_product
                params%p = params%residual + params%beta * params%p
            END IF

            params%alpha = dot_product(params%residual, params%residual) / dot_product(params%p, matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(params%residual, params%residual)
        ELSE IF (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) THEN
            IF (this%preconditioner_type%id == METHOD_PRECOND_GS%id .OR. &
                this%preconditioner_type%id == METHOD_PRECOND_SOR%id) THEN
                STOP "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            END IF

            IF (params%k == 1) THEN
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%p = z_prec
            ELSE IF (params%k /= 1) THEN
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%beta = dot_product(z_prec, params%residual) / params%old_dot_product
                params%p = z_prec + params%beta * params%p
            END IF

            params%alpha = dot_product(z_prec, params%residual) / dot_product(params%p, matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(z_prec, params%residual)
        END IF

    END FUNCTION solve_ConjugateGradient

    !> Conjugate Residual iterative method
    !>
    !> This subroutine implements the Conjugate Residual method for solving linear systems.
    FUNCTION solve_ConjugateResidual(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(size(A, 1)) :: x
        REAL(dp), DIMENSION(size(A, 1)) :: z_prec

        IF (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) THEN
            IF (params%k == 1) THEN
                params%p = params%residual
            ELSE IF (params%k /= 1) THEN
                params%beta = dot_product(params%residual, matmul(A, params%residual)) / params%old_dot_product
                params%p = params%residual + params%beta * params%p
            END IF

            params%alpha = dot_product(params%residual, matmul(A, params%residual)) / &
                           dot_product(matmul(A, params%p), matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(params%residual, matmul(A, params%residual))
        ELSE IF (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) THEN
            IF (this%preconditioner_type%id == METHOD_PRECOND_GS%id .OR. &
                this%preconditioner_type%id == METHOD_PRECOND_SOR%id) THEN
                STOP "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            END IF

            IF (params%k == 1) THEN
                params%p = params%precond(this%preconditioner_type, params%residual)
                z_prec = params%p
            ELSE IF (params%k /= 1) THEN
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%beta = dot_product(z_prec, matmul(A, params%residual)) / params%old_dot_product
                params%p = z_prec + params%beta * params%p
            END IF

            params%alpha = dot_product(z_prec, matmul(A, params%residual)) / &
                           dot_product(matmul(A, params%p), matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(z_prec, matmul(A, params%residual))
        END IF

    END FUNCTION solve_ConjugateResidual

    !> Conjugate Gradient on Normal Equations iterative method
    !>
    !> This subroutine implements the Conjugate Gradient on Normal Equations method (or Craig’s Method) for solving linear systems.
    FUNCTION solve_CGNR(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(size(A, 1)) :: x
        REAL(dp), DIMENSION(size(A, 1)) :: z_prec
        REAL(dp), DIMENSION(size(A, 1)) :: AT_r

        IF (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) THEN
            IF (params%k == 1) THEN
                AT_r = matmul(transpose(A), params%residual)
                params%p = AT_r
            ELSE IF (params%k /= 1) THEN
                AT_r = matmul(transpose(A), params%residual)
                params%beta = dot_product(AT_r, AT_r) / params%old_dot_product
                params%p = AT_r + params%beta * params%p
            END IF

            params%alpha = dot_product(AT_r, AT_r) / &
                           dot_product(matmul(A, params%p), matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(AT_r, AT_r)
        ELSE IF (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) THEN
            IF (this%preconditioner_type%id == METHOD_PRECOND_GS%id .OR. &
                this%preconditioner_type%id == METHOD_PRECOND_SOR%id) THEN
                STOP "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            END IF

            IF (params%k == 1) THEN
                AT_r = matmul(transpose(A), params%residual)
                z_prec = params%precond(this%preconditioner_type, AT_r)
                params%p = z_prec
            ELSE IF (params%k /= 1) THEN
                AT_r = matmul(transpose(A), params%residual)
                z_prec = params%precond(this%preconditioner_type, AT_r)
                params%beta = dot_product(z_prec, AT_r) / params%old_dot_product
                params%p = z_prec + params%beta * params%p
            END IF

            params%alpha = dot_product(z_prec, AT_r) / &
                           dot_product(matmul(A, params%p), matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(z_prec, AT_r)
        END IF

    END FUNCTION solve_CGNR

    !> Conjugate Gradient on Normal Residual iterative method
    !>
    !> This subroutine implements the Conjugate Gradient on Normal Residual method for solving linear systems.
    FUNCTION solve_CGNE(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(size(A, 1)) :: x
        REAL(dp), DIMENSION(size(A, 1)) :: z_prec

        IF (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) THEN
            IF (params%k == 1) THEN
                params%p = matmul(transpose(A), params%residual)
            ELSE IF (params%k /= 1) THEN
                params%beta = dot_product(params%residual, params%residual) / params%old_dot_product
                params%p = matmul(transpose(A), params%residual) + params%beta * params%p
            END IF

            params%alpha = dot_product(params%residual, params%residual) / &
                           dot_product(params%p, params%p)

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(params%residual, params%residual)
        ELSE IF (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) THEN
            IF (this%preconditioner_type%id == METHOD_PRECOND_GS%id .OR. &
                this%preconditioner_type%id == METHOD_PRECOND_SOR%id) THEN
                STOP "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            END IF

            IF (params%k == 1) THEN
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%p = matmul(transpose(A), z_prec)
            ELSE IF (params%k /= 1) THEN
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%beta = dot_product(z_prec, params%residual) / params%old_dot_product
                params%p = matmul(transpose(A), z_prec) + params%beta * params%p
            END IF

            params%alpha = dot_product(z_prec, params%residual) / &
                           dot_product(params%p, params%p)

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(z_prec, params%residual)
        END IF

    END FUNCTION solve_CGNE

    FUNCTION solve_GMRES(this, A, b, x0, params) RESULT(x)
        CLASS(IterativeMethod), INTENT(IN) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        TYPE(IterativeParams), INTENT(INOUT) :: params
        REAL(dp), DIMENSION(size(A, 1)) :: x

        ! in progress

    END FUNCTION solve_GMRES

END MODULE NAFPack_Iterative_methods
