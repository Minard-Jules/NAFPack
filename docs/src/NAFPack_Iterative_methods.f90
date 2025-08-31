!==========================================================
! NAFPack_Iterative_methods.f90
! Module for iterative methods in NAFPack
!==========================================================

!> Module for iterative methods in NAFPack
module NAFPack_Iterative_methods

    use NAFPack_kinds, only: dp, ucs4

    use NAFPack_matrix_decomposition, only: forward, backward, &
                                            Incomplete_Cholesky_decomposition, ILU_decomposition

    use NAFPack_Iterative_types, only: MethodTypeIterative, METHOD_ITERATIVE_NONE, &
                                       METHOD_Jacobi, METHOD_JOR, &
                                       METHOD_GAUSS_SEIDEL, METHOD_SOR, METHOD_SSOR, &
                                       METHOD_SIP_ILU, METHOD_SIP_ICF, &
                                       METHOD_RICHARDSON, &
                                       METHOD_CONJUGATE_GRADIENT, METHOD_CONJUGATE_RESIDUAL, &
                                       METHOD_CGNE, METHOD_CGNR, &
                                       METHOD_GMRES, &
                                       IterativeMethodRequirements, &
                                       Norm_used, NORM_2, NORM_1, NORM_INF, &
                                       relaxation_factor_used, RELAXATION_FACTOR_NONE, &
                                       RELAXATION_FACTOR_OMEGA, RELAXATION_FACTOR_ALPHA

    use NAFPack_Logger_mod, only: Logger, center_with_fill, log_field

    use NAFPack_Preconditioners, only: FILL_LEVEL_USED, FILL_LEVEL_NONE, &
                                       FILL_LEVEL_0, FILL_LEVEL_1, FILL_LEVEL_2, FILL_LEVEL_3, &
                                       FILL_LEVEL_N, &
                                       MethodPreconditioner, METHOD_PRECOND_NONE, &
                                       METHOD_PRECOND_JACOBI, METHOD_PRECOND_GS, &
                                       METHOD_PRECOND_SOR, METHOD_PRECOND_JOR, &
                                       METHOD_PRECOND_ILU, METHOD_PRECOND_ICF, &
                                       METHOD_PRECOND_SSOR, &
                                       Calculate_Gauss_Seidel_preconditioner, &
                                       Calculate_ICF_preconditioner, &
                                       Calculate_ILU_preconditioner, &
                                       Calculate_Jacobi_preconditioner, &
                                       Calculate_JOR_preconditioner, &
                                       Calculate_SOR_preconditioner, &
                                       Calculate_SSOR_preconditioner

    use NAFPack_Iterative_Params, only: IterativeParams, ApplyPreconditioner

    use NAFPack_matrix_properties, only: is_square_matrix, is_SPD, &
                                         is_diagonally_dominant, is_symmetric

    use NAFPack_memory_monitor, only: get_memory_kb

    implicit none(type, external)

    private

    public :: IterativeMethod
    public :: IterativeParams

    public :: MethodTypeIterative
    public :: METHOD_ITERATIVE_NONE
    public :: METHOD_Jacobi, METHOD_JOR
    public :: METHOD_GAUSS_SEIDEL, METHOD_SOR, METHOD_SSOR
    public :: METHOD_SIP_ILU, METHOD_SIP_ICF
    public :: METHOD_RICHARDSON
    public :: METHOD_CONJUGATE_GRADIENT
    public :: METHOD_CONJUGATE_RESIDUAL
    public :: METHOD_CGNR, METHOD_CGNE
    public :: METHOD_GMRES

    public :: MethodPreconditioner
    public :: METHOD_PRECOND_NONE
    public :: METHOD_PRECOND_JACOBI, METHOD_PRECOND_JOR
    public :: METHOD_PRECOND_GS, METHOD_PRECOND_SOR, METHOD_PRECOND_SSOR
    public :: METHOD_PRECOND_ILU, METHOD_PRECOND_ICF

    public :: Norm_used
    public :: NORM_2, NORM_1, NORM_INF

    public :: FILL_LEVEL_USED
    public :: FILL_LEVEL_NONE
    public :: FILL_LEVEL_0, FILL_LEVEL_1, FILL_LEVEL_2, FILL_LEVEL_3
    public :: FILL_LEVEL_N

    !======================
    ! Type definitions
    !======================

    type :: IterativeMethod
        private
        type(MethodTypeIterative) :: method_type = METHOD_ITERATIVE_NONE
        type(MethodPreconditioner) :: preconditioner_type = METHOD_PRECOND_NONE
        type(IterativeMethodRequirements) :: requirements
        type(relaxation_factor_used) :: relaxation_factor = RELAXATION_FACTOR_NONE
        type(relaxation_factor_used) :: relaxation_factor_preconditioner = RELAXATION_FACTOR_NONE
        procedure(solve_interface_Iterative), pass(this), pointer :: solve_method => null()

    contains

        procedure :: set_method => set_method
        procedure :: solve => IterativeMethod_solve
        procedure :: Init_IterativeParams => Init_IterativeParams
        procedure :: Dealocate_IterativeParams => Dealocate_IterativeParams
        procedure :: test_matrix => test_matrix

    end type IterativeMethod

    !======================
    ! Interface
    !======================

    abstract interface
        function solve_interface_Iterative(this, A, b, x0, params) result(x)
            import :: dp
            import :: IterativeParams
            import :: IterativeMethod
            implicit none(type, external)
            class(IterativeMethod), intent(in) :: this
            real(dp), dimension(:, :), intent(in) :: A
            real(dp), dimension(:), intent(in) :: b, x0
            type(IterativeParams), intent(inout) :: params
            real(dp), dimension(size(A, 1)) :: x
        end function solve_interface_Iterative
    end interface

contains

    !======================
    ! Management of iterative methods
    !======================

    subroutine set_method(this, method)
        class(IterativeMethod), intent(inout) :: this
        type(MethodTypeIterative), intent(in) :: method

        this%requirements = IterativeMethodRequirements()

        select case (method%id)
        case (METHOD_Jacobi%id)
            this%solve_method => solve_Jacobi
            this%method_type = METHOD_Jacobi
            this%requirements%needs_square = .true.
            this%requirements%needs_diag_dom = .true.
            this%requirements%needs_SPD = .true.
        case (METHOD_GAUSS_SEIDEL%id)
            this%solve_method => solve_Gauss_Seidel
            this%method_type = METHOD_GAUSS_SEIDEL
            this%requirements%needs_square = .true.
            this%requirements%needs_diag_dom = .true.
            this%requirements%needs_SPD = .true.
        case (METHOD_SOR%id)
            this%solve_method => solve_SOR
            this%method_type = METHOD_SOR
            this%requirements%needs_square = .true.
            this%requirements%needs_diag_dom = .true.
            this%requirements%needs_SPD = .true.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        case (METHOD_JOR%id)
            this%solve_method => solve_JOR
            this%method_type = METHOD_JOR
            this%requirements%needs_square = .true.
            this%requirements%needs_diag_dom = .true.
            this%requirements%needs_SPD = .true.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        case (METHOD_SIP_ILU%id)
            this%solve_method => solve_SIP_ILU
            this%method_type = METHOD_SIP_ILU
            this%requirements%needs_square = .true.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        case (METHOD_SIP_ICF%id)
            this%solve_method => solve_SIP_ICF
            this%method_type = METHOD_SIP_ICF
            this%requirements%needs_square = .true.
            this%requirements%needs_SPD = .true.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        case (METHOD_SSOR%id)
            this%solve_method => solve_SSOR
            this%method_type = METHOD_SSOR
            this%requirements%needs_square = .true.
            this%requirements%needs_SPD = .true.
            this%relaxation_factor = RELAXATION_FACTOR_OMEGA
        case (METHOD_RICHARDSON%id)
            this%solve_method => solve_Richardson
            this%method_type = METHOD_RICHARDSON
            this%requirements%needs_square = .true.
            this%requirements%needs_SPD = .true.
            this%relaxation_factor = RELAXATION_FACTOR_ALPHA
        case (METHOD_CONJUGATE_GRADIENT%id)
            this%solve_method => solve_ConjugateGradient
            this%method_type = METHOD_CONJUGATE_GRADIENT
            this%requirements%needs_square = .true.
            this%requirements%needs_SPD = .true.
        case (METHOD_CONJUGATE_RESIDUAL%id)
            this%solve_method => solve_ConjugateResidual
            this%method_type = METHOD_CONJUGATE_RESIDUAL
            this%requirements%needs_square = .true.
            this%requirements%needs_symetric = .true.
        case (METHOD_CGNR%id)
            this%solve_method => solve_CGNR
            this%method_type = METHOD_CGNR
        case (METHOD_CGNE%id)
            this%solve_method => solve_CGNE
            this%method_type = METHOD_CGNE
        case (METHOD_GMRES%id)
            this%solve_method => solve_GMRES
            this%method_type = METHOD_GMRES
            this%requirements%needs_square = .true.
        case DEFAULT
            stop "ERROR :: Unknown method iterative"
        end select

    end subroutine set_method

    function Init_IterativeParams(this, N, A, x0, max_iter_choice, epsi_tol, omega, Norm_choice, &
                                  fill_level, method_preconditioner, alpha, is_stationary, &
                                  is_strict_mode) result(params)
        class(IterativeMethod), intent(inout) :: this
        integer, intent(in) :: N
        real(dp), dimension(:, :), optional, intent(in) :: A
        real(dp), dimension(:), optional, intent(in) :: x0
        integer, optional, intent(in) :: max_iter_choice
        real(dp), optional, intent(in) :: epsi_tol
        real(dp), optional, intent(in) :: omega
        real(dp), optional, intent(in) :: alpha
        type(Norm_used), optional, intent(in) :: Norm_choice
        type(MethodPreconditioner), optional, intent(in) :: method_preconditioner
        logical, optional, intent(in) :: is_stationary
        logical, optional, intent(in) :: is_strict_mode
        type(FILL_LEVEL_USED), optional, intent(in) :: fill_level
        type(IterativeParams) :: params
        integer :: allocate_status

        allocate (params%x_init(N), STAT=allocate_status)
        if (allocate_status /= 0) stop "ERROR :: Unable to allocate x_init"
        if (present(x0)) then
            params%x_init = x0
        else
            params%x_init = 0.0_dp
        end if

        if (present(max_iter_choice)) params%max_iter = 1000
        if (present(epsi_tol)) params%tol = 1.0e-6_dp
        if (present(omega)) params%omega = omega
        if (present(alpha)) params%alpha = alpha
        if (present(Norm_choice)) params%norm = Norm_choice
        if (present(fill_level)) params%fill_level = fill_level

        allocate (params%residual(N), STAT=allocate_status)
        if (allocate_status /= 0) stop "ERROR :: Unable to allocate residual"

        select case (this%method_type%id)
        case (METHOD_SIP_ILU%id)
            allocate (params%L(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate L"
            allocate (params%U(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate U"
            call ILU_decomposition(A, params%L, params%U)
        case (METHOD_SIP_ICF%id)
            allocate (params%L(N, N), STAT=allocate_status)
            if (allocate_status /= 0) stop "ERROR :: Unable to allocate L"
            call Incomplete_Cholesky_decomposition(A, params%L)
        case (METHOD_CONJUGATE_GRADIENT%id)
            allocate (params%p(N), STAT=allocate_status)
            if (allocate_status /= 0) then
                stop "ERROR :: Unable to allocate optimal descent direction p"
            end if
        case (METHOD_CONJUGATE_RESIDUAL%id)
            allocate (params%p(N), STAT=allocate_status)
            if (allocate_status /= 0) then
                stop "ERROR :: Unable to allocate optimal descent direction p"
            end if
        case (METHOD_CGNR%id)
            allocate (params%p(N), STAT=allocate_status)
            if (allocate_status /= 0) then
                stop "ERROR :: Unable to allocate optimal descent direction p"
            end if
        case (METHOD_CGNE%id)
            allocate (params%p(N), STAT=allocate_status)
            if (allocate_status /= 0) then
                stop "ERROR :: Unable to allocate optimal descent direction p"
            end if
        end select

        if (present(method_preconditioner)) then
            params%precond => ApplyPreconditioner
            select case (method_preconditioner%id)
            case (METHOD_PRECOND_JACOBI%id)
                allocate (params%D(N, N), STAT=allocate_status)
                if (allocate_status /= 0) stop "ERROR :: Unable to allocate D"
                params%D = Calculate_Jacobi_preconditioner(A)
                this%preconditioner_type = METHOD_PRECOND_JACOBI
                this%requirements%needs_diag_dom = .true.
                this%requirements%needs_SPD = .true.
            case (METHOD_PRECOND_GS%id)
                allocate (params%L(N, N), STAT=allocate_status)
                if (allocate_status /= 0) stop "ERROR :: Unable to allocate L"
                params%L = Calculate_Gauss_Seidel_preconditioner(A)
                this%preconditioner_type = METHOD_PRECOND_GS
                this%requirements%needs_diag_dom = .true.
                this%requirements%needs_SPD = .true.
            case (METHOD_PRECOND_SOR%id)
                allocate (params%L(N, N), STAT=allocate_status)
                if (allocate_status /= 0) stop "ERROR :: Unable to allocate L"
                params%L = Calculate_SOR_preconditioner(A, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_SOR
                this%requirements%needs_diag_dom = .true.
                this%requirements%needs_SPD = .true.
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            case (METHOD_PRECOND_JOR%id)
                allocate (params%D(N, N), STAT=allocate_status)
                if (allocate_status /= 0) stop "ERROR :: Unable to allocate D"
                params%D = Calculate_JOR_preconditioner(A, params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_JOR
                this%requirements%needs_diag_dom = .true.
                this%requirements%needs_SPD = .true.
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            case (METHOD_PRECOND_ILU%id)
                allocate (params%L(N, N), STAT=allocate_status)
                if (allocate_status /= 0) stop "ERROR :: Unable to allocate L"
                allocate (params%U(N, N), STAT=allocate_status)
                if (allocate_status /= 0) stop "ERROR :: Unable to allocate U"
                if (params%fill_level%id /= FILL_LEVEL_NONE%id) then
                    call Calculate_ILU_preconditioner(A, params%L, params%U, &
                                                      params%omega, params%alpha, &
                                                      params%fill_level%id)
                else
                    call Calculate_ILU_preconditioner(A, params%L, params%U, &
                                                      params%omega, params%alpha)
                end if
                this%preconditioner_type = METHOD_PRECOND_ILU
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            case (METHOD_PRECOND_ICF%id)
                allocate (params%L(N, N), STAT=allocate_status)
                if (allocate_status /= 0) stop "ERROR :: Unable to allocate L"
                if (params%fill_level%id /= FILL_LEVEL_NONE%id) then
                    params%L = Calculate_ICF_preconditioner(A, &
                                                            params%omega, params%alpha, &
                                                            params%fill_level%id)
                else
                    params%L = Calculate_ICF_preconditioner(A, &
                                                            params%omega, params%alpha)
                end if
                this%preconditioner_type = METHOD_PRECOND_ICF
                this%requirements%needs_SPD = .true.
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            case (METHOD_PRECOND_SSOR%id)
                allocate (params%L(N, N), STAT=allocate_status)
                if (allocate_status /= 0) stop "ERROR :: Unable to allocate L"
                allocate (params%D(N, N), STAT=allocate_status)
                if (allocate_status /= 0) stop "ERROR :: Unable to allocate D"
                call Calculate_SSOR_preconditioner(A, params%L, params%D, &
                                                   params%omega, params%alpha)
                this%preconditioner_type = METHOD_PRECOND_SSOR
                this%requirements%needs_SPD = .true.
                this%relaxation_factor_preconditioner = RELAXATION_FACTOR_OMEGA
            case DEFAULT
                stop "ERROR :: Unknown method "
            end select
        end if

        if (present(is_stationary)) then
            if (is_stationary) then
                params%is_stationary = .true.
                this%requirements%needs_SPD = .false.
            else
                params%is_stationary = .false.
            end if
        end if

        if (present(is_strict_mode)) then
            if (is_strict_mode) then
                params%strict_mode = .true.
            else
                params%strict_mode = .false.
            end if
        end if

    end function Init_IterativeParams

    subroutine test_matrix(this, A, params, verbose)
        class(IterativeMethod), intent(inout) :: this
        real(dp), dimension(:, :), intent(in) :: A
        type(IterativeParams), intent(in) :: params
        type(Logger), optional, intent(inout) :: verbose
        character(KIND=ucs4, LEN=100) :: msg
        logical :: show_matrix_test

        show_matrix_test = .false.
        if (present(verbose)) then
            if (verbose%show_matrix_test) show_matrix_test = .true.
        end if

        if (show_matrix_test) then
            call verbose%write(center_with_fill("Testing matrix properties for method:"// &
                                                trim(this%method_type%name), &
                                                100, fill_char="="), box_style="top")
            call verbose%write(ucs4_"", box_style="middle")
        end if

        if (this%requirements%needs_square) then
            if (show_matrix_test) then
                msg = ucs4_"Checking if the matrix is square..."
                call verbose%log_info(msg)
            end if
            if (.not. is_square_matrix(A)) then
                write (msg, '(2A)') trim(this%method_type%name), " requires a square matrix"
                if (params%strict_mode) then
                    if (show_matrix_test) call verbose%log_error(msg)
                    stop
                else
                    if (show_matrix_test) call verbose%log_warning(msg)
                end if
            end if
        end if

        if (this%requirements%needs_SPD) then
            if (show_matrix_test) then
                msg = ucs4_"Checking if the matrix is symmetric positive definite (SPD)..."
                call verbose%log_info(msg)
            end if
            if (.not. is_SPD(A)) then
                write (msg, '(2A)') trim(this%method_type%name), &
                    " method requires a symmetric positive definite matrix."
                if (params%strict_mode) then
                    if (show_matrix_test) call verbose%log_error(msg)
                    stop
                else
                    if (show_matrix_test) call verbose%log_warning(msg)
                end if
            end if
        end if

        if (this%requirements%needs_diag_dom) then
            if (show_matrix_test) then
                msg = ucs4_"Checking if the matrix is diagonally dominant..."
                call verbose%log_info(msg)
            end if
            if (.not. is_diagonally_dominant(A)) then
                write (msg, '(2A)') trim(this%method_type%name), &
                    " method requires a diagonally dominant matrix."
                if (params%strict_mode) then
                    if (show_matrix_test) call verbose%log_error(msg)
                    stop
                else
                    if (show_matrix_test) call verbose%log_warning(msg)
                end if
            end if
        end if

        if (this%requirements%needs_symetric) then
            if (show_matrix_test) then
                msg = ucs4_"Checking if the matrix is symmetric..."
                call verbose%log_info(msg)
            end if
            if (.not. is_symmetric(A)) then
                write (msg, '(2A)') trim(this%method_type%name), &
                    " method requires a symmetric matrix."
                if (params%strict_mode) then
                    if (show_matrix_test) call verbose%log_error(msg)
                    stop
                else
                    if (show_matrix_test) call verbose%log_warning(msg)
                end if
            end if
        end if

        if (show_matrix_test) then
            call verbose%write(center_with_fill("", width=100, fill_char="="), box_style="bottom")
            call verbose%write(ucs4_"")
        end if

    end subroutine test_matrix

    subroutine Dealocate_IterativeParams(this, params, success)
        class(IterativeMethod), intent(inout) :: this
        type(IterativeParams), intent(inout) :: params
        logical, optional, intent(out) :: success
        integer :: deallocate_status

        if (present(success)) success = .true.

        if (allocated(params%x_init)) deallocate (params%x_init, STAT=deallocate_status)
        if (allocated(params%L)) deallocate (params%L, STAT=deallocate_status)
        if (allocated(params%U)) deallocate (params%U, STAT=deallocate_status)
        if (allocated(params%D)) deallocate (params%D, STAT=deallocate_status)
        if (allocated(params%residual)) deallocate (params%residual, STAT=deallocate_status)
        if (allocated(params%p)) deallocate (params%p, STAT=deallocate_status)
        params%norm = NORM_2
        params%fill_level = FILL_LEVEL_NONE

        if (deallocate_status /= 0 .and. present(success)) success = .false.

        this%preconditioner_type = METHOD_PRECOND_NONE
        this%method_type = METHOD_ITERATIVE_NONE
        this%relaxation_factor = RELAXATION_FACTOR_NONE
        this%relaxation_factor_preconditioner = RELAXATION_FACTOR_NONE

        this%requirements = IterativeMethodRequirements()

    end subroutine Dealocate_IterativeParams

    !======================
    ! Solve the system
    !======================

    function IterativeMethod_solve(this, A, b, params, verbose) result(x)

        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b
        type(IterativeParams), intent(inout) :: params
        type(Logger), optional, intent(inout) :: verbose
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: x0, x_new
        integer :: k, N, frequency
        integer :: start_system_clock, end_system_clock, rate
        real(dp) :: elapsed_time
        logical :: show_info_solver
        logical :: show_iteration
        logical :: show_final

        N = size(A, 1)
        x0 = params%x_init

        params%residual = b - matmul(A, x0)
        params%norm_initial_residual = params%norm_function(params%residual)

        show_info_solver = .false.
        show_iteration = .false.
        show_final = .false.
        frequency = 1
        if (present(verbose)) then
            if (verbose%show_info_solver) show_info_solver = .true.
            if (verbose%show_iteration) show_iteration = .true.
            if (verbose%show_final) show_final = .true.
            frequency = verbose%frequency
        end if

        if (show_info_solver) call log_solver_info(this, params, verbose, N)
        if (show_iteration) call log_iteration_header(verbose)

        call system_clock(start_system_clock, count_rate=rate)

        do k = 1, params%max_iter

            params%k = k
            if (k == params%max_iter) then
                exit
            end if

            x_new = this%solve_method(A, b, x0, params)

            params%residual = b - matmul(A, x_new)
            params%norm_residual = params%norm_function(params%residual)
            if (show_iteration .and. mod(k, frequency) == 0) then
                call system_clock(end_system_clock)
                elapsed_time = real(end_system_clock - start_system_clock, dp) / real(rate, dp)
                call log_iteration_step(verbose, k, params, elapsed_time)
            end if

            if (params%norm_residual < params%tol) exit

            x0 = x_new

        end do

        if (show_iteration) then
            call verbose%write(center_with_fill("", width=100, fill_char="="), box_style="bottom")
        end if

        if (show_final) then
            call system_clock(end_system_clock)
            elapsed_time = real(end_system_clock - start_system_clock, dp) / real(rate, dp)
            call log_final_result(verbose, k, params, x_new, elapsed_time, N)
        end if

        x = x_new

    end function IterativeMethod_solve

    !======================
    ! Log solver information
    !======================

    subroutine log_solver_info(solver, params, verbose, N)
        type(IterativeMethod), intent(in) :: solver
        type(IterativeParams), intent(in) :: params
        type(Logger), intent(inout) :: verbose
        integer, intent(in) :: N
        character(10) :: date, time
        character(KIND=ucs4, LEN=100) :: msg

        call verbose%write(ucs4_"")
        call verbose%write(center_with_fill("Starting system solver", width=100, fill_char="="), &
                           box_style="top")
        call verbose%write(ucs4_"", box_style="middle")

        ! call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date, TIME=time)
        write (date, '(A)') date(:4)//"-"//date(5:6)//"-"//date(7:8)
        write (time, '(A)') time(:2)//":"//time(3:4)//":"//time(5:6)
        call log_field(verbose, "Date and time", trim(date)//" "//trim(time))

        if (solver%relaxation_factor%id == RELAXATION_FACTOR_OMEGA%id) then
            call log_field(verbose, &
                           "Method used", &
                           trim(solver%method_type%name)//" "//trim(solver%method_type%name2))
            write (msg, '(A, T40, 3A, ES0.4)') "Relaxation factor used", &
                ": ", &
                trim(solver%relaxation_factor%name), &
                " = ", params%omega
            call verbose%log_info(msg)
        else if (solver%relaxation_factor%id == RELAXATION_FACTOR_ALPHA%id) then
            call log_field(verbose, "Method used", trim(solver%method_type%name))
            write (msg, '(A, T40, 3A, ES0.4)') "Relaxation factor used", &
                ": ", &
                trim(solver%relaxation_factor%name), &
                " = ", params%alpha
            call verbose%log_info(msg)
        else
            call log_field(verbose, "Method used", trim(solver%method_type%name))
        end if

        if (solver%preconditioner_type%id /= METHOD_PRECOND_NONE%id) then
            if (solver%relaxation_factor_preconditioner%id == RELAXATION_FACTOR_OMEGA%id) then
                call log_field(verbose, &
                               "Preconditioner used", &
                               trim(solver%preconditioner_type%name))
                write (msg, '(A, T40, 3A, ES0.4)') "Relaxation factor used", ": ", &
                    trim(solver%relaxation_factor_preconditioner%name), " = ", params%omega
                call verbose%log_info(msg)
            else
                call log_field(verbose, &
                               "Preconditioner used", &
                               trim(solver%preconditioner_type%name))
            end if
        end if

        if (solver%preconditioner_type%id == METHOD_PRECOND_ILU%id .or. &
            solver%preconditioner_type%id == METHOD_PRECOND_ICF%id) then
            if (params%fill_level%id /= FILL_LEVEL_NONE%id) then
                call log_field(verbose, "Fill level used", trim(params%fill_level%name))
            else
                call log_field(verbose, &
                               "Fill level used", &
                               "basic "//trim(solver%preconditioner_type%name))
            end if
        else if (solver%method_type%id == METHOD_SIP_ILU%id .or. &
                 solver%method_type%id == METHOD_SIP_ICF%id) then
            if (params%fill_level%id /= FILL_LEVEL_NONE%id) then
                call log_field(verbose, "Fill level used", trim(params%fill_level%name))
            else
                call log_field(verbose, &
                               "Fill level used", &
                               "basic "//trim(solver%method_type%name2))
            end if
        end if

        write (msg, '(A,T40,A,I0,A,I0)') "System size", ": ", N, " x ", N
        call verbose%log_info(msg)

        write (msg, '(A,T40,A,I0,A)') "Memory used", ": ", get_memory_kb(), " KB"
        call verbose%log_info(msg)

        call log_field(verbose, "Norme used", trim(params%norm%name))

        write (msg, '(A,T40,A, ES0.2)') "Convergence criterion (Tolerance)", &
            ": ||r|| < ", &
            params%tol
        call verbose%log_info(msg)

        call log_field(verbose, "Max iterations", params%max_iter)

        write (msg, '(A,T36,A,ES0.7)') "Initial residual norm", &
            ": ||r0|| = ", &
            params%norm_initial_residual
        call verbose%log_detail(msg)

        call verbose%write(center_with_fill("", width=100, fill_char="="), box_style="bottom")
        call verbose%write(ucs4_"")

    end subroutine log_solver_info

    subroutine log_iteration_header(verbose)
        class(Logger), intent(inout) :: verbose
        character(KIND=ucs4, LEN=100) :: msg

        call verbose%write(center_with_fill("Iterations", width=100, fill_char="="), &
                           box_style="top")
        call verbose%write(ucs4_"", box_style="middle")

        write (msg, '(A,5X,A,5X,A,5X,A)') "Iter", &
            "Residual norm (||r||)", &
            "Relative residual norm (||r||/||r0||)", &
            "Time (s)"
        call verbose%log_info(msg)
        call verbose%write(repeat(ucs4_"-", 100), box_style="middle")
    end subroutine log_iteration_header

    subroutine log_iteration_step(verbose, k, params, elapsed_time)
        class(Logger), intent(inout) :: verbose
        integer, intent(in) :: k
        real(dp), intent(in) :: elapsed_time
        type(IterativeParams), intent(in) :: params
        character(KIND=ucs4, LEN=100) :: msg
        integer :: end_system_clock

        call system_clock(end_system_clock)
        write (msg, '(T2,I0,T15,ES0.7,T48,ES0.7,T79,ES0.7)') k, &
            params%norm_residual, &
            params%norm_residual / params%norm_initial_residual, &
            elapsed_time
        call verbose%log_time(msg)
    end subroutine log_iteration_step

    subroutine log_final_result(verbose, k, params, x_new, elapsed_time, N)
        class(Logger), intent(inout) :: verbose
        integer, intent(in) :: N
        integer, intent(in) :: k
        real(dp), intent(in) :: elapsed_time
        type(IterativeParams), intent(in) :: params
        real(dp), dimension(:), intent(in) :: x_new
        character(KIND=ucs4, LEN=100) :: msg
        integer :: end_system_clock
        integer :: i

        call verbose%write(ucs4_"")
        call verbose%write(center_with_fill("Results", width=100, fill_char="="), box_style="top")
        call verbose%write(ucs4_" ", box_style="middle")

        if (k < params%max_iter) then
            call log_field(verbose, "Status", "CONVERGED")

            write (msg, '(A,T40,A,ES0.7,A,ES0.1,A)') "Final residual", &
                ": ||r|| = ", &
                params%norm_residual, &
                " < ", &
                params%tol, &
                " (convergence achieved)"
            call verbose%log_info(msg)

            write (msg, '(A,T36,A,ES0.7)') "Relative residual norm", ": ||r||/||r0|| = ", &
                params%norm_residual / params%norm_initial_residual
            call verbose%log_detail(msg)

            call log_field(verbose, "Total iterations", k)

            write (msg, '(A,T40,A,ES0.7)') "Solution", ": x = ["
            if (N < 6) then
                do i = 1, size(x_new)
                    write (msg, '(2A,ES0.7)') trim(msg), " ", x_new(i)
                end do
            else
                write (msg, '(2A,ES0.7,3X,ES0.7,3X,ES0.7,3X,A,ES0.7,3X,ES0.7,3X,ES0.7,A)') &
                    trim(msg), &
                    " ", &
                    x_new(1), &
                    x_new(2), &
                    x_new(3), &
                    x_new(N - 3), &
                    x_new(N - 2), &
                    x_new(N - 1)
            end if
            write (msg, '(2A,ES0.7)') trim(msg), "]"
            call verbose%log_info(msg)

            call system_clock(end_system_clock)
            write (msg, '(A,T40,A,ES0.7,A,I0)') "Solver completed in ", ": ", &
                elapsed_time, " seconds"
            call verbose%log_info(msg)

            write (msg, '(A,T40,A,I0,A)') "Memory used", ": ", get_memory_kb(), " KB"
            call verbose%log_info(msg)
        else
            call log_field(verbose, "Status", "NOT CONVERGED")

            write (msg, '(A,T40,A,ES0.7,A,ES0.1,A)') "Final residual", &
                ": ||r|| = ", &
                params%norm_residual, &
                " >= ", params%tol, " (convergence not achieved)"
            call verbose%log_info(msg)
        end if
        call verbose%write(center_with_fill("", width=100, fill_char="="), box_style="bottom")
    end subroutine log_final_result

    !======================
    ! Solve methods
    !======================

    !> Jacobi iterative method
    !>
    !> This subroutine implements the Jacobi method for solving linear systems.
    function solve_Jacobi(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        integer :: i, N

        N = size(A, 1)

        ! forward
        do i = 1, N
            x(i) = b(i) - &
                   dot_product(A(i, 1:i - 1), x0(1:i - 1)) - &
                   dot_product(A(i, i + 1:N), x0(i + 1:N))
            x(i) = x(i) / A(i, i)
        end do

    end function solve_Jacobi

    !> Gauss-Seidel iterative method
    !>
    !> This subroutine implements the Gauss-Seidel method for solving linear systems.
    function solve_Gauss_Seidel(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        integer :: i, N

        N = size(A, 1)

        ! forward
        do i = 1, N
            x(i) = b(i) - &
                   dot_product(A(i, 1:i - 1), x(1:i - 1)) - &
                   dot_product(A(i, i + 1:N), x0(i + 1:N))
            x(i) = x(i) / A(i, i)
        end do

    end function solve_Gauss_Seidel

    !> Successive Over-Relaxation (SOR) iterative method
    !>
    !> This subroutine implements the SOR method for solving linear systems.
    function solve_SOR(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        integer :: i, N

        N = size(A, 1)

        ! forward
        do i = 1, N
            x(i) = b(i) - &
                   dot_product(A(i, 1:i - 1), x(1:i - 1)) - &
                   dot_product(A(i, i + 1:N), x0(i + 1:N))
            x(i) = params%omega * (x(i) / A(i, i) - x0(i)) + x0(i)
        end do

    end function solve_SOR

    !> Jacobi over-relaxation (JOR) iterative method
    !>
    !> This subroutine implements the Jacobi over-relaxation method for solving linear systems.
    function solve_JOR(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        integer :: i, N

        N = size(A, 1)

        ! forward
        do i = 1, N
            x(i) = b(i) - &
                   dot_product(A(i, 1:i - 1), x0(1:i - 1)) - &
                   dot_product(A(i, i + 1:N), x0(i + 1:N))
            x(i) = params%omega * (x(i) / A(i, i) - x0(i)) + x0(i)
        end do

    end function solve_JOR

    !> strongly implicit procedure (SIP) method (or stone's method)
    !>
    !> This subroutine implements the SIP method for solving linear systems.
    !> It uses the incomplete LU decomposition of the matrix A.
    function solve_SIP_ILU(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: y, z

        if (.not. allocated(params%L) .or. .not. allocated(params%U)) then
            stop "ERROR :: Incomplete LU decomposition not initialized"
        end if

        y = forward(params%L, params%residual)

        z = backward(params%U, y)

        x = x0 + params%omega * z

    end function solve_SIP_ILU

    !> strongly implicit procedure (SIP) method (or stone's method)
    !>
    !> This subroutine implements the SIP method for solving linear systems.
    !> It uses the incomplete Cholesky decomposition of the matrix A.
    function solve_SIP_ICF(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: y, z

        if (.not. allocated(params%L)) then
            stop "ERROR :: Incomplete LU decomposition not initialized"
        end if

        y = forward(params%L, params%residual)

        z = backward(transpose(params%L), y)

        x = x0 + params%omega * z

    end function solve_SIP_ICF

    !> Symmetric successive Over-Relaxation (SSOR) iterative method
    !>
    !> This subroutine implements the SSOR method for solving linear systems.
    function solve_SSOR(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: x_tmp
        integer :: i, N

        N = size(A, 1)

        ! forward
        do i = 1, N
            x_tmp(i) = b(i) - &
                       dot_product(A(i, 1:i - 1), x_tmp(1:i - 1)) - &
                       dot_product(A(i, i + 1:N), x0(i + 1:N))
            x_tmp(i) = params%omega * (x_tmp(i) / A(i, i) - x0(i)) + x0(i)
        end do

        ! backward
        do i = N, 1, -1
            x(i) = b(i) - &
                   dot_product(A(i, 1:i - 1), x_tmp(1:i - 1)) - &
                   dot_product(A(i, i + 1:N), x(i + 1:N))
            x(i) = params%omega * (x(i) / A(i, i) - x_tmp(i)) + x_tmp(i)
        end do

    end function solve_SSOR

    !> Richardson iterative method
    !>
    !> This subroutine implements the Richardson method for solving linear systems.
    function solve_Richardson(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: z_prec

        if (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) then
            if (.not. params%is_stationary) then
                params%alpha = dot_product(params%residual, params%residual) / &
                               dot_product(params%residual, matmul(A, params%residual))
            end if
            x = x0 + params%alpha * params%residual
        else if (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) then
            z_prec = params%precond(this%preconditioner_type, params%residual)
            if (.not. params%is_stationary) then
                params%alpha = dot_product(params%residual, z_prec) / &
                               dot_product(z_prec, matmul(A, z_prec))
            end if
            x = x0 + params%alpha * z_prec
        end if

    end function solve_Richardson

    !> Conjugate Gradient iterative method
    !>
    !> This subroutine implements the Conjugate Gradient method for solving linear systems.
    function solve_ConjugateGradient(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: z_prec

        if (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) then
            if (params%k == 1) then
                params%p = params%residual
            else if (params%k /= 1) then
                params%beta = dot_product(params%residual, params%residual) / &
                              params%old_dot_product
                params%p = params%residual + params%beta * params%p
            end if

            params%alpha = dot_product(params%residual, params%residual) / &
                           dot_product(params%p, matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(params%residual, params%residual)
        else if (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) then
            if (this%preconditioner_type%id == METHOD_PRECOND_GS%id .or. &
                this%preconditioner_type%id == METHOD_PRECOND_SOR%id) then
                stop "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            end if

            if (params%k == 1) then
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%p = z_prec
            else if (params%k /= 1) then
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%beta = dot_product(z_prec, params%residual) / params%old_dot_product
                params%p = z_prec + params%beta * params%p
            end if

            params%alpha = dot_product(z_prec, params%residual) / &
                           dot_product(params%p, matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(z_prec, params%residual)
        end if

    end function solve_ConjugateGradient

    !> Conjugate Residual iterative method
    !>
    !> This subroutine implements the Conjugate Residual method for solving linear systems.
    function solve_ConjugateResidual(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: z_prec

        if (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) then
            if (params%k == 1) then
                params%p = params%residual
            else if (params%k /= 1) then
                params%beta = dot_product(params%residual, matmul(A, params%residual)) / &
                              params%old_dot_product
                params%p = params%residual + params%beta * params%p
            end if

            params%alpha = dot_product(params%residual, matmul(A, params%residual)) / &
                           dot_product(matmul(A, params%p), matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(params%residual, matmul(A, params%residual))
        else if (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) then
            if (this%preconditioner_type%id == METHOD_PRECOND_GS%id .or. &
                this%preconditioner_type%id == METHOD_PRECOND_SOR%id) then
                stop "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            end if

            if (params%k == 1) then
                params%p = params%precond(this%preconditioner_type, params%residual)
                z_prec = params%p
            else if (params%k /= 1) then
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%beta = dot_product(z_prec, matmul(A, params%residual)) / &
                params%old_dot_product
                params%p = z_prec + params%beta * params%p
            end if

            params%alpha = dot_product(z_prec, matmul(A, params%residual)) / &
                           dot_product(matmul(A, params%p), matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(z_prec, matmul(A, params%residual))
        end if

    end function solve_ConjugateResidual

    !> Conjugate Gradient on Normal Equations iterative method
    !>
    !> This subroutine implements the Conjugate Gradient on Normal Equations method (or Craig’s Method) for solving linear systems.
    function solve_CGNR(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: z_prec
        real(dp), dimension(size(A, 1)) :: AT_r

        if (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) then
            if (params%k == 1) then
                AT_r = matmul(transpose(A), params%residual)
                params%p = AT_r
            else if (params%k /= 1) then
                AT_r = matmul(transpose(A), params%residual)
                params%beta = dot_product(AT_r, AT_r) / params%old_dot_product
                params%p = AT_r + params%beta * params%p
            end if

            params%alpha = dot_product(AT_r, AT_r) / &
                           dot_product(matmul(A, params%p), matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(AT_r, AT_r)
        else if (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) then
            if (this%preconditioner_type%id == METHOD_PRECOND_GS%id .or. &
                this%preconditioner_type%id == METHOD_PRECOND_SOR%id) then
                stop "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            end if

            if (params%k == 1) then
                AT_r = matmul(transpose(A), params%residual)
                z_prec = params%precond(this%preconditioner_type, AT_r)
                params%p = z_prec
            else if (params%k /= 1) then
                AT_r = matmul(transpose(A), params%residual)
                z_prec = params%precond(this%preconditioner_type, AT_r)
                params%beta = dot_product(z_prec, AT_r) / params%old_dot_product
                params%p = z_prec + params%beta * params%p
            end if

            params%alpha = dot_product(z_prec, AT_r) / &
                           dot_product(matmul(A, params%p), matmul(A, params%p))

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(z_prec, AT_r)
        end if

    end function solve_CGNR

    !> Conjugate Gradient on Normal Residual iterative method
    !>
    !> This subroutine implements the Conjugate Gradient on Normal Residual method for solving linear systems.
    function solve_CGNE(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x
        real(dp), dimension(size(A, 1)) :: z_prec

        if (this%preconditioner_type%id == METHOD_PRECOND_NONE%id) then
            if (params%k == 1) then
                params%p = matmul(transpose(A), params%residual)
            else if (params%k /= 1) then
                params%beta = dot_product(params%residual, params%residual) / params%old_dot_product
                params%p = matmul(transpose(A), params%residual) + params%beta * params%p
            end if

            params%alpha = dot_product(params%residual, params%residual) / &
                           dot_product(params%p, params%p)

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(params%residual, params%residual)
        else if (this%preconditioner_type%id /= METHOD_PRECOND_NONE%id) then
            if (this%preconditioner_type%id == METHOD_PRECOND_GS%id .or. &
                this%preconditioner_type%id == METHOD_PRECOND_SOR%id) then
                stop "ERROR :: Preconditioner Gauss-Seidel and SOR not supported for Conjugate Gradient method"
            end if

            if (params%k == 1) then
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%p = matmul(transpose(A), z_prec)
            else if (params%k /= 1) then
                z_prec = params%precond(this%preconditioner_type, params%residual)
                params%beta = dot_product(z_prec, params%residual) / params%old_dot_product
                params%p = matmul(transpose(A), z_prec) + params%beta * params%p
            end if

            params%alpha = dot_product(z_prec, params%residual) / &
                           dot_product(params%p, params%p)

            x = x0 + params%alpha * params%p

            params%old_dot_product = dot_product(z_prec, params%residual)
        end if

    end function solve_CGNE

    function solve_GMRES(this, A, b, x0, params) result(x)
        class(IterativeMethod), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: b, x0
        type(IterativeParams), intent(inout) :: params
        real(dp), dimension(size(A, 1)) :: x

        ! in progress

    end function solve_GMRES

end module NAFPack_Iterative_methods
