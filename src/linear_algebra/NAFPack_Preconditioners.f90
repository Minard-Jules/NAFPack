module NAFPack_Preconditioners

    use NAFPack_kinds, only: dp
    USE NAFPack_constant, only: TOL_CONVERGENCE_dp
    use NAFPack_matricielle, only: Diag
    use NAFPack_matrix_decomposition, only: Incomplete_Cholesky_decomposition, ILU_decomposition

    implicit none(type, external)

    private

    public :: MethodPreconditioner
    public :: METHOD_PRECOND_NONE
    public :: METHOD_PRECOND_JACOBI, METHOD_PRECOND_JOR
    public :: METHOD_PRECOND_GS, METHOD_PRECOND_SOR, METHOD_PRECOND_SSOR
    public :: METHOD_PRECOND_ILU, METHOD_PRECOND_ICF

    public :: FILL_LEVEL_USED
    public :: FILL_LEVEL_NONE
    public :: FILL_LEVEL_0, FILL_LEVEL_1, FILL_LEVEL_2, FILL_LEVEL_3
    public :: FILL_LEVEL_N

    public :: Calculate_Jacobi_preconditioner
    public :: Calculate_Gauss_Seidel_preconditioner
    public :: Calculate_SOR_preconditioner
    public :: Calculate_JOR_preconditioner
    public :: Calculate_ILU_preconditioner
    public :: Calculate_ICF_preconditioner
    public :: Calculate_SSOR_preconditioner

    type :: MethodPreconditioner
        integer :: id
        character(LEN=64) :: name
    end type MethodPreconditioner

    type :: Fill_level_used
        integer :: id
        character(LEN=64) :: name
        integer :: value
    end type Fill_level_used

    type(MethodPreconditioner), parameter :: METHOD_PRECOND_NONE = &
                                             MethodPreconditioner(0, "None")
    type(MethodPreconditioner), parameter :: METHOD_PRECOND_JACOBI = &
                                             MethodPreconditioner(1, "Jacobi")
    type(MethodPreconditioner), parameter :: METHOD_PRECOND_GS = &
                                             MethodPreconditioner(2, "Gauss-Seidel")
    type(MethodPreconditioner), parameter :: METHOD_PRECOND_SOR = &
                                             MethodPreconditioner(3, "Successive Over-Relaxation")
    type(MethodPreconditioner), parameter :: METHOD_PRECOND_JOR = &
                                             MethodPreconditioner(4, "Jacobi Over-Relaxation")
    type(MethodPreconditioner), parameter :: METHOD_PRECOND_ILU = &
                                             MethodPreconditioner(5, "ILU")
    type(MethodPreconditioner), parameter :: METHOD_PRECOND_ICF = &
                                             MethodPreconditioner(6, "ICF")
    type(MethodPreconditioner), parameter :: METHOD_PRECOND_SSOR = &
                                             MethodPreconditioner(7, "SSOR")

    type(Fill_level_used), parameter :: FILL_LEVEL_NONE = Fill_level_used(-1, "None", -huge(1))
    type(Fill_level_used), parameter :: FILL_LEVEL_0 = Fill_level_used(0, "Level 0", 0)
    type(Fill_level_used), parameter :: FILL_LEVEL_1 = Fill_level_used(1, "Level 1", 1)
    type(Fill_level_used), parameter :: FILL_LEVEL_2 = Fill_level_used(2, "Level 2", 2)
    type(Fill_level_used), parameter :: FILL_LEVEL_3 = Fill_level_used(3, "Level 3", 3)
    type(Fill_level_used) :: FILL_LEVEL_N = Fill_level_used(3, "Level N", 0)

contains

    function Calculate_Jacobi_preconditioner(A) result(D)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 2)) :: D
        integer :: N, i

        N = size(A, 1)

        D = 0.d0

        if (any(Diag(A) < TOL_CONVERGENCE_dp)) stop "ERROR :: Zero diagonal in Jacobi preconditioner"
        forall (i=1:N) D(i, i) = 1.d0 / A(i, i)

    end function Calculate_Jacobi_preconditioner

    function Calculate_Gauss_Seidel_preconditioner(A) result(L)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 2)) :: L
        integer :: N, i, j

        N = size(A, 1)

        L = 0.d0

        if (any(Diag(A) < TOL_CONVERGENCE_dp)) stop "ERROR :: Zero diagonal in Gauss-Seidel preconditioner"
        forall (i=1:size(A, 1), j=1:size(A, 2), i >= j) L(i, j) = A(i, j)

    end function Calculate_Gauss_Seidel_preconditioner

    function Calculate_SOR_preconditioner(A, omega, alpha) result(L)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), intent(in) :: omega, alpha
        real(dp), dimension(size(A, 1), size(A, 2)) :: L
        integer :: N, i

        N = size(A, 1)

        L = 0.d0

        if (any(Diag(A) < TOL_CONVERGENCE_dp)) stop "ERROR :: Zero diagonal in SOR preconditioner"
        do i = 1, size(A, 1)
            L(i, i) = 1.d0 / omega * A(i, i)
            L(i, 1:i - 1) = A(i, 1:i - 1)
        end do

        L = alpha * L

    end function Calculate_SOR_preconditioner

    function Calculate_JOR_preconditioner(A, omega, alpha) result(D)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), intent(in) :: omega, alpha
        real(dp), dimension(size(A, 1), size(A, 2)) :: D
        integer :: N, i

        N = size(A, 1)

        D = 0.d0

        if (any(Diag(A) < TOL_CONVERGENCE_dp)) stop "ERROR :: Zero diagonal in JOR preconditioner"
        forall (i=1:size(A, 1)) D(i, i) = omega / A(i, i)

        D = D / alpha

    end function Calculate_JOR_preconditioner

    subroutine Calculate_ILU_preconditioner(A, L, U, omega, alpha, fill_level)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), intent(in) :: omega, alpha
        real(dp), dimension(size(A, 1), size(A, 2)), intent(out) :: L, U
        integer, optional, intent(in) :: fill_level
        integer :: N

        N = size(A, 1)

        L = 0.d0
        U = 0.d0

        if (present(fill_level)) then
            call ILU_decomposition(A, L, U, fill_level)
        else
            call ILU_decomposition(A, L, U)
        end if

        L = alpha / omega * L

    end subroutine Calculate_ILU_preconditioner

    function Calculate_ICF_preconditioner(A, omega, alpha, fill_level) result(L)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), intent(in) :: omega, alpha
        real(dp), dimension(size(A, 1), size(A, 2)) :: L
        integer, optional, intent(in) :: fill_level
        integer :: N

        N = size(A, 1)

        L = 0.d0

        if (present(fill_level)) then
            call Incomplete_Cholesky_decomposition(A, L, fill_level)
        else
            call Incomplete_Cholesky_decomposition(A, L)
        end if

        L = alpha / omega * L

    end function Calculate_ICF_preconditioner

    subroutine Calculate_SSOR_preconditioner(A, L, D, omega, alpha)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), intent(in) :: omega, alpha
        real(dp), dimension(size(A, 1), size(A, 2)), intent(out) :: L, D
        integer :: N, i

        N = size(A, 1)

        L = 0.d0
        D = 0.d0

        do i = 1, size(A, 1)
            L(i, i) = 1.d0 / omega * A(i, i)
            L(i, 1:i - 1) = A(i, 1:i - 1)

            D(i, i) = A(i, i)
        end do

        L = (alpha * omega) / (2 - omega) * L

    end subroutine Calculate_SSOR_preconditioner

end module NAFPack_Preconditioners
