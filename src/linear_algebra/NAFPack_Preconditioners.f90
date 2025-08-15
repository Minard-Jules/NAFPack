MODULE NAFPack_Preconditioners

    USE NAFPack_constant
    USE NAFPack_Iterative_types
    USE NAFPack_matricielle
    USE NAFPack_matrix_decomposition

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE

    PUBLIC :: MethodPreconditioner
    PUBLIC :: METHOD_PRECOND_NONE
    PUBLIC :: METHOD_PRECOND_JACOBI, METHOD_PRECOND_JOR
    PUBLIC :: METHOD_PRECOND_GS, METHOD_PRECOND_SOR, METHOD_PRECOND_SSOR
    PUBLIC :: METHOD_PRECOND_ILU, METHOD_PRECOND_ICF

    PUBLIC :: FILL_LEVEL_USED
    PUBLIC :: FILL_LEVEL_NONE
    PUBLIC :: FILL_LEVEL_0, FILL_LEVEL_1, FILL_LEVEL_2, FILL_LEVEL_3
    PUBLIC :: FILL_LEVEL_N

    PUBLIC :: Calculate_Jacobi_preconditioner
    PUBLIC :: Calculate_Gauss_Seidel_preconditioner
    PUBLIC :: Calculate_SOR_preconditioner
    PUBLIC :: Calculate_JOR_preconditioner
    PUBLIC :: Calculate_ILU_preconditioner
    PUBLIC :: Calculate_ICF_preconditioner
    PUBLIC :: Calculate_SSOR_preconditioner

    TYPE :: MethodPreconditioner
        INTEGER :: id
        CHARACTER(LEN=64) :: name
    END TYPE MethodPreconditioner

    TYPE :: Fill_level_used
        INTEGER :: id
        CHARACTER(LEN=64) :: name
        INTEGER :: VALUE
    END TYPE Fill_level_used

    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_NONE = MethodPreconditioner(0, "None")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_JACOBI = MethodPreconditioner(1, "Jacobi")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_GS = MethodPreconditioner(2, "Gauss-Seidel")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_SOR = MethodPreconditioner(3, "Successive Over-Relaxation")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_JOR = MethodPreconditioner(4, "Jacobi Over-Relaxation")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_ILU = MethodPreconditioner(5, "ILU")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_ICF = MethodPreconditioner(6, "ICF")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_SSOR = MethodPreconditioner(7, "SSOR")

    TYPE(Fill_level_used), PARAMETER :: FILL_LEVEL_NONE = Fill_level_used(-1, "None", -huge(1))
    TYPE(Fill_level_used), PARAMETER :: FILL_LEVEL_0 = Fill_level_used(0, "Level 0", 0)
    TYPE(Fill_level_used), PARAMETER :: FILL_LEVEL_1 = Fill_level_used(1, "Level 1", 1)
    TYPE(Fill_level_used), PARAMETER :: FILL_LEVEL_2 = Fill_level_used(2, "Level 2", 2)
    TYPE(Fill_level_used), PARAMETER :: FILL_LEVEL_3 = Fill_level_used(3, "Level 3", 3)
    TYPE(Fill_level_used) :: FILL_LEVEL_N = Fill_level_used(3, "Level N", 0)

CONTAINS

    FUNCTION Calculate_Jacobi_preconditioner(A) RESULT(D)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(size(A, 1), size(A, 2)) :: D
        INTEGER :: N, i

        N = size(A, 1)

        D = 0.d0

        IF (any(Diag(A) < epsi)) STOP "ERROR :: Zero diagonal in Jacobi preconditioner"
        FORALL (i=1:N) D(i, i) = 1.d0 / A(i, i)

    END FUNCTION Calculate_Jacobi_preconditioner

    FUNCTION Calculate_Gauss_Seidel_preconditioner(A) RESULT(L)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(size(A, 1), size(A, 2)) :: L
        INTEGER :: N, i, j

        N = size(A, 1)

        L = 0.d0

        IF (any(Diag(A) < epsi)) STOP "ERROR :: Zero diagonal in Gauss-Seidel preconditioner"
        FORALL (i=1:size(A, 1), j=1:size(A, 2), i >= j) L(i, j) = A(i, j)

    END FUNCTION Calculate_Gauss_Seidel_preconditioner

    FUNCTION Calculate_SOR_preconditioner(A, omega, alpha) RESULT(L)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(size(A, 1), size(A, 2)) :: L
        INTEGER :: N, i

        N = size(A, 1)

        L = 0.d0

        IF (any(Diag(A) < epsi)) STOP "ERROR :: Zero diagonal in SOR preconditioner"
        DO i = 1, size(A, 1)
            L(i, i) = 1.d0 / omega * A(i, i)
            L(i, 1:i - 1) = A(i, 1:i - 1)
        END DO

        L = alpha * L

    END FUNCTION Calculate_SOR_preconditioner

    FUNCTION Calculate_JOR_preconditioner(A, omega, alpha) RESULT(D)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(size(A, 1), size(A, 2)) :: D
        INTEGER :: N, i

        N = size(A, 1)

        D = 0.d0

        IF (any(Diag(A) < epsi)) STOP "ERROR :: Zero diagonal in JOR preconditioner"
        FORALL (i=1:size(A, 1)) D(i, i) = omega / A(i, i)

        D = D / alpha

    END FUNCTION Calculate_JOR_preconditioner

    SUBROUTINE Calculate_ILU_preconditioner(A, L, U, omega, alpha, fill_level)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(size(A, 1), size(A, 2)), INTENT(OUT) :: L, U
        INTEGER, OPTIONAL, INTENT(IN) :: fill_level
        INTEGER :: N

        N = size(A, 1)

        L = 0.d0
        U = 0.d0

        IF (present(fill_level)) THEN
            CALL ILU_decomposition(A, L, U, fill_level)
        ELSE
            CALL ILU_decomposition(A, L, U)
        END IF

        L = alpha / omega * L

    END SUBROUTINE Calculate_ILU_preconditioner

    FUNCTION Calculate_ICF_preconditioner(A, omega, alpha, fill_level) RESULT(L)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(size(A, 1), size(A, 2)) :: L
        INTEGER, OPTIONAL, INTENT(IN) :: fill_level
        INTEGER :: N

        N = size(A, 1)

        L = 0.d0

        IF (present(fill_level)) THEN
            CALL Incomplete_Cholesky_decomposition(A, L, fill_level)
        ELSE
            CALL Incomplete_Cholesky_decomposition(A, L)
        END IF

        L = alpha / omega * L

    END FUNCTION Calculate_ICF_preconditioner

    SUBROUTINE Calculate_SSOR_preconditioner(A, L, D, omega, alpha)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(size(A, 1), size(A, 2)), INTENT(OUT) :: L, D
        INTEGER :: N, i

        N = size(A, 1)

        L = 0.d0
        D = 0.d0

        DO i = 1, size(A, 1)
            L(i, i) = 1.d0 / omega * A(i, i)
            L(i, 1:i - 1) = A(i, 1:i - 1)

            D(i, i) = A(i, i)
        END DO

        L = (alpha * omega) / (2 - omega) * L

    END SUBROUTINE Calculate_SSOR_preconditioner

END MODULE NAFPack_Preconditioners
