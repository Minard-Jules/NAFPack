MODULE NAFPack_Preconditioners

    USE NAFPack_constant
    USE NAFPack_Iterative_types
    USE NAFPack_matricielle
    USE NAFPack_matrix_decomposition
    
    IMPLICIT NONE

    PRIVATE

    PUBLIC :: MethodPreconditioner
    PUBLIC :: METHOD_PRECOND_NONE
    PUBLIC :: METHOD_PRECOND_JACOBI, METHOD_PRECOND_JOR
    PUBLIC :: METHOD_PRECOND_GS, METHOD_PRECOND_SOR, METHOD_PRECOND_SSOR
    PUBLIC :: METHOD_PRECOND_ILU, METHOD_PRECOND_ICF

    PUBLIC :: Calculate_Jacobi_preconditioner
    PUBLIC :: Calculate_Gauss_Seidel_preconditioner
    PUBLIC :: Calculate_SOR_preconditioner
    PUBLIC :: Calculate_JOR_preconditioner
    PUBLIC :: Calculate_ILU_preconditioner
    PUBLIC :: Calculate_ICF_preconditioner
    PUBLIC :: Calculate_SSOR_preconditioner

    TYPE :: MethodPreconditioner
        INTEGER :: value
        CHARACTER(LEN=64) :: name
    END TYPE MethodPreconditioner

    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_NONE = MethodPreconditioner(0, "None")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_JACOBI = MethodPreconditioner(1, "Jacobi")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_GS = MethodPreconditioner(2, "Gauss-Seidel")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_SOR = MethodPreconditioner(3, "Successive Over-Relaxation")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_JOR = MethodPreconditioner(4, "Jacobi Over-Relaxation")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_ILU = MethodPreconditioner(5, "ILU")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_ICF = MethodPreconditioner(6, "ICF")
    TYPE(MethodPreconditioner), PARAMETER :: METHOD_PRECOND_SSOR = MethodPreconditioner(7, "SSOR")

    CONTAINS

    FUNCTION Calculate_Jacobi_preconditioner(A) RESULT(D)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: D
        INTEGER :: N, i

        N = SIZE(A, 1)

        D = 0.d0

        IF (ANY(Diag(A) < epsi)) STOP "ERROR :: Zero diagonal in Jacobi preconditioner"
        FORALL(i = 1:N) D(i, i) = 1.d0 / A(i, i)

    END FUNCTION Calculate_Jacobi_preconditioner

    FUNCTION Calculate_Gauss_Seidel_preconditioner(A) RESULT(L)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: L
        INTEGER :: N, i, j

        N = SIZE(A, 1)

        L = 0.d0

        IF (ANY(Diag(A) < epsi)) STOP "ERROR :: Zero diagonal in Gauss-Seidel preconditioner"
        FORALL(i=1:SIZE(A,1), j=1:SIZE(A,2), i >= j) L(i, j) = A(i, j)

    END FUNCTION Calculate_Gauss_Seidel_preconditioner

    FUNCTION Calculate_SOR_preconditioner(A, omega, alpha) RESULT(L)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: L
        INTEGER :: N, i

        N = SIZE(A, 1)

        L = 0.d0

        IF (ANY(Diag(A) < epsi)) STOP "ERROR :: Zero diagonal in SOR preconditioner"
        DO i = 1, SIZE(A,1)
            L(i,i) = 1.d0/omega * A(i,i)
            L(i,1:i-1) = A(i,1:i-1)
        END DO

        L = alpha * L

    END FUNCTION Calculate_SOR_preconditioner

    FUNCTION Calculate_JOR_preconditioner(A, omega, alpha) RESULT(D)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: D
        INTEGER :: N, i

        N = SIZE(A, 1)

        D = 0.d0

        IF (ANY(Diag(A) < epsi)) STOP "ERROR :: Zero diagonal in JOR preconditioner"
        FORALL(i=1:SIZE(A, 1)) D(i,i) = omega / A(i,i)

        D = D / alpha

    END FUNCTION Calculate_JOR_preconditioner

    SUBROUTINE Calculate_ILU_preconditioner(A, L, U, omega, alpha)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)), INTENT(OUT) :: L, U
        INTEGER :: N

        N = SIZE(A, 1)

        L = 0.d0
        U = 0.d0

        CALL ILU_decomposition(A, L, U)
        
        L = alpha / omega * L

    END SUBROUTINE Calculate_ILU_preconditioner

    FUNCTION Calculate_ICF_preconditioner(A, omega, alpha) RESULT(L)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: L
        INTEGER :: N

        N = SIZE(A, 1)

        L = 0.d0

        CALL Incomplete_Cholesky_decomposition(A, L)

        L = alpha / omega * L

    END FUNCTION Calculate_ICF_preconditioner

    SUBROUTINE Calculate_SSOR_preconditioner(A, L, D, U, omega, alpha)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), INTENT(IN) :: omega, alpha
        REAL(dp), DIMENSION(SIZE(A, 1), SIZE(A, 2)), INTENT(OUT) :: L, D, U
        INTEGER :: N, i

        N = SIZE(A, 1)

        L = 0.d0
        D = 0.d0
        U = 0.d0

        DO i = 1, SIZE(A,1)
            L(i,i) = 1.d0/omega * A(i,i)
            L(i,1:i-1) = A(i,1:i-1)
            
            D(i,i) = A(i,i)

            U(i,i) = 1.d0/omega * A(i,i)
            U(i,i+1:) = A(i,i+1:)
        END DO

        L = (alpha * omega)/(2-omega) * L

    END SUBROUTINE Calculate_SSOR_preconditioner

END MODULE NAFPack_Preconditioners
