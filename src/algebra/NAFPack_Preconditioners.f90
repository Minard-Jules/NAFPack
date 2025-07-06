MODULE NAFPack_Preconditioners

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    
    IMPLICIT NONE

    PRIVATE

    PUBLIC :: preconditioner, preconditioner_use
    PUBLIC :: ILU_preconditioner

    TYPE :: preconditioner
        LOGICAL :: Jacobi = .FALSE.
        LOGICAL :: Gauss_Seidel = .FALSE.
        LOGICAL :: SOR = .FALSE.
        LOGICAL :: ILU = .FALSE.
        LOGICAL :: ICF = .FALSE.
        LOGICAL :: SSOR = .FALSE.

        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: L_ilu, U_ilu

        CONTAINS

        PROCEDURE :: Allocate_ILU => Allocate_ILU_preconditioner
        PROCEDURE :: DEALLOCATE_ILU => DEALLOCATE_ILU_preconditioner
        PROCEDURE :: Calculate_ILU => Calculate_ILU

    END TYPE preconditioner

    TYPE(preconditioner) :: preconditioner_use

    CONTAINS

    SUBROUTINE Allocate_ILU_preconditioner(this, N)
        CLASS(preconditioner), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: N
        INTEGER :: alloc_stat

        this%ILU = .TRUE.
        ALLOCATE(this%L_ilu(N,N), this%U_ilu(N,N), stat=alloc_stat)
        IF (alloc_stat /= 0) STOP "ERROR :: Memory allocation failed for ILU preconditioner"
    END SUBROUTINE Allocate_ILU_preconditioner

    SUBROUTINE DEALLOCATE_ILU_preconditioner(this)
        CLASS(preconditioner), INTENT(INOUT) :: this

        IF (this%ILU) THEN
            DEALLOCATE(this%L_ilu, this%U_ilu)
            this%ILU = .FALSE.
        END IF
    END SUBROUTINE DEALLOCATE_ILU_preconditioner

    SUBROUTINE Calculate_ILU(this, A)
        CLASS(preconditioner), INTENT(INOUT) :: this
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A

        IF (this%ILU) CALL ILU_decomposition(A, this%L_ilu, this%U_ilu)

    END SUBROUTINE Calculate_ILU

    SUBROUTINE ILU_preconditioner(r, z, A)
        REAL(dp), DIMENSION(:), INTENT(IN) :: r
        REAL(dp), DIMENSION(:), INTENT(OUT) :: z
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION((SIZE(r))) :: y

        IF(.NOT. preconditioner_use%ILU) THEN
            CALL preconditioner_use%Allocate_ILU(SIZE(r))
            CALL preconditioner_use%Calculate_ILU(A)
        END IF

        y = forward(preconditioner_use%L_ilu, r)
        z = backward(preconditioner_use%U_ilu, y)

    END SUBROUTINE ILU_preconditioner

END MODULE NAFPack_Preconditioners
