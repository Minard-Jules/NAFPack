MODULE NAFPack_Preconditioners

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    
    IMPLICIT NONE

    REAL(dp), ALLOCATABLE :: L_ilu(:,:), U_ilu(:,:)

CONTAINS

    SUBROUTINE Set_ILU(L, U)
        REAL(dp), INTENT(IN) :: L(:,:), U(:,:)
        
        L_ilu = L
        U_ilu = U

    END SUBROUTINE

    SUBROUTINE Apply_ILU(r, z)
        REAL(dp), DIMENSION(:), INTENT(IN) :: r
        REAL(dp), DIMENSION(:), INTENT(OUT) :: z
        REAL(dp), DIMENSION((SIZE(r))) :: y

        y = forward(L_ilu, r)

        z = backward(U_ilu, y)

    END SUBROUTINE

END MODULE NAFPack_Preconditioners
