MODULE NAFPack_Direct_types

    IMPLICIT NONE

    PRIVATE 

    PUBLIC :: MethodTypeDirect, MethodQR
    PUBLIC :: METHOD_DIRECT_NONE
    PUBLIC :: METHOD_Gauss, METHOD_Gauss_JORDAN
    PUBLIC :: METHOD_LU, METHOD_LDU
    PUBLIC :: METHOD_CHOLESKY, METHOD_LDL_Cholesky, METHOD_QR
    PUBLIC :: METHOD_TDMA, METHOD_FADDEEV_LEVERRIER
    PUBLIC :: QR_HOUSEHOLDER, QR_GIVENS, QR_GRAM_SCHMIDT
    PUBLIC :: QR_GRAM_SCHMIDT_Modified

    PUBLIC :: DirectMethodRequirements

    TYPE :: MethodTypeDirect
        INTEGER :: value
        CHARACTER(LEN=64) :: name
    END TYPE MethodTypeDirect

    TYPE :: MethodQR
        INTEGER :: value
        CHARACTER(LEN=64) :: name
    END TYPE MethodQR

    TYPE :: DirectMethodRequirements
        LOGICAL :: needs_SPD = .FALSE.
        LOGICAL :: needs_non_zero_diag = .FALSE.
        LOGICAL :: needs_square = .FALSE.
        LOGICAL :: needs_tridiagonal = .FALSE.
        LOGICAL :: needs_symmetric = .FALSE.
    END TYPE DirectMethodRequirements

    TYPE(MethodTypeDirect), PARAMETER :: METHOD_DIRECT_NONE = MethodTypeDirect(0, "None")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_Gauss = MethodTypeDirect(1, "Gauss")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_Gauss_JORDAN = MethodTypeDirect(2, "Gauss-Jordan")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_LU = MethodTypeDirect(3, "LU")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_LDU = MethodTypeDirect(4, "LDU")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_CHOLESKY = MethodTypeDirect(5, "Cholesky")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_LDL_Cholesky = MethodTypeDirect(6, "LDL-Cholesky")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_QR = MethodTypeDirect(7, "QR")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_TDMA = MethodTypeDirect(8, "TDMA")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_FADDEEV_LEVERRIER = MethodTypeDirect(9, "Faddeev-Leverrier")

    TYPE(MethodQR), PARAMETER :: QR_HOUSEHOLDER = MethodQR(1, "Householder")
    TYPE(MethodQR), PARAMETER :: QR_GIVENS = MethodQR(2, "Givens")
    TYPE(MethodQR), PARAMETER :: QR_GRAM_SCHMIDT = MethodQR(3, "Gram-Schmidt")
    TYPE(MethodQR), PARAMETER :: QR_GRAM_SCHMIDT_Modified = MethodQR(4, "Gram_Schmidt_Modified")
    
END MODULE NAFPack_Direct_types