MODULE NAFPack_Direct_types

    IMPLICIT NONE

    PRIVATE 

    PUBLIC :: MethodTypeDirect, MethodQR
    PUBLIC :: METHOD_DIRECT_NONE
    PUBLIC :: METHOD_Gauss
    PUBLIC :: METHOD_LU, METHOD_LDU
    PUBLIC :: METHOD_CHOLESKY, METHOD_LDL_Cholesky, METHOD_QR
    PUBLIC :: METHOD_TDMA, METHOD_FADDEEV_LEVERRIER
    PUBLIC :: QR_HOUSEHOLDER, QR_GIVENS, QR_GRAM_SCHMIDT
    PUBLIC :: QR_GRAM_SCHMIDT_Modified

    TYPE :: MethodTypeDirect
        INTEGER :: value
        CHARACTER(LEN=64) :: name
    END TYPE MethodTypeDirect

    TYPE :: MethodQR
        INTEGER :: value
        CHARACTER(LEN=64) :: name
    END TYPE MethodQR

    TYPE(MethodTypeDirect), PARAMETER :: METHOD_DIRECT_NONE = MethodTypeDirect(0, "None")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_Gauss = MethodTypeDirect(1, "Gauss")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_LU = MethodTypeDirect(2, "LU")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_LDU = MethodTypeDirect(3, "LDU")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_CHOLESKY = MethodTypeDirect(4, "Cholesky")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_LDL_Cholesky = MethodTypeDirect(5, "LDL-Cholesky")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_QR = MethodTypeDirect(6, "QR")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_TDMA = MethodTypeDirect(7, "TDMA")
    TYPE(MethodTypeDirect), PARAMETER :: METHOD_FADDEEV_LEVERRIER = MethodTypeDirect(8, "Faddeev-Leverrier")

    TYPE(MethodQR), PARAMETER :: QR_HOUSEHOLDER = MethodQR(1, "Householder")
    TYPE(MethodQR), PARAMETER :: QR_GIVENS = MethodQR(2, "Givens")
    TYPE(MethodQR), PARAMETER :: QR_GRAM_SCHMIDT = MethodQR(3, "Gram-Schmidt")
    TYPE(MethodQR), PARAMETER :: QR_GRAM_SCHMIDT_Modified = MethodQR(4, "Gram_Schmidt_Modified")
    
END MODULE NAFPack_Direct_types