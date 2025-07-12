MODULE NAFPack_linalg

    USE NAFPack_Direct_types
    USE NAFPack_Direct_method

    USE NAFPack_Iterative_types
    USE NAFPack_Iterative_methods

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: linalg, DirectMethod, IterativeMethod

    PUBLIC :: MethodTypeDirect, MethodQR
    PUBLIC :: METHOD_DIRECT_NONE
    PUBLIC :: METHOD_Gauss
    PUBLIC :: METHOD_LU, METHOD_LDU
    PUBLIC :: METHOD_CHOLESKY, METHOD_LDL_Cholesky, METHOD_QR
    PUBLIC :: METHOD_TDMA, METHOD_FADDEEV_LEVERRIER
    PUBLIC :: QR_HOUSEHOLDER, QR_GIVENS, QR_GRAM_SCHMIDT
    PUBLIC :: QR_GRAM_SCHMIDT_Modified

    PUBLIC :: IterativeParams
    PUBLIC :: METHOD_ITERATIVE_NONE
    PUBLIC :: METHOD_Jacobi, METHOD_JOR
    PUBLIC :: METHOD_GAUSS_SEIDEL, METHOD_SOR, METHOD_SSOR
    PUBLIC :: METHOD_SIP_ILU, METHOD_SIP_ICF
    PUBLIC :: METHOD_RICHARDSON

    PUBLIC :: MethodPreconditioner
    PUBLIC :: METHOD_PRECOND_NONE
    PUBLIC :: METHOD_PRECOND_JACOBI, METHOD_PRECOND_JOR
    PUBLIC :: METHOD_PRECOND_GS, METHOD_PRECOND_SOR, METHOD_PRECOND_SSOR
    PUBLIC :: METHOD_PRECOND_ILU, METHOD_PRECOND_ICF

    TYPE :: linalg

        TYPE(DirectMethod) :: direct
        TYPE(IterativeMethod) :: iterative

    END TYPE linalg

    CONTAINS

END MODULE NAFPack_linalg