module NAFPack_linalg

    use NAFPack_Direct_types
    use NAFPack_Direct_method

    use NAFPack_Iterative_types
    use NAFPack_Iterative_methods

    implicit none(type, external)

    private

    public :: linalg, DirectMethod, IterativeMethod

    public :: MethodTypeDirect, MethodQR
    public :: METHOD_DIRECT_NONE
    public :: METHOD_Gauss, METHOD_Gauss_JORDAN
    public :: METHOD_LU, METHOD_LDU
    public :: METHOD_CHOLESKY, METHOD_LDL_Cholesky, METHOD_QR
    public :: METHOD_TDMA, METHOD_FADDEEV_LEVERRIER
    public :: QR_HOUSEHOLDER, QR_GIVENS, QR_GRAM_SCHMIDT
    public :: QR_GRAM_SCHMIDT_Modified

    public :: IterativeParams
    public :: METHOD_ITERATIVE_NONE
    public :: METHOD_Jacobi, METHOD_JOR
    public :: METHOD_GAUSS_SEIDEL, METHOD_SOR, METHOD_SSOR
    public :: METHOD_SIP_ILU, METHOD_SIP_ICF
    public :: METHOD_RICHARDSON
    public :: METHOD_CONJUGATE_GRADIENT
    public :: METHOD_CONJUGATE_RESIDUAL
    public :: METHOD_CGNE, METHOD_CGNR
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

    type :: linalg

        type(DirectMethod) :: direct
        type(IterativeMethod) :: iterative

    end type linalg

contains

end module NAFPack_linalg
