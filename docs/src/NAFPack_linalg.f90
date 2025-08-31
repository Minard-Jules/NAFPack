module NAFPack_linalg

    use NAFPack_Direct_types, only: MethodTypeDirect, METHOD_DIRECT_NONE, &
                                    METHOD_CHOLESKY, METHOD_LDL_Cholesky, &
                                    METHOD_FADDEEV_LEVERRIER, &
                                    METHOD_Gauss, METHOD_Gauss_JORDAN, &
                                    METHOD_LU, METHOD_LDU, &
                                    METHOD_QR, METHOD_TDMA, &
                                    DirectMethodRequirements, MethodQR, &
                                    QR_HOUSEHOLDER, QR_GIVENS, &
                                    QR_GRAM_SCHMIDT, QR_GRAM_SCHMIDT_Modified
    use NAFPack_Direct_method, only: DirectMethod

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
    use NAFPack_Iterative_Params, only: IterativeParams
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
    use NAFPack_Iterative_methods, only: IterativeMethod

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
