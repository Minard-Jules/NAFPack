MODULE NAFPack_Iterative_types

    USE NAFPack_constant

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE

    PUBLIC :: MethodTypeIterative
    PUBLIC :: METHOD_ITERATIVE_NONE
    PUBLIC :: METHOD_Jacobi, METHOD_JOR
    PUBLIC :: METHOD_GAUSS_SEIDEL, METHOD_SOR, METHOD_SSOR
    PUBLIC :: METHOD_SIP_ILU, METHOD_SIP_ICF
    PUBLIC :: METHOD_RICHARDSON
    PUBLIC :: METHOD_CONJUGATE_GRADIENT
    PUBLIC :: METHOD_CONJUGATE_RESIDUAL
    PUBLIC :: METHOD_CGNE, METHOD_CGNR
    PUBLIC :: METHOD_GMRES

    PUBLIC :: IterativeMethodRequirements

    PUBLIC :: Norm_used
    PUBLIC :: NORM_2, NORM_1, NORM_INF

    PUBLIC :: relaxation_factor_used
    PUBLIC :: RELAXATION_FACTOR_NONE, RELAXATION_FACTOR_OMEGA, RELAXATION_FACTOR_ALPHA

    INTEGER, PARAMETER :: CK = selected_char_kind('ISO_10646')
    CHARACTER(KIND=ucs4, LEN=4), PARAMETER :: NONE = "None"
    CHARACTER(KIND=ucs4, LEN=1), PARAMETER :: omega = char(int(z'03C9'), ucs4)
    CHARACTER(KIND=ucs4, LEN=1), PARAMETER :: alpha = char(int(z'03B1'), ucs4)

    TYPE :: MethodTypeIterative
        INTEGER :: id
        CHARACTER(LEN=64) :: name
        CHARACTER(LEN=64) :: name2 = ""
    END TYPE MethodTypeIterative

    TYPE :: IterativeMethodRequirements
        LOGICAL :: needs_SPD = .FALSE.
        LOGICAL :: needs_diag_dom = .FALSE.
        LOGICAL :: needs_square = .FALSE.
        LOGICAL :: needs_symetric = .FALSE.
    END TYPE IterativeMethodRequirements

    TYPE :: Norm_used
        INTEGER :: id
        CHARACTER(LEN=64) :: name
    END TYPE Norm_used

    TYPE :: relaxation_factor_used
        INTEGER :: id
        CHARACTER(KIND=ucs4, LEN=64) :: name
    END TYPE relaxation_factor_used

    TYPE(MethodTypeIterative), PARAMETER :: METHOD_ITERATIVE_NONE = MethodTypeIterative(0, "None")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_Jacobi = MethodTypeIterative(1, "Jacobi")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_GAUSS_SEIDEL = MethodTypeIterative(2, "Gauss-Seidel")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_SOR = MethodTypeIterative(3, "Successive Over-Relaxation")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_JOR = MethodTypeIterative(4, "Jacobi Over-Relaxation")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_SIP_ILU = MethodTypeIterative(5, "Strongly Implicit Procedure", "ILU")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_SIP_ICF = MethodTypeIterative(6, "Strongly Implicit Procedure", "ICF")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_SSOR = MethodTypeIterative(7, "Symmetric Successive Over-Relaxation")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_RICHARDSON = MethodTypeIterative(8, "Richardson")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_CONJUGATE_GRADIENT = MethodTypeIterative(9, "Conjugate Gradient")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_CONJUGATE_RESIDUAL = MethodTypeIterative(10, "Conjugate Residual")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_CGNE = MethodTypeIterative(11, "Conjugate Gradient on Normal Equations")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_CGNR = MethodTypeIterative(12, "Conjugate Gradient on Normal Residual")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_GMRES = MethodTypeIterative(13, "Generalized Minimal Residual")

    TYPE(Norm_used), PARAMETER :: NORM_2 = Norm_used(1, "Norm L2 or Euclidean")
    TYPE(Norm_used), PARAMETER :: NORM_1 = Norm_used(2, "Norm L1 or Manhattan")
    TYPE(Norm_used), PARAMETER :: NORM_INF = Norm_used(3, "Norm LInfini or Maximum")

    TYPE(relaxation_factor_used), PARAMETER :: RELAXATION_FACTOR_NONE = relaxation_factor_used(0, NONE)
    TYPE(relaxation_factor_used), PARAMETER :: RELAXATION_FACTOR_OMEGA = relaxation_factor_used(1, omega)
    TYPE(relaxation_factor_used), PARAMETER :: RELAXATION_FACTOR_ALPHA = relaxation_factor_used(2, alpha)

END MODULE NAFPack_Iterative_types
