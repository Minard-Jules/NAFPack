MODULE NAFPack_Iterative_types

    USE NAFPack_constant

    IMPLICIT NONE

    PRIVATE 

    PUBLIC :: MethodTypeIterative
    PUBLIC :: METHOD_ITERATIVE_NONE
    PUBLIC :: METHOD_Jacobi, METHOD_JOR
    PUBLIC :: METHOD_GAUSS_SEIDEL, METHOD_SOR, METHOD_SSOR
    PUBLIC :: METHOD_SIP_ILU, METHOD_SIP_ICF
    PUBLIC :: METHOD_RICHARDSON
    PUBLIC :: METHOD_CONJUGATE_GRADIENT

    TYPE :: MethodTypeIterative
        INTEGER :: value
        CHARACTER(LEN=64) :: name
    END TYPE MethodTypeIterative

    TYPE(MethodTypeIterative), PARAMETER :: METHOD_ITERATIVE_NONE = MethodTypeIterative(0, "None")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_Jacobi = MethodTypeIterative(1, "Jacobi")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_GAUSS_SEIDEL = MethodTypeIterative(2, "Gauss-Seidel")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_SOR = MethodTypeIterative(3, "Successive Over-Relaxation")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_JOR = MethodTypeIterative(4, "Jacobi Over-Relaxation")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_SIP_ILU = MethodTypeIterative(5, "Strongly Implicit Procedure with ILU")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_SIP_ICF = MethodTypeIterative(6, "Strongly Implicit Procedure with ICF")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_SSOR = MethodTypeIterative(7, "Symmetric Successive Over-Relaxation")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_RICHARDSON = MethodTypeIterative(8, "Richardson")
    TYPE(MethodTypeIterative), PARAMETER :: METHOD_CONJUGATE_GRADIENT = MethodTypeIterative(9, "Conjugate Gradient")

END MODULE NAFPack_Iterative_types