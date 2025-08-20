module NAFPack_Iterative_types

    use NAFPack_constant

    implicit none(type, external)

    private

    public :: MethodTypeIterative
    public :: METHOD_ITERATIVE_NONE
    public :: METHOD_Jacobi, METHOD_JOR
    public :: METHOD_GAUSS_SEIDEL, METHOD_SOR, METHOD_SSOR
    public :: METHOD_SIP_ILU, METHOD_SIP_ICF
    public :: METHOD_RICHARDSON
    public :: METHOD_CONJUGATE_GRADIENT
    public :: METHOD_CONJUGATE_RESIDUAL
    public :: METHOD_CGNE, METHOD_CGNR
    public :: METHOD_GMRES

    public :: IterativeMethodRequirements

    public :: Norm_used
    public :: NORM_2, NORM_1, NORM_INF

    public :: relaxation_factor_used
    public :: RELAXATION_FACTOR_NONE, RELAXATION_FACTOR_OMEGA, RELAXATION_FACTOR_ALPHA

    integer, parameter :: CK = selected_char_kind('ISO_10646')
    character(KIND=ucs4, LEN=4), parameter :: none = "None"
    character(KIND=ucs4, LEN=1), parameter :: omega = char(int(z'03C9'), ucs4)
    character(KIND=ucs4, LEN=1), parameter :: alpha = char(int(z'03B1'), ucs4)

    type :: MethodTypeIterative
        integer :: id
        character(LEN=64) :: name
        character(LEN=64) :: name2 = ""
    end type MethodTypeIterative

    type :: IterativeMethodRequirements
        logical :: needs_SPD = .false.
        logical :: needs_diag_dom = .false.
        logical :: needs_square = .false.
        logical :: needs_symetric = .false.
    end type IterativeMethodRequirements

    type :: Norm_used
        integer :: id
        character(LEN=64) :: name
    end type Norm_used

    type :: relaxation_factor_used
        integer :: id
        character(KIND=ucs4, LEN=64) :: name
    end type relaxation_factor_used

    type(MethodTypeIterative), parameter :: METHOD_ITERATIVE_NONE = MethodTypeIterative(0, "None")
    type(MethodTypeIterative), parameter :: METHOD_Jacobi = MethodTypeIterative(1, "Jacobi")
    type(MethodTypeIterative), parameter :: METHOD_GAUSS_SEIDEL = MethodTypeIterative(2, "Gauss-Seidel")
    type(MethodTypeIterative), parameter :: METHOD_SOR = MethodTypeIterative(3, "Successive Over-Relaxation")
    type(MethodTypeIterative), parameter :: METHOD_JOR = MethodTypeIterative(4, "Jacobi Over-Relaxation")
    type(MethodTypeIterative), parameter :: METHOD_SIP_ILU = MethodTypeIterative(5, "Strongly Implicit Procedure", "ILU")
    type(MethodTypeIterative), parameter :: METHOD_SIP_ICF = MethodTypeIterative(6, "Strongly Implicit Procedure", "ICF")
    type(MethodTypeIterative), parameter :: METHOD_SSOR = MethodTypeIterative(7, "Symmetric Successive Over-Relaxation")
    type(MethodTypeIterative), parameter :: METHOD_RICHARDSON = MethodTypeIterative(8, "Richardson")
    type(MethodTypeIterative), parameter :: METHOD_CONJUGATE_GRADIENT = MethodTypeIterative(9, "Conjugate Gradient")
    type(MethodTypeIterative), parameter :: METHOD_CONJUGATE_RESIDUAL = MethodTypeIterative(10, "Conjugate Residual")
    type(MethodTypeIterative), parameter :: METHOD_CGNE = MethodTypeIterative(11, "Conjugate Gradient on Normal Equations")
    type(MethodTypeIterative), parameter :: METHOD_CGNR = MethodTypeIterative(12, "Conjugate Gradient on Normal Residual")
    type(MethodTypeIterative), parameter :: METHOD_GMRES = MethodTypeIterative(13, "Generalized Minimal Residual")

    type(Norm_used), parameter :: NORM_2 = Norm_used(1, "Norm L2 or Euclidean")
    type(Norm_used), parameter :: NORM_1 = Norm_used(2, "Norm L1 or Manhattan")
    type(Norm_used), parameter :: NORM_INF = Norm_used(3, "Norm LInfini or Maximum")

    type(relaxation_factor_used), parameter :: RELAXATION_FACTOR_NONE = relaxation_factor_used(0, none)
    type(relaxation_factor_used), parameter :: RELAXATION_FACTOR_OMEGA = relaxation_factor_used(1, omega)
    type(relaxation_factor_used), parameter :: RELAXATION_FACTOR_ALPHA = relaxation_factor_used(2, alpha)

end module NAFPack_Iterative_types
