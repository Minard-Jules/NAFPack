module NAFPack_Direct_types

    implicit none(type, external)

    private

    public :: MethodTypeDirect, MethodQR
    public :: METHOD_DIRECT_NONE
    public :: METHOD_Gauss, METHOD_Gauss_JORDAN
    public :: METHOD_LU, METHOD_LDU
    public :: METHOD_CHOLESKY, METHOD_LDL_Cholesky, METHOD_QR
    public :: METHOD_TDMA, METHOD_FADDEEV_LEVERRIER
    public :: QR_HOUSEHOLDER, QR_GIVENS, QR_GRAM_SCHMIDT
    public :: QR_GRAM_SCHMIDT_Modified

    public :: DirectMethodRequirements

    type :: MethodTypeDirect
        integer :: id
        character(LEN=64) :: name
    end type MethodTypeDirect

    type :: MethodQR
        integer :: id
        character(LEN=64) :: name
    end type MethodQR

    type :: DirectMethodRequirements
        logical :: needs_SPD = .false.
        logical :: needs_non_zero_diag = .false.
        logical :: needs_square = .false.
        logical :: needs_tridiagonal = .false.
        logical :: needs_symmetric = .false.
    end type DirectMethodRequirements

    type(MethodTypeDirect), parameter :: METHOD_DIRECT_NONE = MethodTypeDirect(0, "None")
    type(MethodTypeDirect), parameter :: METHOD_Gauss = MethodTypeDirect(1, "Gauss")
    type(MethodTypeDirect), parameter :: METHOD_Gauss_JORDAN = MethodTypeDirect(2, "Gauss-Jordan")
    type(MethodTypeDirect), parameter :: METHOD_LU = MethodTypeDirect(3, "LU")
    type(MethodTypeDirect), parameter :: METHOD_LDU = MethodTypeDirect(4, "LDU")
    type(MethodTypeDirect), parameter :: METHOD_CHOLESKY = MethodTypeDirect(5, "Cholesky")
    type(MethodTypeDirect), parameter :: METHOD_LDL_Cholesky = MethodTypeDirect(6, "LDL-Cholesky")
    type(MethodTypeDirect), parameter :: METHOD_QR = MethodTypeDirect(7, "QR")
    type(MethodTypeDirect), parameter :: METHOD_TDMA = MethodTypeDirect(8, "TDMA")
    type(MethodTypeDirect), parameter :: METHOD_FADDEEV_LEVERRIER = MethodTypeDirect(9, "Faddeev-Leverrier")

    type(MethodQR), parameter :: QR_HOUSEHOLDER = MethodQR(1, "Householder")
    type(MethodQR), parameter :: QR_GIVENS = MethodQR(2, "Givens")
    type(MethodQR), parameter :: QR_GRAM_SCHMIDT = MethodQR(3, "Gram-Schmidt")
    type(MethodQR), parameter :: QR_GRAM_SCHMIDT_Modified = MethodQR(4, "Gram_Schmidt_Modified")

end module NAFPack_Direct_types
