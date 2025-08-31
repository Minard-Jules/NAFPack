!> Module for defining constants used in NAFPack
!>
!> This module includes mathematical constants,
!> and other parameters that are used throughout the NAFPack library.
module NAFPack_constant

    use NAFPack_kinds, only: dp

    implicit none(type, external)

    public

    !> \( \pi \) constant
    real(dp), parameter :: pi = acos(-1.0_dp)
    !> Imaginary unit \( i^2 = -1 \)
    complex(dp), parameter :: im = (0.0_dp, 1.0_dp)
    !> Integer infinity
    integer, parameter :: int_inf = huge(1)

    !> Error codes for better error handling
    integer, parameter :: NAF_SUCCESS = 0
    integer, parameter :: NAF_ERROR_DIMENSION = 1
    integer, parameter :: NAF_ERROR_SINGULAR = 2
    integer, parameter :: NAF_ERROR_CONVERGENCE = 3
    integer, parameter :: NAF_ERROR_MEMORY = 4
    integer, parameter :: NAF_ERROR_INVALID_METHOD = 5

    ! Numerical tolerances
    real(dp), parameter :: TOL_PIVOT = 1.0e-14_dp
    real(dp), parameter :: TOL_CONVERGENCE = 1.0e-12_dp
    real(dp), parameter :: TOL_RESIDUAL = 1.0e-10_dp
    real(dp), parameter :: TOL_TEST = 1.0e-10_dp


    ! Performance settings
    integer, parameter :: MAX_ITERATION = 10000
    ! integer, parameter :: block_size = 64
    ! logical, parameter :: use_openmp = .true.
    ! logical, parameter :: use_blas = .true.
    ! logical, parameter :: use_lapack = .true.

end module NAFPack_constant
