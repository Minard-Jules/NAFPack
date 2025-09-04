!> Module for defining constants used in NAFPack
!>
!> This module includes mathematical constants,
!> and other parameters that are used throughout the NAFPack library.
module NAFPack_constant

    use NAFPack_kinds, only: sp, dp, qp

    implicit none(type, external)

    public

    !> \( \pi \) constant
    real(sp), parameter :: pi_sp = acos(-1.0_sp)
    real(dp), parameter :: pi_dp = acos(-1.0_dp)
    real(qp), parameter :: pi_qp = acos(-1.0_qp)
    !> Imaginary unit \( i^2 = -1 \)
    complex(sp), parameter :: im_sp = (0.0_sp, 1.0_sp)
    complex(dp), parameter :: im_dp = (0.0_dp, 1.0_dp)
    complex(qp), parameter :: im_qp = (0.0_qp, 1.0_qp)
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
