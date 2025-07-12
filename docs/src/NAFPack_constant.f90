!> Module for defining constants used in NAFPack
!>
!> This module includes mathematical constants, colors for terminal output,
!> and other parameters that are used throughout the NAFPack library.
MODULE NAFPack_constant

    USE, INTRINSIC :: iso_fortran_env, ONLY: sp=>real32, dp=>real64, isp=>int32, idp=>int64
    
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: pi, im, epsi, kmax, sp, dp, epsi_test, idp, isp
    PUBLIC :: status_len
    PUBLIC :: red_color, green_color, yellow_color, blue_color, white_color, cyan_color, purple_color, reset_color
    PUBLIC :: NAF_SUCCESS, NAF_ERROR_DIMENSION, NAF_ERROR_SINGULAR, NAF_ERROR_CONVERGENCE, NAF_ERROR_MEMORY, &
              NAF_ERROR_INVALID_METHOD
    PUBLIC :: TOL_PIVOT, TOL_CONVERGENCE

    !> \( \pi \) constant
    REAL(dp), PARAMETER :: pi = ACOS(-1.d0)
    !> Imaginary unit \( i^2 = -1 \)
    COMPLEX(dp), PARAMETER :: im = (0.d0, 1.d0)
    !> Small \( \epsilon \) value
    REAL(dp), PARAMETER :: epsi = 1.d-12, epsi_test = 1.d-6
    !> Maximum number of iterations for iterative methods
    INTEGER, PARAMETER :: kmax = 10000
    !> len of status messages
    CHARACTER(LEN = *), PARAMETER :: status_len = REPEAT(" ", 15)

    !> red colors for terminal output
    CHARACTER(len=10) :: red_color    = CHAR(27)//"[31m"
    !> green colors for terminal output
    CHARACTER(len=10) :: green_color  = CHAR(27)//"[32m"
    !> yellow colors for terminal output
    CHARACTER(len=10) :: yellow_color = CHAR(27)//"[33m"
    !> blue colors for terminal output
    CHARACTER(len=10) :: blue_color  =  CHAR(27)//"[34m"
    !> purple colors for terminal output
    CHARACTER(len=10) :: purple_color = CHAR(27)//"[35m"
    !> cyan colors for terminal output
    CHARACTER(len=10) :: cyan_color   = CHAR(27)//"[36m"
    !> white colors for terminal output
    CHARACTER(len=10) :: white_color  = CHAR(27)//"[37m"
    !> reset colors for terminal output
    CHARACTER(len=10) :: reset_color  = CHAR(27)//"[0m"
    
    !> Error codes for better error handling
    INTEGER, PARAMETER :: NAF_SUCCESS = 0
    INTEGER, PARAMETER :: NAF_ERROR_DIMENSION = 1
    INTEGER, PARAMETER :: NAF_ERROR_SINGULAR = 2
    INTEGER, PARAMETER :: NAF_ERROR_CONVERGENCE = 3
    INTEGER, PARAMETER :: NAF_ERROR_MEMORY = 4
    INTEGER, PARAMETER :: NAF_ERROR_INVALID_METHOD = 5
    
    !> Improved tolerance parameters
    REAL(dp), PARAMETER :: TOL_PIVOT = 1.0e-14_dp
    REAL(dp), PARAMETER :: TOL_CONVERGENCE = 1.0e-12_dp

END MODULE NAFPack_constant