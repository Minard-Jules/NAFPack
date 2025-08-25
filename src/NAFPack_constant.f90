!> Module for defining constants used in NAFPack
!>
!> This module includes mathematical constants, colors for terminal output,
!> and other parameters that are used throughout the NAFPack library.
module NAFPack_constant

    use, intrinsic :: iso_fortran_env, only: output_unit, &
                                             sp => real32, &
                                             dp => real64, &
                                             isp => int32, &
                                             idp => int64

    implicit none(type, external)

    private
    public :: pi, im, epsi, kmax, sp, dp, epsi_test, idp, isp, int_inf, output_unit
    public :: status_len
    public :: red_color, green_color, &
              yellow_color, blue_color, &
              white_color, cyan_color, &
              purple_color, reset_color
    public :: red_color_ucs4, green_color_ucs4, yellow_color_ucs4, blue_color_ucs4, &
              white_color_ucs4, cyan_color_ucs4, purple_color_ucs4, reset_color_ucs4
    public :: ascii, ucs4
    public :: NAF_SUCCESS, NAF_ERROR_DIMENSION, &
              NAF_ERROR_SINGULAR, NAF_ERROR_CONVERGENCE, &
              NAF_ERROR_MEMORY, NAF_ERROR_INVALID_METHOD
    public :: TOL_PIVOT, TOL_CONVERGENCE

    !> \( \pi \) constant
    real(dp), parameter :: pi = acos(-1.d0)
    !> Imaginary unit \( i^2 = -1 \)
    complex(dp), parameter :: im = (0.d0, 1.d0)
    !> Small \( \epsilon \) value
    real(dp), parameter :: epsi = 1.d-12, epsi_test = 1.d-6
    !>
    integer, parameter :: int_inf = huge(1)
    !> Maximum number of iterations for iterative methods
    integer, parameter :: kmax = 10000
    !> len of status messages
    character(LEN=*), parameter :: status_len = repeat(" ", 15)

    integer, parameter :: ascii = selected_char_kind('ascii')
    integer, parameter :: ucs4 = selected_char_kind('ISO_10646')

    !> red colors for terminal output
    character(len=10) :: red_color = char(27)//"[31m"
    character(KIND=ucs4, LEN=10), parameter :: red_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[31m"
    !> green colors for terminal output
    character(len=10) :: green_color = char(27)//"[32m"
    character(KIND=ucs4, LEN=10), parameter :: green_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[32m"
    !> yellow colors for terminal output
    character(len=10) :: yellow_color = char(27)//"[33m"
    character(KIND=ucs4, LEN=10), parameter :: yellow_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[33m"
    !> blue colors for terminal output
    character(len=10) :: blue_color = char(27)//"[34m"
    character(KIND=ucs4, LEN=10), parameter :: blue_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[34m"
    !> purple colors for terminal output
    character(len=10) :: purple_color = char(27)//"[35m"
    character(KIND=ucs4, LEN=10), parameter :: purple_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[35m"
    !> cyan colors for terminal output
    character(len=10) :: cyan_color = char(27)//"[36m"
    character(KIND=ucs4, LEN=10), parameter :: cyan_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[36m"
    !> white colors for terminal output
    character(len=10) :: white_color = char(27)//"[37m"
    character(KIND=ucs4, LEN=10), parameter :: white_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[37m"
    !> reset colors for terminal output
    character(len=10) :: reset_color = char(27)//"[0m"
    character(KIND=ucs4, LEN=10), parameter :: reset_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[0m"

    !> Error codes for better error handling
    integer, parameter :: NAF_SUCCESS = 0
    integer, parameter :: NAF_ERROR_DIMENSION = 1
    integer, parameter :: NAF_ERROR_SINGULAR = 2
    integer, parameter :: NAF_ERROR_CONVERGENCE = 3
    integer, parameter :: NAF_ERROR_MEMORY = 4
    integer, parameter :: NAF_ERROR_INVALID_METHOD = 5

    !> Improved tolerance parameters
    real(dp), parameter :: TOL_PIVOT = 1.0e-14_dp
    real(dp), parameter :: TOL_CONVERGENCE = 1.0e-12_dp

end module NAFPack_constant
