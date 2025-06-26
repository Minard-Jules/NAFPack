MODULE NAFPack_constant

    USE, INTRINSIC :: iso_fortran_env, ONLY: sp=>real32, dp=>real64, isp=>int32, idp=>int64
    
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: pi, im, epsi, kmax, sp, dp, epsi_test, idp, isp
    PUBLIC :: status_len
    PUBLIC :: red_color, green_color, yellow_color, blue_color, white_color, cyan_color, purple_color, reset_color

    REAL(dp), PARAMETER :: pi = ACOS(-1.d0)
    COMPLEX(dp), PARAMETER :: im = (0.d0, 1.d0)
    REAL(dp), PARAMETER :: epsi = 1.d-12, epsi_test = 1.d-6
    INTEGER, PARAMETER :: kmax = 10000
    CHARACTER(LEN = *), PARAMETER :: status_len = REPEAT(" ", 15)

    CHARACTER(len=10) :: red_color    = CHAR(27)//"[31m"
    CHARACTER(len=10) :: green_color  = CHAR(27)//"[32m"
    CHARACTER(len=10) :: yellow_color = CHAR(27)//"[33m"
    CHARACTER(len=10) :: blue_color  =  CHAR(27)//"[34m"
    CHARACTER(len=10) :: purple_color = CHAR(27)//"[35m"
    CHARACTER(len=10) :: cyan_color   = CHAR(27)//"[36m"
    CHARACTER(len=10) :: white_color  = CHAR(27)//"[37m"
    CHARACTER(len=10) :: reset_color  = CHAR(27)//"[0m"
    
END MODULE NAFPack_constant