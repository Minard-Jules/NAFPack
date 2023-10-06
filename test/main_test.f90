PROGRAM test

    USE NAFPack_constant
    USE test_NAFPack_linear_algebra
    USE test_NAFPack_fft

    IMPLICIT NONE
    
    LOGICAL :: stat = .FALSE.

    WRITE(*,'(A)') purple_color,"Test linear systeme :", reset_color
    PRINT*, " "
    PRINT*, " "
    CALL test_linear_system(stat)
    PRINT*, " "

    WRITE(*,'(A)') purple_color,"Test linear algebre :", reset_color
    PRINT*, " "
    PRINT*, " "
    CALL test_linear_algebra(stat)
    PRINT*, " "

    WRITE(*,'(A)') purple_color,"Test FFT :", reset_color
    PRINT*, " "
    PRINT*, " "
    CALL test_FFT(stat)
    PRINT*, " "

    IF(stat)THEN
        WRITE(*,'(A)') red_color,"Test failed", reset_color
    ELSE
        WRITE(*,'(A)') green_color,"Test success", reset_color
    END IF

END PROGRAM test