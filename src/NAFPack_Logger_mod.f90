!> 
MODULE NAFPack_Logger_mod

    USE NAFPack_constant

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: Logger
    
    TYPE :: Logger
        LOGICAL :: to_terminal = .TRUE.
        LOGICAL :: to_file = .FALSE.
        INTEGER :: frequency = 10
        CHARACTER(LEN=100) :: filename = "log.txt"
        INTEGER :: file_unit = 99
        INTEGER :: step = 0
        LOGICAL :: show_iteration = .TRUE.
        LOGICAL :: show_final = .TRUE.

        CONTAINS
        
        PROCEDURE :: log => log_message
        PROCEDURE :: init => init_logger
        PROCEDURE :: close => close_logger
    END TYPE Logger

    CONTAINS

    SUBROUTINE init_logger(this)
        CLASS(Logger), INTENT(INOUT) :: this
        IF (this%to_file) THEN
        OPEN(UNIT=this%file_unit, FILE=TRIM(this%filename), STATUS='REPLACE', ACTION='WRITE')
        END IF
    END SUBROUTINE init_logger

    SUBROUTINE log_message(this, msg)
        CLASS(Logger), INTENT(INOUT) :: this
        CHARACTER(LEN=*), INTENT(IN) :: msg

        IF (MOD(this%step, this%frequency) /= 0 .AND. this%show_iteration) RETURN

        IF (this%to_terminal) PRINT *, msg
        IF (this%to_file) WRITE(this%file_unit, '(A)') TRIM(msg)
    END SUBROUTINE log_message

    SUBROUTINE close_logger(this)
        CLASS(Logger), INTENT(INOUT) :: this
        IF (this%to_file) CLOSE(this%file_unit)
    END SUBROUTINE close_logger

END MODULE NAFPack_Logger_mod