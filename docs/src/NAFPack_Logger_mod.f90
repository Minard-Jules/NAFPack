!>
MODULE NAFPack_Logger_mod

    USE NAFPack_constant

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE

    PUBLIC :: Logger
    PUBLIC :: Format_file
    PUBLIC :: FORMAT_FILE_BIN, FORMAT_FILE_TXT, FORMAT_FILE_CSV, FORMAT_FILE_LOG, FORMAT_FILE_TSV
    PUBLIC :: FORMAT_FILE_JSON, FORMAT_FILE_XML, FORMAT_FILE_YAML

    PUBLIC :: log_field

    PUBLIC :: center_with_fill

    TYPE :: Format_file
        INTEGER :: id = 1
        CHARACTER(LEN=10) :: format_name = "txt"
        CHARACTER(LEN=100) :: format_description = "Text file format"
    END TYPE Format_file

    TYPE(Format_file), PARAMETER :: FORMAT_FILE_BIN = Format_file(0, "binary", "Binary file format")
    TYPE(Format_file), PARAMETER :: FORMAT_FILE_TXT = Format_file(1, "txt", "Text file format")
    TYPE(Format_file), PARAMETER :: FORMAT_FILE_CSV = Format_file(2, "csv", "Comma-separated values format")
    TYPE(Format_file), PARAMETER :: FORMAT_FILE_LOG = Format_file(3, "log", "Log file format")
    TYPE(Format_file), PARAMETER :: FORMAT_FILE_TSV = Format_file(4, "tsv", "Tab-separated values format")
    TYPE(Format_file), PARAMETER :: FORMAT_FILE_JSON = Format_file(5, "json", "JSON file format")
    TYPE(Format_file), PARAMETER :: FORMAT_FILE_XML = Format_file(6, "xml", "XML file format")
    TYPE(Format_file), PARAMETER :: FORMAT_FILE_YAML = Format_file(7, "yaml", "YAML file format")

    TYPE :: Logger
        INTEGER :: verbosity_level = 1
        LOGICAL :: to_terminal = .TRUE.
        LOGICAL :: to_file = .FALSE.
        INTEGER :: frequency = 10
        CHARACTER(LEN=100) :: filename = "Log"
        TYPE(Format_file) :: file_format = FORMAT_FILE_LOG
        INTEGER :: file_unit = 99
        CHARACTER(LEN=100) :: message = "Default log message"

        LOGICAL :: show_Logger_initialization = .TRUE.
        LOGICAL :: show_matrix_test = .TRUE.
        LOGICAL :: show_info_solver = .TRUE.
        LOGICAL :: show_iteration = .TRUE.
        LOGICAL :: show_final = .TRUE.

    CONTAINS

        PROCEDURE :: init => init_logger
        PROCEDURE :: log_info
        PROCEDURE :: log_detail
        PROCEDURE :: log_warning
        PROCEDURE :: log_error
        PROCEDURE :: log_time
        PROCEDURE :: WRITE => write_output
        PROCEDURE :: CLOSE => close_logger
    END TYPE Logger

    INTERFACE log_field
        MODULE PROCEDURE log_field_str, &
            log_field_real, &
            log_field_int, &
            log_field_ucs4, &
            log_field_logical
    END INTERFACE log_field

CONTAINS

!==========================================================================

    SUBROUTINE init_logger(this)
        CLASS(Logger), INTENT(INOUT) :: this
        CHARACTER(KIND=ucs4, LEN=100) :: msg

        IF (this%to_file) THEN
            OPEN (UNIT=this%file_unit, FILE=trim(this%filename)//"."//this%file_format%format_name, STATUS='REPLACE', &
                  ACTION='WRITE', ENCODING='UTF-8')
        END IF

        IF (this%to_terminal) THEN
            OPEN (output_unit, encoding='UTF-8')
        END IF

        IF (this%show_Logger_initialization) THEN
            CALL this%WRITE(center_with_fill("NAFPack Logger initialized", width=100, fill_char="="), box_style="top")
            CALL this%WRITE(ucs4_"", box_style="middle")

            CALL log_field(this, "Verbosity level", this%verbosity_level)
            CALL log_field(this, "Output to terminal", this%to_terminal)
            CALL log_field(this, "Output to file", this%to_file)
            IF (this%to_file) THEN
                CALL log_field(this, "File unit", this%file_unit)
                CALL log_field(this, "File format", this%file_format%format_name)
                CALL log_field(this, "File name", trim(this%filename)//"."//trim(this%file_format%format_name))
            END IF

            CALL this%WRITE(center_with_fill("", width=100, fill_char="="), box_style="bottom")
            CALL this%WRITE(ucs4_"")
        END IF

    END SUBROUTINE init_logger

!==========================================================================

    SUBROUTINE log_info(this, msg)
        CLASS(Logger), INTENT(INOUT) :: this
        CHARACTER(KIND=ucs4, LEN=*), INTENT(IN) :: msg

        IF (this%verbosity_level >= 2) CALL this%WRITE(msg, ucs4_"INFO", blue_color_ucs4)

    END SUBROUTINE log_info

    SUBROUTINE log_detail(this, msg)
        CLASS(Logger), INTENT(INOUT) :: this
        CHARACTER(KIND=ucs4, LEN=*), INTENT(IN) :: msg

        IF (this%verbosity_level >= 3) CALL this%WRITE(ucs4_"    "//msg, ucs4_"DETAIL", green_color_ucs4)

    END SUBROUTINE log_detail

    SUBROUTINE log_warning(this, msg)
        CLASS(Logger), INTENT(INOUT) :: this
        CHARACTER(KIND=ucs4, LEN=*), INTENT(IN) :: msg

        IF (this%verbosity_level >= 1) CALL this%WRITE(msg, ucs4_"WARNING", yellow_color_ucs4)

    END SUBROUTINE log_warning

    SUBROUTINE log_error(this, msg)
        CLASS(Logger), INTENT(INOUT) :: this
        CHARACTER(KIND=ucs4, LEN=*), INTENT(IN) :: msg

        IF (this%verbosity_level >= 1) CALL this%WRITE(msg, ucs4_"ERROR", red_color_ucs4)

    END SUBROUTINE log_error

    SUBROUTINE log_time(this, msg)
        CLASS(Logger), INTENT(INOUT) :: this
        CHARACTER(KIND=ucs4, LEN=*), INTENT(IN) :: msg
        CHARACTER(LEN=10) :: time
        CHARACTER(KIND=ucs4, LEN=10) :: time_ucs4

        CALL date_and_time(TIME=time)
        WRITE (time_ucs4, '(A)') time(:2)//":"//time(3:4)//":"//time(5:6)

        IF (this%verbosity_level >= 2) CALL this%WRITE(msg, time_ucs4, purple_color_ucs4)

    END SUBROUTINE log_time

!==========================================================================

    SUBROUTINE write_output(this, msg, name_level, color_level, box_style)
        CLASS(Logger), INTENT(IN) :: this
        CHARACTER(KIND=ucs4, LEN=*), INTENT(IN) :: msg
        CHARACTER(KIND=ucs4, LEN=*), OPTIONAL, INTENT(IN) :: name_level
        CHARACTER(KIND=ucs4, LEN=*), OPTIONAL, INTENT(IN) :: color_level
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: box_style
        CHARACTER(KIND=ucs4, LEN=100) :: info_char
        CHARACTER(LEN=4) :: box_char = " "

        IF (present(box_style)) THEN
            SELECT CASE (trim(adjustl(box_style)))
            CASE ("top")
                box_char = "╔"
            CASE ("bottom")
                box_char = "╚"
            CASE ("middle")
                box_char = "║ "
            CASE ("None")
                box_char = " "
            CASE DEFAULT
                box_char = " "
            END SELECT
        END IF

        info_char = ""
        IF (this%to_file) THEN
            IF (present(name_level)) THEN
                info_char = ucs4_"["//trim(name_level)//ucs4_"] "
                WRITE (this%file_unit, '(A, T15, "║ ", A)') trim(info_char), trim(msg)
            ELSE
                IF (present(box_style)) THEN
                    WRITE (this%file_unit, '(T15, A, A)') trim(box_char), trim(msg)
                ELSE
                    WRITE (this%file_unit, '(A)') trim(msg)
                END IF
            END IF
        END IF

        IF (this%to_terminal) THEN
            IF (present(name_level)) THEN
                IF (present(color_level)) THEN
                    info_char = ucs4_"["//trim(color_level)//trim(name_level)//trim(reset_color_ucs4)//ucs4_"] "
                ELSE
                    info_char = ucs4_"["//trim(name_level)//ucs4_"] "
                END IF
                WRITE (output_unit, '(A, T24, "║ ", A)') trim(info_char), trim(msg)
            ELSE
                IF (present(box_style)) THEN
                    WRITE (output_unit, '(T15, A, A)') trim(box_char), trim(msg)
                ELSE
                    WRITE (output_unit, '(A)') trim(msg)
                END IF
            END IF

        END IF

    END SUBROUTINE write_output

!==========================================================================

    SUBROUTINE close_logger(this)
        CLASS(Logger), INTENT(INOUT) :: this

        IF (this%to_file) CLOSE (this%file_unit)

    END SUBROUTINE close_logger

!==========================================================================

    SUBROUTINE log_field_str(verbose, label, VALUE)
        TYPE(Logger), INTENT(INOUT) :: verbose
        CHARACTER(*), INTENT(IN) :: label, VALUE
        CHARACTER(KIND=ucs4, LEN=100) :: msg

        WRITE (msg, '(A, T40, 2A)') trim(label), ": ", trim(VALUE)
        CALL verbose%log_info(msg)

    END SUBROUTINE log_field_str

    SUBROUTINE log_field_ucs4(verbose, label, VALUE)
        TYPE(Logger), INTENT(INOUT) :: verbose
        CHARACTER(*), INTENT(IN) :: label
        CHARACTER(KIND=ucs4, LEN=*), INTENT(IN) :: VALUE
        CHARACTER(KIND=ucs4, LEN=100) :: msg

        WRITE (msg, '(A, T40, 2A)') trim(label), ": ", trim(VALUE)
        CALL verbose%log_info(msg)

    END SUBROUTINE log_field_ucs4

    SUBROUTINE log_field_int(verbose, label, VALUE)
        TYPE(Logger), INTENT(INOUT) :: verbose
        CHARACTER(*), INTENT(IN) :: label
        INTEGER, INTENT(IN) :: VALUE
        CHARACTER(KIND=ucs4, LEN=100) :: msg

        WRITE (msg, '(A, T40, A, I0)') trim(label), ": ", VALUE
        CALL verbose%log_info(msg)

    END SUBROUTINE log_field_int

    SUBROUTINE log_field_real(verbose, label, VALUE)
        TYPE(Logger), INTENT(INOUT) :: verbose
        CHARACTER(*), INTENT(IN) :: label
        REAL(dp), INTENT(IN) :: VALUE
        CHARACTER(KIND=ucs4, LEN=100) :: msg

        WRITE (msg, '(A, T40, A, ES0.7)') trim(label), ": ", VALUE
        CALL verbose%log_info(msg)

    END SUBROUTINE log_field_real

    SUBROUTINE log_field_logical(verbose, label, VALUE)
        TYPE(Logger), INTENT(INOUT) :: verbose
        CHARACTER(*), INTENT(IN) :: label
        LOGICAL, INTENT(IN) :: VALUE
        CHARACTER(KIND=ucs4, LEN=100) :: msg

        WRITE (msg, '(A, T40, A, L)') trim(label), ": ", VALUE
        CALL verbose%log_info(msg)

    END SUBROUTINE log_field_logical

    FUNCTION center_with_fill(text, width, fill_char) RESULT(centered_text)
        CHARACTER(LEN=*), INTENT(IN) :: text
        INTEGER, INTENT(IN) :: width
        CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: fill_char
        CHARACTER(LEN=1) :: fill
        CHARACTER(KIND=ucs4, LEN=width) :: centered_text
        INTEGER :: text_len, padding, left_padding, right_padding, i

        IF (present(fill_char)) THEN
            fill = fill_char
        ELSE
            fill = " "
        END IF

        text_len = len_trim(text)
        IF (text_len >= width) THEN
            centered_text = text(1:width)
            RETURN
        END IF

        ! Calculate the total padding required
        padding = width - text_len
        IF (trim(text) == "") THEN
            left_padding = padding / 2
            right_padding = padding - left_padding - mod(padding, 2)
        ELSE
            text_len = text_len + 1
            left_padding = padding / 2 - 1
            right_padding = padding - left_padding - mod(padding, 2)
        END IF

        ! Initialize the result
        centered_text = repeat(' ', width)

        ! Fill with fill on the left
        DO i = 1, left_padding
            centered_text(i:i) = fill
        END DO

        ! Place the text in the center
        IF (trim(text) == "") THEN
            centered_text(left_padding + 1:left_padding + text_len) = trim(text)
        ELSE
            centered_text(left_padding + 1:left_padding + text_len) = " "//trim(text)//" "
        END IF

        ! Fill with fill on the right
        DO i = right_padding + text_len, width
            centered_text(i:i) = fill
        END DO

    END FUNCTION center_with_fill

END MODULE NAFPack_Logger_mod
