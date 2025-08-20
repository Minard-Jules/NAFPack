!>
module NAFPack_Logger_mod

    use NAFPack_constant

    implicit none(type, external)

    private

    public :: Logger
    public :: Format_file
    public :: FORMAT_FILE_BIN, FORMAT_FILE_TXT, FORMAT_FILE_CSV, FORMAT_FILE_LOG, FORMAT_FILE_TSV
    public :: FORMAT_FILE_JSON, FORMAT_FILE_XML, FORMAT_FILE_YAML

    public :: log_field

    public :: center_with_fill

    type :: Format_file
        integer :: id = 1
        character(LEN=10) :: format_name = "txt"
        character(LEN=100) :: format_description = "Text file format"
    end type Format_file

    type(Format_file), parameter :: FORMAT_FILE_BIN = Format_file(0, "binary", "Binary file format")
    type(Format_file), parameter :: FORMAT_FILE_TXT = Format_file(1, "txt", "Text file format")
    type(Format_file), parameter :: FORMAT_FILE_CSV = Format_file(2, "csv", "Comma-separated values format")
    type(Format_file), parameter :: FORMAT_FILE_LOG = Format_file(3, "log", "Log file format")
    type(Format_file), parameter :: FORMAT_FILE_TSV = Format_file(4, "tsv", "Tab-separated values format")
    type(Format_file), parameter :: FORMAT_FILE_JSON = Format_file(5, "json", "JSON file format")
    type(Format_file), parameter :: FORMAT_FILE_XML = Format_file(6, "xml", "XML file format")
    type(Format_file), parameter :: FORMAT_FILE_YAML = Format_file(7, "yaml", "YAML file format")

    type :: Logger
        integer :: verbosity_level = 1
        logical :: to_terminal = .true.
        logical :: to_file = .false.
        integer :: frequency = 10
        character(LEN=100) :: filename = "Log"
        type(Format_file) :: file_format = FORMAT_FILE_LOG
        integer :: file_unit = 99
        character(LEN=100) :: message = "Default log message"

        logical :: show_Logger_initialization = .true.
        logical :: show_matrix_test = .true.
        logical :: show_info_solver = .true.
        logical :: show_iteration = .true.
        logical :: show_final = .true.

    contains

        procedure :: init => init_logger
        procedure :: log_info
        procedure :: log_detail
        procedure :: log_warning
        procedure :: log_error
        procedure :: log_time
        procedure :: write => write_output
        procedure :: close => close_logger
    end type Logger

    interface log_field
        module procedure log_field_str, &
            log_field_real, &
            log_field_int, &
            log_field_ucs4, &
            log_field_logical
    end interface log_field

contains

!==========================================================================

    subroutine init_logger(this)
        class(Logger), intent(INOUT) :: this
        character(KIND=ucs4, LEN=100) :: msg

        if (this%to_file) then
            open (UNIT=this%file_unit, FILE=trim(this%filename)//"."//this%file_format%format_name, STATUS='REPLACE', &
                  ACTION='WRITE', ENCODING='UTF-8')
        end if

        if (this%to_terminal) then
            open (output_unit, encoding='UTF-8')
        end if

        if (this%show_Logger_initialization) then
            call this%write(center_with_fill("NAFPack Logger initialized", width=100, fill_char="="), box_style="top")
            call this%write(ucs4_"", box_style="middle")

            call log_field(this, "Verbosity level", this%verbosity_level)
            call log_field(this, "Output to terminal", this%to_terminal)
            call log_field(this, "Output to file", this%to_file)
            if (this%to_file) then
                call log_field(this, "File unit", this%file_unit)
                call log_field(this, "File format", this%file_format%format_name)
                call log_field(this, "File name", trim(this%filename)//"."//trim(this%file_format%format_name))
            end if

            call this%write(center_with_fill("", width=100, fill_char="="), box_style="bottom")
            call this%write(ucs4_"")
        end if

    end subroutine init_logger

!==========================================================================

    subroutine log_info(this, msg)
        class(Logger), intent(INOUT) :: this
        character(KIND=ucs4, LEN=*), intent(IN) :: msg

        if (this%verbosity_level >= 2) call this%write(msg, ucs4_"INFO", blue_color_ucs4)

    end subroutine log_info

    subroutine log_detail(this, msg)
        class(Logger), intent(INOUT) :: this
        character(KIND=ucs4, LEN=*), intent(IN) :: msg

        if (this%verbosity_level >= 3) call this%write(ucs4_"    "//msg, ucs4_"DETAIL", green_color_ucs4)

    end subroutine log_detail

    subroutine log_warning(this, msg)
        class(Logger), intent(INOUT) :: this
        character(KIND=ucs4, LEN=*), intent(IN) :: msg

        if (this%verbosity_level >= 1) call this%write(msg, ucs4_"WARNING", yellow_color_ucs4)

    end subroutine log_warning

    subroutine log_error(this, msg)
        class(Logger), intent(INOUT) :: this
        character(KIND=ucs4, LEN=*), intent(IN) :: msg

        if (this%verbosity_level >= 1) call this%write(msg, ucs4_"ERROR", red_color_ucs4)

    end subroutine log_error

    subroutine log_time(this, msg)
        class(Logger), intent(INOUT) :: this
        character(KIND=ucs4, LEN=*), intent(IN) :: msg
        character(LEN=10) :: time
        character(KIND=ucs4, LEN=10) :: time_ucs4

        call date_and_time(TIME=time)
        write (time_ucs4, '(A)') time(:2)//":"//time(3:4)//":"//time(5:6)

        if (this%verbosity_level >= 2) call this%write(msg, time_ucs4, purple_color_ucs4)

    end subroutine log_time

!==========================================================================

    subroutine write_output(this, msg, name_level, color_level, box_style)
        class(Logger), intent(IN) :: this
        character(KIND=ucs4, LEN=*), intent(IN) :: msg
        character(KIND=ucs4, LEN=*), optional, intent(IN) :: name_level
        character(KIND=ucs4, LEN=*), optional, intent(IN) :: color_level
        character(LEN=*), optional, intent(IN) :: box_style
        character(KIND=ucs4, LEN=100) :: info_char
        character(LEN=4) :: box_char = " "

        if (present(box_style)) then
            select case (trim(adjustl(box_style)))
            case ("top")
                box_char = "╔"
            case ("bottom")
                box_char = "╚"
            case ("middle")
                box_char = "║ "
            case ("None")
                box_char = " "
            case DEFAULT
                box_char = " "
            end select
        end if

        info_char = ""
        if (this%to_file) then
            if (present(name_level)) then
                info_char = ucs4_"["//trim(name_level)//ucs4_"] "
                write (this%file_unit, '(A, T15, "║ ", A)') trim(info_char), trim(msg)
            else
                if (present(box_style)) then
                    write (this%file_unit, '(T15, A, A)') trim(box_char), trim(msg)
                else
                    write (this%file_unit, '(A)') trim(msg)
                end if
            end if
        end if

        if (this%to_terminal) then
            if (present(name_level)) then
                if (present(color_level)) then
                    info_char = ucs4_"["//trim(color_level)//trim(name_level)//trim(reset_color_ucs4)//ucs4_"] "
                else
                    info_char = ucs4_"["//trim(name_level)//ucs4_"] "
                end if
                write (output_unit, '(A, T24, "║ ", A)') trim(info_char), trim(msg)
            else
                if (present(box_style)) then
                    write (output_unit, '(T15, A, A)') trim(box_char), trim(msg)
                else
                    write (output_unit, '(A)') trim(msg)
                end if
            end if

        end if

    end subroutine write_output

!==========================================================================

    subroutine close_logger(this)
        class(Logger), intent(INOUT) :: this

        if (this%to_file) close (this%file_unit)

    end subroutine close_logger

!==========================================================================

    subroutine log_field_str(verbose, label, value)
        type(Logger), intent(INOUT) :: verbose
        character(*), intent(IN) :: label, value
        character(KIND=ucs4, LEN=100) :: msg

        write (msg, '(A, T40, 2A)') trim(label), ": ", trim(value)
        call verbose%log_info(msg)

    end subroutine log_field_str

    subroutine log_field_ucs4(verbose, label, value)
        type(Logger), intent(INOUT) :: verbose
        character(*), intent(IN) :: label
        character(KIND=ucs4, LEN=*), intent(IN) :: value
        character(KIND=ucs4, LEN=100) :: msg

        write (msg, '(A, T40, 2A)') trim(label), ": ", trim(value)
        call verbose%log_info(msg)

    end subroutine log_field_ucs4

    subroutine log_field_int(verbose, label, value)
        type(Logger), intent(INOUT) :: verbose
        character(*), intent(IN) :: label
        integer, intent(IN) :: value
        character(KIND=ucs4, LEN=100) :: msg

        write (msg, '(A, T40, A, I0)') trim(label), ": ", value
        call verbose%log_info(msg)

    end subroutine log_field_int

    subroutine log_field_real(verbose, label, value)
        type(Logger), intent(INOUT) :: verbose
        character(*), intent(IN) :: label
        real(dp), intent(IN) :: value
        character(KIND=ucs4, LEN=100) :: msg

        write (msg, '(A, T40, A, ES0.7)') trim(label), ": ", value
        call verbose%log_info(msg)

    end subroutine log_field_real

    subroutine log_field_logical(verbose, label, value)
        type(Logger), intent(INOUT) :: verbose
        character(*), intent(IN) :: label
        logical, intent(IN) :: value
        character(KIND=ucs4, LEN=100) :: msg

        write (msg, '(A, T40, A, L)') trim(label), ": ", value
        call verbose%log_info(msg)

    end subroutine log_field_logical

    function center_with_fill(text, width, fill_char) result(centered_text)
        character(LEN=*), intent(IN) :: text
        integer, intent(IN) :: width
        character(LEN=1), optional, intent(IN) :: fill_char
        character(LEN=1) :: fill
        character(KIND=ucs4, LEN=width) :: centered_text
        integer :: text_len, padding, left_padding, right_padding, i

        if (present(fill_char)) then
            fill = fill_char
        else
            fill = " "
        end if

        text_len = len_trim(text)
        if (text_len >= width) then
            centered_text = text(1:width)
            return
        end if

        ! Calculate the total padding required
        padding = width - text_len
        if (trim(text) == "") then
            left_padding = padding / 2
            right_padding = padding - left_padding - mod(padding, 2)
        else
            text_len = text_len + 1
            left_padding = padding / 2 - 1
            right_padding = padding - left_padding - mod(padding, 2)
        end if

        ! Initialize the result
        centered_text = repeat(' ', width)

        ! Fill with fill on the left
        do i = 1, left_padding
            centered_text(i:i) = fill
        end do

        ! Place the text in the center
        if (trim(text) == "") then
            centered_text(left_padding + 1:left_padding + text_len) = trim(text)
        else
            centered_text(left_padding + 1:left_padding + text_len) = " "//trim(text)//" "
        end if

        ! Fill with fill on the right
        do i = right_padding + text_len, width
            centered_text(i:i) = fill
        end do

    end function center_with_fill

end module NAFPack_Logger_mod
