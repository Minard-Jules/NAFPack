submodule(NAFPack_ANSI) NAFPack_ANSI_ASCII

    character(len=1, kind=ascii), parameter :: esc = char(27)
    character(len=2, kind=ascii), parameter :: CSI = esc//"["
    character(len=1, kind=ascii), parameter :: final_character = "m"
    character(len=1, kind=ascii), parameter :: delimiter = ";"

contains

    module subroutine init_colors_ascii(this)
        class(ColorsAscii), intent(out) :: this
        this%reset = create_ansi_ascii(Ansi_Constants%STYLE_RESET)
        this%bold = create_ansi_ascii(Ansi_Constants%STYLE_BOLD)
        this%faint = create_ansi_ascii(Ansi_Constants%STYLE_FAINT)
        this%italic = create_ansi_ascii(Ansi_Constants%STYLE_ITALIC)
        this%underline = create_ansi_ascii(Ansi_Constants%STYLE_UNDERLINE)
        this%blink = create_ansi_ascii(Ansi_Constants%STYLE_BLINK)
        this%reverse = create_ansi_ascii(Ansi_Constants%STYLE_REVERSE)
        this%hidden = create_ansi_ascii(Ansi_Constants%STYLE_HIDDEN)
        this%strikethrough = create_ansi_ascii(Ansi_Constants%STYLE_STRIKETHROUGH)
        this%red = create_ansi_ascii(Ansi_Constants%FG_RED)
        this%green = create_ansi_ascii(Ansi_Constants%FG_GREEN)
        this%yellow = create_ansi_ascii(Ansi_Constants%FG_YELLOW)
        this%blue = create_ansi_ascii(Ansi_Constants%FG_BLUE)
        this%magenta = create_ansi_ascii(Ansi_Constants%FG_MAGENTA)
        this%cyan = create_ansi_ascii(Ansi_Constants%FG_CYAN)
        this%white = create_ansi_ascii(Ansi_Constants%FG_WHITE)
        this%bright_red = create_ansi_ascii(Ansi_Constants%FG_BRIGHT_RED)
        this%bright_green = create_ansi_ascii(Ansi_Constants%FG_BRIGHT_GREEN)
        this%bright_yellow = create_ansi_ascii(Ansi_Constants%FG_BRIGHT_YELLOW)
        this%bright_blue = create_ansi_ascii(Ansi_Constants%FG_BRIGHT_BLUE)
        this%bright_magenta = create_ansi_ascii(Ansi_Constants%FG_BRIGHT_MAGENTA)
        this%bright_cyan = create_ansi_ascii(Ansi_Constants%FG_BRIGHT_CYAN)
        this%bright_white = create_ansi_ascii(Ansi_Constants%FG_BRIGHT_WHITE)
    end subroutine init_colors_ascii

    pure module function create_ansi_ascii(ansi_code) result(ansi_string)
        type(AnsiCode), intent(in) :: ansi_code
        character(len=:, kind=ascii), allocatable :: ansi_string
        integer(i8) :: style
        integer(i16) :: fg, bg

        ansi_string = CSI//"0" !Always reset

        if (ansi_code%use_style) then
            style = ansi_code%style
            ansi_string = ansi_string//delimiter//to_str_ascii(style)
        end if

        if (ansi_code%use_fg) then
            fg = int(ansi_code%fg, kind=i16)
            ansi_string = ansi_string//delimiter//"38;5;"//to_str_ascii(fg)
        end if

        if (ansi_code%use_bg) then
            bg = int(ansi_code%bg, kind=i16)
            ansi_string = ansi_string//delimiter//"48;5;"//to_str_ascii(bg)
        end if

        if (.not. (ansi_code%use_style .or. ansi_code%use_fg .or. ansi_code%use_bg)) then
            style = ansi_code%style
            ansi_string = ansi_string//delimiter//to_str_ascii(style)
        end if

        ansi_string = ansi_string//"m"

    end function create_ansi_ascii

    pure module function apply_style_ascii(text, style) result(styled_text)
        character(*, kind=ascii), intent(in) :: text
        type(AnsiCode), intent(in) :: style
        character(:, kind=ascii), allocatable :: styled_text

        styled_text = create_ansi_ascii(style)//text//create_ansi_ascii(Ansi_Constants%STYLE_RESET)
    end function apply_style_ascii

    pure module function colorize_text_ascii(text, color) result(colored_text)
        character(*, kind=ascii), intent(in) :: text
        integer(i16), intent(in) :: color
        character(:, kind=ascii), allocatable :: colored_text

        colored_text = apply_style_ascii(text, set_ansi_code(fg_color=color))
    end function colorize_text_ascii

    subroutine cursor_position_ascii(row, col)
        integer, intent(in) :: row, col
        character(:, kind=ascii), allocatable :: sequence
        character(:, kind=ascii), allocatable :: temp

        write (temp, '(I0,A,I0)') row, delimiter, col
        sequence = CSI//temp//'H'
        write (output_unit, '(A)', advance='no') sequence
    end subroutine cursor_position_ascii

    subroutine clear_screen_ascii()
        write (*, '(A)', advance='no') CSI//'2J'
    end subroutine clear_screen_ascii

    subroutine clear_line_ascii()
        write (*, '(A)', advance='no') CSI//'2K'
    end subroutine clear_line_ascii

    subroutine save_cursor_ascii()
        write (*, '(A)', advance='no') CSI//'s'
    end subroutine save_cursor_ascii

    subroutine restore_cursor_ascii()
        write (*, '(A)', advance='no') CSI//'u'
    end subroutine restore_cursor_ascii

end submodule NAFPack_ANSI_ASCII
