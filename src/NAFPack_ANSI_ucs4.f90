submodule(NAFPack_ANSI) NAFPack_ANSI_ucs4

    character(len=1, kind=ucs4), parameter :: esc = char(27)
    character(len=2, kind=ucs4), parameter :: CSI = esc//ucs4_"["
    character(len=1, kind=ucs4), parameter :: final_character = ucs4_"m"
    character(len=1, kind=ucs4), parameter :: delimiter = ucs4_";"

contains

    module subroutine init_colors_ucs4(this)
        class(ColorsUcs4), intent(out) :: this
        this%reset         = create_ansi_ucs4(Ansi_Constants%STYLE_RESET)
        this%bold          = create_ansi_ucs4(Ansi_Constants%STYLE_BOLD)
        this%faint         = create_ansi_ucs4(Ansi_Constants%STYLE_FAINT)
        this%italic        = create_ansi_ucs4(Ansi_Constants%STYLE_ITALIC)
        this%underline     = create_ansi_ucs4(Ansi_Constants%STYLE_UNDERLINE)
        this%blink         = create_ansi_ucs4(Ansi_Constants%STYLE_BLINK)
        this%reverse       = create_ansi_ucs4(Ansi_Constants%STYLE_REVERSE)
        this%hidden        = create_ansi_ucs4(Ansi_Constants%STYLE_HIDDEN)
        this%strikethrough = create_ansi_ucs4(Ansi_Constants%STYLE_STRIKETHROUGH)
        this%red           = create_ansi_ucs4(Ansi_Constants%FG_RED)
        this%green         = create_ansi_ucs4(Ansi_Constants%FG_GREEN)
        this%yellow        = create_ansi_ucs4(Ansi_Constants%FG_YELLOW)
        this%blue          = create_ansi_ucs4(Ansi_Constants%FG_BLUE)
        this%magenta       = create_ansi_ucs4(Ansi_Constants%FG_MAGENTA)
        this%cyan          = create_ansi_ucs4(Ansi_Constants%FG_CYAN)
        this%white         = create_ansi_ucs4(Ansi_Constants%FG_WHITE)
        this%bright_red    = create_ansi_ucs4(Ansi_Constants%FG_BRIGHT_RED)
        this%bright_green  = create_ansi_ucs4(Ansi_Constants%FG_BRIGHT_GREEN)
        this%bright_yellow = create_ansi_ucs4(Ansi_Constants%FG_BRIGHT_YELLOW)
        this%bright_blue   = create_ansi_ucs4(Ansi_Constants%FG_BRIGHT_BLUE)
        this%bright_magenta= create_ansi_ucs4(Ansi_Constants%FG_BRIGHT_MAGENTA)
        this%bright_cyan   = create_ansi_ucs4(Ansi_Constants%FG_BRIGHT_CYAN)
        this%bright_white  = create_ansi_ucs4(Ansi_Constants%FG_BRIGHT_WHITE)
    end subroutine init_colors_ucs4

    pure module function create_ansi_ucs4(ansi_code) result(ansi_string)
        type(AnsiCode), intent(in) :: ansi_code
        character(len=:, kind=ucs4), allocatable :: ansi_string
        integer(i8) :: style
        integer(i16) :: fg, bg

        ansi_string = CSI//ucs4_"0" !Always reset

        if (ansi_code%use_style) then
            style = ansi_code%style
            ansi_string = ansi_string//delimiter//to_str_ucs4(style)
        end if

        if (ansi_code%use_fg) then
            fg = int(ansi_code%fg, kind=i16)
            ansi_string = ansi_string//delimiter//ucs4_"38;5;"//to_str_ucs4(fg)
        end if

        if (ansi_code%use_bg) then
            bg = int(ansi_code%bg, kind=i16)
            ansi_string = ansi_string//delimiter//ucs4_"48;5;"//to_str_ucs4(bg)
        end if

        if (.not. (ansi_code%use_style .or. ansi_code%use_fg .or. ansi_code%use_bg)) then
            style = ansi_code%style
            ansi_string = ansi_string//delimiter//to_str_ucs4(style)
        end if

        ansi_string = ansi_string//ucs4_"m"

    end function create_ansi_ucs4

    pure module function apply_style_ucs4(text, style) result(styled_text)
        character(*, kind=ucs4), intent(in) :: text
        type(AnsiCode), intent(in) :: style
        character(:, kind=ucs4), allocatable :: styled_text

        styled_text = create_ansi_ucs4(style)//text//create_ansi_ucs4(Ansi_Constants%STYLE_RESET)
    end function apply_style_ucs4

    pure module function colorize_text_ucs4(text, color) result(colored_text)
        character(*, kind=ucs4), intent(in) :: text
        integer(i16), intent(in) :: color
        character(:, kind=ucs4), allocatable :: colored_text

        colored_text = apply_style_ucs4(text, set_ansi_code(fg_color=color))
    end function colorize_text_ucs4

    subroutine cursor_position_ucs4(row, col)
        integer, intent(in) :: row, col
        character(:, kind=ucs4), allocatable :: sequence
        character(:, kind=ucs4), allocatable :: temp

        write (temp, '(I0,A,I0)') row, delimiter, col
        sequence = CSI//temp//ucs4_'H'
        write (output_unit, '(A)', advance='no') sequence
    end subroutine cursor_position_ucs4

    subroutine clear_screen_ucs4()
        write (*, '(A)', advance='no') CSI//ucs4_'2J'
    end subroutine clear_screen_ucs4

    subroutine clear_line_ucs4()
        write (*, '(A)', advance='no') CSI//ucs4_'2K'
    end subroutine clear_line_ucs4

    subroutine save_cursor_ucs4()
        write (*, '(A)', advance='no') CSI//ucs4_'s'
    end subroutine save_cursor_ucs4

    subroutine restore_cursor_ucs4()
        write (*, '(A)', advance='no') CSI//ucs4_'u'
    end subroutine restore_cursor_ucs4

end submodule NAFPack_ANSI_ucs4
