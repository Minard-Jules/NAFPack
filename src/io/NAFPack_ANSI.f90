module NAFPack_ANSI

    use, intrinsic :: iso_fortran_env, only: output_unit
    use NAFPack_kinds, only: i8, i16, ascii, ucs4
    use NAFPack_io_utils, only: to_str_ascii, to_str_ucs4

    implicit none(type, external)

    private

    public :: output_unit
    public :: AnsiCode
    public :: Ansi_Constants
    public :: ColorsAscii, ColorsUcs4

    public :: set_ansi_code, reset_ansi_code

    public :: create_ansi_ascii, create_ansi_ucs4
    public :: apply_style_ascii, apply_style_ucs4
    public :: colorize_text_ascii, colorize_text_ucs4
    public :: cursor_position_ascii, cursor_position_ucs4
    public :: clear_screen_ascii, clear_screen_ucs4
    public :: clear_line_ascii, clear_line_ucs4
    public :: save_cursor_ascii, save_cursor_ucs4
    public :: restore_cursor_ascii, restore_cursor_ucs4

    type :: AnsiCode
        private
        !> Style descriptor
        logical :: use_style = .false.
        integer(i8) :: style = -1_i8
        !> Foreground color
        logical :: use_fg = .false.
        integer(i16) :: fg = -1_i16
        !> Background color
        logical :: use_bg = .false.
        integer(i16) :: bg = -1_i16
    end type AnsiCode

    type :: AnsiConstants 
        ! Styles
        type(AnsiCode) :: &
            STYLE_RESET = AnsiCode(use_style=.true., style=0_i8), &
            STYLE_BOLD = AnsiCode(use_style=.true., style=1_i8), &
            STYLE_FAINT = AnsiCode(use_style=.true., style=2_i8), &
            STYLE_ITALIC = AnsiCode(use_style=.true., style=3_i8), &
            STYLE_UNDERLINE = AnsiCode(use_style=.true., style=4_i8), &
            STYLE_BLINK = AnsiCode(use_style=.true., style=5_i8), &
            STYLE_REVERSE = AnsiCode(use_style=.true., style=7_i8), &
            STYLE_HIDDEN = AnsiCode(use_style=.true., style=8_i8), &
            STYLE_STRIKETHROUGH = AnsiCode(use_style=.true., style=9_i8)

        ! fg colors
        type(AnsiCode) :: &
            FG_BLACK = AnsiCode(use_fg=.true., fg=0_i16), &
            FG_RED = AnsiCode(use_fg=.true., fg=1_i16), &
            FG_GREEN = AnsiCode(use_fg=.true., fg=2_i16), &
            FG_YELLOW = AnsiCode(use_fg=.true., fg=3_i16), &
            FG_BLUE = AnsiCode(use_fg=.true., fg=4_i16), &
            FG_MAGENTA = AnsiCode(use_fg=.true., fg=5_i16), &
            FG_CYAN = AnsiCode(use_fg=.true., fg=6_i16), &
            FG_WHITE = AnsiCode(use_fg=.true., fg=7_i16)

        ! fg bright colors
        type(AnsiCode) :: &
            FG_BRIGHT_BLACK = AnsiCode(use_fg=.true., fg=8_i16), &
            FG_BRIGHT_RED = AnsiCode(use_fg=.true., fg=9_i16), &
            FG_BRIGHT_GREEN = AnsiCode(use_fg=.true., fg=10_i16), &
            FG_BRIGHT_YELLOW = AnsiCode(use_fg=.true., fg=11_i16), &
            FG_BRIGHT_BLUE = AnsiCode(use_fg=.true., fg=12_i16), &
            FG_BRIGHT_MAGENTA = AnsiCode(use_fg=.true., fg=13_i16), &
            FG_BRIGHT_CYAN = AnsiCode(use_fg=.true., fg=14_i16), &
            FG_BRIGHT_WHITE = AnsiCode(use_fg=.true., fg=15_i16)

        ! bg colors
        type(AnsiCode) :: &
            BG_BLACK = AnsiCode(use_bg=.true., bg=0_i16), &
            BG_RED = AnsiCode(use_bg=.true., bg=1_i16), &
            BG_GREEN = AnsiCode(use_bg=.true., bg=2_i16), &
            BG_YELLOW = AnsiCode(use_bg=.true., bg=3_i16), &
            BG_BLUE = AnsiCode(use_bg=.true., bg=4_i16), &
            BG_MAGENTA = AnsiCode(use_bg=.true., bg=5_i16), &
            BG_CYAN = AnsiCode(use_bg=.true., bg=6_i16), &
            BG_WHITE = AnsiCode(use_bg=.true., bg=7_i16)

        ! bg bright colors
        type(AnsiCode) :: &
            BG_BRIGHT_BLACK = AnsiCode(use_bg=.true., bg=8_i16), &
            BG_BRIGHT_RED = AnsiCode(use_bg=.true., bg=9_i16), &
            BG_BRIGHT_GREEN = AnsiCode(use_bg=.true., bg=10_i16), &
            BG_BRIGHT_YELLOW = AnsiCode(use_bg=.true., bg=11_i16), &
            BG_BRIGHT_BLUE = AnsiCode(use_bg=.true., bg=12_i16), &
            BG_BRIGHT_MAGENTA = AnsiCode(use_bg=.true., bg=13_i16), &
            BG_BRIGHT_CYAN = AnsiCode(use_bg=.true., bg=14_i16), &
            BG_BRIGHT_WHITE = AnsiCode(use_bg=.true., bg=15_i16)
    end type AnsiConstants

    type(AnsiConstants), parameter :: Ansi_Constants = AnsiConstants()

    type :: ColorsAscii
        character(len=:), allocatable :: &
            reset, bold, faint, italic, underline, blink, reverse, hidden, strikethrough, &
            red, green, yellow, blue, magenta, cyan, white, &
            bright_red, bright_green, bright_yellow, bright_blue, bright_magenta, bright_cyan, &
            bright_white
        contains
            procedure :: init => init_colors_ascii
    end type ColorsAscii

    interface
        module subroutine init_colors_ascii(this)
            class(ColorsAscii), intent(out) :: this
        end subroutine init_colors_ascii
    end interface

    type :: ColorsUcs4
        character(len=:, kind=ucs4), allocatable :: &
            reset, bold, faint, italic, underline, blink, reverse, hidden, strikethrough, &
            red, green, yellow, blue, magenta, cyan, white, &
            bright_red, bright_green, bright_yellow, bright_blue, bright_magenta, bright_cyan, &
            bright_white
        contains
            procedure :: init => init_colors_ucs4
    end type ColorsUcs4

    interface
        module subroutine init_colors_ucs4(this)
            class(ColorsUcs4), intent(out) :: this
        end subroutine init_colors_ucs4
    end interface

    interface
        pure module function create_ansi_ascii(ansi_code) result(ansi_string)
            type(AnsiCode), intent(in) :: ansi_code
            character(len=:, kind=ascii), allocatable :: ansi_string
        end function create_ansi_ascii

        pure module function apply_style_ascii(text, style) result(styled_text)
            character(*, kind=ascii), intent(in) :: text
            type(AnsiCode), intent(in) :: style
            character(:, kind=ascii), allocatable :: styled_text
        end function apply_style_ascii

        pure module function colorize_text_ascii(text, color) result(colored_text)
            character(*, kind=ascii), intent(in) :: text
            integer(i16), intent(in) :: color
            character(:, kind=ascii), allocatable :: colored_text
        end function colorize_text_ascii

        module subroutine cursor_position_ascii(row, col)
            integer, intent(in) :: row, col
        end subroutine cursor_position_ascii

        module subroutine clear_screen_ascii()
        end subroutine clear_screen_ascii

        module subroutine clear_line_ascii()
        end subroutine clear_line_ascii

        module subroutine save_cursor_ascii()
        end subroutine save_cursor_ascii

        module subroutine restore_cursor_ascii()
        end subroutine restore_cursor_ascii
    end interface

    interface
        pure module function create_ansi_ucs4(ansi_code) result(ansi_string)
            type(AnsiCode), intent(in) :: ansi_code
            character(len=:, kind=ucs4), allocatable :: ansi_string
        end function create_ansi_ucs4

        pure module function apply_style_ucs4(text, style) result(styled_text)
            character(*, kind=ucs4), intent(in) :: text
            type(AnsiCode), intent(in) :: style
            character(:, kind=ucs4), allocatable :: styled_text
        end function apply_style_ucs4

        pure module function colorize_text_ucs4(text, color) result(colored_text)
            character(*, kind=ucs4), intent(in) :: text
            integer(i16), intent(in) :: color
            character(:, kind=ucs4), allocatable :: colored_text
        end function colorize_text_ucs4

        module subroutine cursor_position_ucs4(row, col)
            integer, intent(in) :: row, col
        end subroutine cursor_position_ucs4

        module subroutine clear_screen_ucs4()
        end subroutine clear_screen_ucs4

        module subroutine clear_line_ucs4()
        end subroutine clear_line_ucs4

        module subroutine save_cursor_ucs4()
        end subroutine save_cursor_ucs4

        module subroutine restore_cursor_ucs4()
        end subroutine restore_cursor_ucs4
    end interface

contains

    pure module function set_ansi_code(style, fg_color, bg_color) result(ansi_code)
        integer(i8), optional, intent(in) :: style
        integer(i16), optional, intent(in) :: fg_color, bg_color
        type(AnsiCode) :: ansi_code

        ansi_code%use_style = .false.
        ansi_code%use_fg = .false.
        ansi_code%use_bg = .false.

        if (present(style)) then
            ansi_code%use_style = style >= 0 .and. style < 10
            if (ansi_code%use_style) ansi_code%style = style
        end if

        if (present(fg_color)) then
            ansi_code%use_fg = fg_color >= 0 .and. fg_color <= 255
            if (ansi_code%use_fg) ansi_code%fg = fg_color
        end if

        if (present(bg_color)) then
            ansi_code%use_bg = bg_color >= 0 .and. bg_color <= 255
            if (ansi_code%use_bg) ansi_code%bg = bg_color
        end if

        if (.not. (ansi_code%use_style .or. ansi_code%use_fg .or. ansi_code%use_bg)) then
            ansi_code = Ansi_Constants%STYLE_RESET
        end if

    end function set_ansi_code

    pure module function reset_ansi_code() result(ansi_code)
        type(AnsiCode) :: ansi_code
        ansi_code = Ansi_Constants%STYLE_RESET
    end function reset_ansi_code
end module NAFPack_ANSI
