module NAFPack_terminal

    use iso_fortran_env, only: output_unit

    
    use NAFPack_kinds, only: ascii, ucs4

    ! ASCII colors for the terminal
    use NAFPack_terminal_colors, only: &
        red_color, green_color, yellow_color, blue_color, &
        white_color, cyan_color, purple_color, reset_color

    ! UCS4 colors for the terminal (Unicode)
    use NAFPack_terminal_colors, only: &
        red_color_ucs4, green_color_ucs4, yellow_color_ucs4, blue_color_ucs4, &
        white_color_ucs4, cyan_color_ucs4, purple_color_ucs4, reset_color_ucs4

    implicit none(type, external)

    interface print_colored_message
        module procedure print_colored_message_ucs4
        module procedure print_colored_message_ascii
    end interface print_colored_message

    private

    public :: output_unit, print_colored_message

contains

    subroutine print_colored_message_ascii(color, message)
        character(kind=ascii, len=*), intent(in) :: message
        character(kind=ascii, len=*), intent(in) :: color

        open (output_unit, encoding='DEFAULT')

        write (output_unit, '(A)', advance='no') trim(color)
        write (output_unit, '(A)', advance='no') trim(message)
        write (output_unit, '(A)', advance='yes') trim(reset_color)

        close (output_unit)

    end subroutine print_colored_message_ascii

    subroutine print_colored_message_ucs4(color, message)
        character(kind=ucs4, len=*), intent(in) :: message
        character(kind=ucs4, len=*), intent(in) :: color

        open (output_unit, encoding='UTF-8')

        write (output_unit, '(A)', advance='no') trim(color)
        write (output_unit, '(A)', advance='no') trim(message)
        write (output_unit, '(A)', advance='yes') trim(reset_color_ucs4)

        close (output_unit)

    end subroutine print_colored_message_ucs4

end module NAFPack_terminal
