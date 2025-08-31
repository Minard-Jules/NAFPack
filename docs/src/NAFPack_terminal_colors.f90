module NAFPack_terminal_colors
    
    use NAFPack_kinds, only: ascii, ucs4

    implicit none(type, external)

    public

    character(len=10), parameter :: red_color = char(27)//"[31m"
    character(len=10), parameter :: green_color = char(27)//"[32m"
    character(len=10), parameter :: yellow_color = char(27)//"[33m"
    character(len=10), parameter :: blue_color = char(27)//"[34m"
    character(len=10), parameter :: purple_color = char(27)//"[35m"
    character(len=10), parameter :: cyan_color = char(27)//"[36m"
    character(len=10), parameter :: white_color = char(27)//"[37m"
    character(len=10), parameter :: reset_color = char(27)//"[0m"

    character(KIND=ucs4, LEN=10), parameter :: red_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[31m"
    character(KIND=ucs4, LEN=10), parameter :: green_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[32m"
    character(KIND=ucs4, LEN=10), parameter :: yellow_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[33m"
    character(KIND=ucs4, LEN=10), parameter :: blue_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[34m"
    character(KIND=ucs4, LEN=10), parameter :: purple_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[35m"
    character(KIND=ucs4, LEN=10), parameter :: cyan_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[36m"
    character(KIND=ucs4, LEN=10), parameter :: white_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[37m"
    character(KIND=ucs4, LEN=10), parameter :: reset_color_ucs4 = char(27, KIND=ucs4)//ucs4_"[0m"

end module NAFPack_terminal_colors