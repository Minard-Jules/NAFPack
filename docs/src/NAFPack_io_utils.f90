module NAFPack_io_utils

    use NAFPack_kinds, only: ascii, ucs4, &
                             sp, dp, qp, &
                             i8, i16, isp, idp

    implicit none(type, external)

    private

    public :: to_str_ascii, to_str_ucs4

    interface to_str_ascii
        module procedure int8_to_str_ascii
        module procedure int16_to_str_ascii
    end interface to_str_ascii

    interface to_str_ucs4
        module procedure int8_to_str_ucs4
        module procedure int16_to_str_ucs4
    end interface to_str_ucs4

contains

    pure function int8_to_str_ascii(value) result(str)
        integer(i8), intent(in) :: value
        character(len=:, kind=ascii), allocatable :: str
        character(len=64, kind=ascii) :: buffer
        integer :: nlen

        write (buffer, '(I0)') value

        nlen = len_trim(buffer)

        allocate (character(len=nlen, kind=ascii) :: str)

        str = buffer(1:nlen)
    end function int8_to_str_ascii

    pure function int16_to_str_ascii(value) result(str)
        integer(i16), intent(in) :: value
        character(len=:, kind=ascii), allocatable :: str
        character(len=64, kind=ascii) :: buffer
        integer :: nlen

        write (buffer, '(I0)') value

        nlen = len_trim(buffer)

        allocate (character(len=nlen, kind=ascii) :: str)

        str = buffer(1:nlen)
    end function int16_to_str_ascii

    pure function int8_to_str_ucs4(value) result(str)
        integer(i8), intent(in) :: value
        character(len=:, kind=ucs4), allocatable :: str
        character(len=64, kind=ucs4) :: buffer
        integer :: nlen

        write (buffer, '(I0)') value

        nlen = len_trim(buffer)

        allocate (character(len=nlen, kind=ucs4) :: str)

        str = buffer(1:nlen)
    end function int8_to_str_ucs4

    pure function int16_to_str_ucs4(value) result(str)
        integer(i16), intent(in) :: value
        character(len=:, kind=ucs4), allocatable :: str
        character(len=64, kind=ucs4) :: buffer
        integer :: nlen

        write (buffer, '(I0)') value

        nlen = len_trim(buffer)

        allocate (character(len=nlen, kind=ucs4) :: str)

        str = buffer(1:nlen)
    end function int16_to_str_ucs4

end module NAFPack_io_utils
