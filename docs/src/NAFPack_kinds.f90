module NAFPack_kinds

    use iso_fortran_env, only: int8, int16, int32, int64, &
                               real32, real64, real128

    implicit none(type, external)

    integer, parameter :: ascii = selected_char_kind('ascii')
    integer, parameter :: ucs4 = selected_char_kind('ISO_10646')

    integer, parameter :: sp = real32
    integer, parameter :: dp = real64
    integer, parameter :: qp = real128

    integer, parameter :: i8 = int8
    integer, parameter :: i16 = int16
    integer, parameter :: isp = int32
    integer, parameter :: idp = int64

    public

end module NAFPack_kinds
