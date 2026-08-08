submodule(NAFPack_utils) NAFPack_utils_integer

    implicit none(type, external)

contains

    pure module function present_arg_integer_i8(arg, default_value) result(val)
        integer(i8), intent(in), optional :: arg
        integer(i8), intent(in) :: default_value
        integer(i8) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_integer_i8

    pure module function present_arg_integer_i16(arg, default_value) result(val)
        integer(i16), intent(in), optional :: arg
        integer(i16), intent(in) :: default_value
        integer(i16) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_integer_i16

    pure module function present_arg_integer_isp(arg, default_value) result(val)
        integer(isp), intent(in), optional :: arg
        integer(isp), intent(in) :: default_value
        integer(isp) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_integer_isp

    pure module function present_arg_integer_idp(arg, default_value) result(val)
        integer(idp), intent(in), optional :: arg
        integer(idp), intent(in) :: default_value
        integer(idp) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_integer_idp

end submodule NAFPack_utils_integer
