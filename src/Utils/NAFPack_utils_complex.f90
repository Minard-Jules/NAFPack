submodule(NAFPack_utils) NAFPack_utils_complex

    implicit none(type, external)

contains

    pure module function present_arg_complex_sp(arg, default_value) result(val)
        complex(sp), intent(in), optional :: arg
        complex(sp), intent(in) :: default_value
        complex(sp) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_complex_sp

    pure module function present_arg_complex_dp(arg, default_value) result(val)
        complex(dp), intent(in), optional :: arg
        complex(dp), intent(in) :: default_value
        complex(dp) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_complex_dp

    pure module function present_arg_complex_qp(arg, default_value) result(val)
        complex(qp), intent(in), optional :: arg
        complex(qp), intent(in) :: default_value
        complex(qp) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_complex_qp

end submodule NAFPack_utils_complex
