submodule(NAFPack_utils) NAFPack_utils_real

    implicit none(type, external)

contains

    pure module function present_arg_real_sp(arg, default_value) result(val)
        real(sp), intent(in), optional :: arg
        real(sp), intent(in) :: default_value
        real(sp) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_real_sp

    pure module function present_arg_real_dp(arg, default_value) result(val)
        real(dp), intent(in), optional :: arg
        real(dp), intent(in) :: default_value
        real(dp) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_real_dp

    pure module function present_arg_real_qp(arg, default_value) result(val)
        real(qp), intent(in), optional :: arg
        real(qp), intent(in) :: default_value
        real(qp) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_real_qp
end submodule NAFPack_utils_real
