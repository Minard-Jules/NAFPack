module NAFPack_utils

    use NAFPack_kinds, only: dp, sp, qp, i8, i16, isp, idp, ascii, ucs4
    use NAFPack_plan_type, only: DecimationMethod
    use NAFPack_plan_type, only: FFTAlgorithm
    use NAFPack_loop_method_type, only: LoopMethod
    use NAFPack_implementation_type, only: ImplementationType

    implicit none(type, external)

    private
    public :: present_arg

    interface present_arg
        procedure present_arg_logical

        procedure present_arg_integer_i8
        procedure present_arg_integer_i16
        procedure present_arg_integer_isp
        procedure present_arg_integer_idp

        procedure present_arg_real_sp
        procedure present_arg_real_dp
        procedure present_arg_real_qp

        procedure present_arg_complex_sp
        procedure present_arg_complex_dp
        procedure present_arg_complex_qp

        procedure present_arg_character_ascii
        procedure present_arg_character_ucs4

        procedure present_arg_type_DecimationMethod
        procedure present_arg_type_FFTAlgorithm
        procedure present_arg_type_LoopMethod
        procedure present_arg_type_ImplementationType
    end interface present_arg

    interface
        pure module function present_arg_integer_i8(arg, default_value) result(val)
            integer(i8), intent(in), optional :: arg
            integer(i8), intent(in) :: default_value
            integer(i8) :: val
        end function present_arg_integer_i8

        pure module function present_arg_integer_i16(arg, default_value) result(val)
            integer(i16), intent(in), optional :: arg
            integer(i16), intent(in) :: default_value
            integer(i16) :: val
        end function present_arg_integer_i16

        pure module function present_arg_integer_isp(arg, default_value) result(val)
            integer(isp), intent(in), optional :: arg
            integer(isp), intent(in) :: default_value
            integer(isp) :: val
        end function present_arg_integer_isp

        pure module function present_arg_integer_idp(arg, default_value) result(val)
            integer(idp), intent(in), optional :: arg
            integer(idp), intent(in) :: default_value
            integer(idp) :: val
        end function present_arg_integer_idp
    end interface

    interface
        pure module function present_arg_real_sp(arg, default_value) result(val)
            real(sp), intent(in), optional :: arg
            real(sp), intent(in) :: default_value
            real(sp) :: val
        end function present_arg_real_sp

        pure module function present_arg_real_dp(arg, default_value) result(val)
            real(dp), intent(in), optional :: arg
            real(dp), intent(in) :: default_value
            real(dp) :: val
        end function present_arg_real_dp

        pure module function present_arg_real_qp(arg, default_value) result(val)
            real(qp), intent(in), optional :: arg
            real(qp), intent(in) :: default_value
            real(qp) :: val
        end function present_arg_real_qp
    end interface

    interface
        pure module function present_arg_complex_sp(arg, default_value) result(val)
            complex(sp), intent(in), optional :: arg
            complex(sp), intent(in) :: default_value
            complex(sp) :: val
        end function present_arg_complex_sp

        pure module function present_arg_complex_dp(arg, default_value) result(val)
            complex(dp), intent(in), optional :: arg
            complex(dp), intent(in) :: default_value
            complex(dp) :: val
        end function present_arg_complex_dp

        pure module function present_arg_complex_qp(arg, default_value) result(val)
            complex(qp), intent(in), optional :: arg
            complex(qp), intent(in) :: default_value
            complex(qp) :: val
        end function present_arg_complex_qp
    end interface

    interface
        pure module function present_arg_type_DecimationMethod(arg, default_value) result(val)
            type(DecimationMethod), intent(in), optional :: arg
            type(DecimationMethod), intent(in) :: default_value
            type(DecimationMethod) :: val
        end function present_arg_type_DecimationMethod

        pure module function present_arg_type_FFTAlgorithm(arg, default_value) result(val)
            type(FFTAlgorithm), intent(in), optional :: arg
            type(FFTAlgorithm), intent(in) :: default_value
            type(FFTAlgorithm) :: val
        end function present_arg_type_FFTAlgorithm

        pure module function present_arg_type_LoopMethod(arg, default_value) result(val)
            type(LoopMethod), intent(in), optional :: arg
            type(LoopMethod), intent(in) :: default_value
            type(LoopMethod) :: val
        end function present_arg_type_LoopMethod

        pure module function present_arg_type_ImplementationType(arg, default_value) result(val)
            type(ImplementationType), intent(in), optional :: arg
            type(ImplementationType), intent(in) :: default_value
            type(ImplementationType) :: val
        end function present_arg_type_ImplementationType
    end interface

contains

    pure function present_arg_logical(arg, default_value) result(val)
        logical, intent(in), optional :: arg
        logical, intent(in) :: default_value
        logical :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_logical

    pure function present_arg_character_ascii(arg, default_value) result(val)
        character(len=*, kind=ascii), intent(in), optional :: arg
        character(len=*, kind=ascii), intent(in) :: default_value
        character(len=:, kind=ascii), allocatable :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_character_ascii

    pure function present_arg_character_ucs4(arg, default_value) result(val)
        character(len=*, kind=ucs4), intent(in), optional :: arg
        character(len=*, kind=ucs4), intent(in) :: default_value
        character(len=:, kind=ucs4), allocatable :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_character_ucs4

end module NAFPack_utils
