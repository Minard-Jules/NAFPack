submodule(NAFPack_utils) NAFPack_utils_type

    implicit none(type, external)

contains

    pure module function present_arg_type_DecimationMethod(arg, default_value) result(val)
        type(DecimationMethod), intent(in), optional :: arg
        type(DecimationMethod), intent(in) :: default_value
        type(DecimationMethod) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_type_DecimationMethod

    pure module function present_arg_type_FFTAlgorithm(arg, default_value) result(val)
        type(FFTAlgorithm), intent(in), optional :: arg
        type(FFTAlgorithm), intent(in) :: default_value
        type(FFTAlgorithm) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_type_FFTAlgorithm

    pure module function present_arg_type_LoopMethod(arg, default_value) result(val)
        type(LoopMethod), intent(in), optional :: arg
        type(LoopMethod), intent(in) :: default_value
        type(LoopMethod) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_type_LoopMethod

    pure module function present_arg_type_ImplementationType(arg, default_value) result(val)
        type(ImplementationType), intent(in), optional :: arg
        type(ImplementationType), intent(in) :: default_value
        type(ImplementationType) :: val

        if (present(arg)) then
            val = arg
        else
            val = default_value
        end if
    end function present_arg_type_ImplementationType

end submodule NAFPack_utils_type
