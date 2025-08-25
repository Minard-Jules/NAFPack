module FFTW3

    use, intrinsic :: iso_c_binding, only: C_INT32_T, C_INT, &
                                           C_INTPTR_T, C_PTR, C_FUNPTR, &
                                           C_FLOAT, C_DOUBLE, &
                                           C_FLOAT_COMPLEX, C_DOUBLE_COMPLEX, &
                                           C_CHAR, C_SIZE_T

    implicit none(type, external)

    public

    include '../../include/fftw-3.3.10/fftw3.f03'

end module FFTW3
