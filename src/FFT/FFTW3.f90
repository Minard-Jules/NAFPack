module FFTW3

    use, intrinsic :: iso_c_binding, only: &
        c_int32_t, c_int, c_intptr_t, c_ptr, c_funptr, c_float, c_double, &
        c_float_complex, c_double_complex, c_char, c_size_t

    implicit none(type, external)

    public

    include '../../include/fftw-3.3.10/fftw3.f03'

end module FFTW3
