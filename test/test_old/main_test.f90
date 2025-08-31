program test

    use NAFPack_constant
    use test_NAFPack_linear_algebra
    use test_NAFPack_fft
    use NAFPack_ANSI, only: ColorsAscii


    implicit none(type, external)

    logical :: stat = .false.
    type(ColorsAscii) :: colors

    CALL colors%init()

    write (*, '(A)') colors%magenta, "Test linear systeme :", colors%reset
    print*," "
    print*," "
    call test_linear_system(stat)
    print*," "

    write (*, '(A)') colors%magenta, "Test linear algebre :", colors%reset
    print*," "
    print*," "
    call test_linear_algebra(stat)
    print*," "

    write (*, '(A)') colors%magenta, "Test FFT :", colors%reset
    print*," "
    print*," "
    call test_FFT(stat)
    print*," "

    if (stat) then
        write (*, '(A)') colors%red, "Test failed", colors%reset
    else
        write (*, '(A)') colors%green, "Test success", colors%reset
    end if

end program test
