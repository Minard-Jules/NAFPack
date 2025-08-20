program test

    use NAFPack_constant
    use test_NAFPack_linear_algebra
    use test_NAFPack_fft

    implicit none(type, external)

    logical :: stat = .false.

    write (*, '(A)') purple_color, "Test linear systeme :", reset_color
    print*," "
    print*," "
    call test_linear_system(stat)
    print*," "

    write (*, '(A)') purple_color, "Test linear algebre :", reset_color
    print*," "
    print*," "
    call test_linear_algebra(stat)
    print*," "

    write (*, '(A)') purple_color, "Test FFT :", reset_color
    print*," "
    print*," "
    call test_FFT(stat)
    print*," "

    if (stat) then
        write (*, '(A)') red_color, "Test failed", reset_color
    else
        write (*, '(A)') green_color, "Test success", reset_color
    end if

end program test
