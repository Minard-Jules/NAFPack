module NAFPack_loop_method

    use NAFPack_kinds, only: isp
    use NAFPack_utils, only: present_arg
    use NAFPack_loop_method_type, only: LoopMethod, default_loop_method, empty_loop_method

    implicit none(type, external)

    private
    public :: init_loop_method, count_true_methods, check_loop_method

contains

    pure function init_loop_method( &
        use_do_classic, &
        use_openmp_simd, &
        use_array_syntax, &
        use_do_concurrent, &
        use_openmp, &
        use_mpi, &
        num_threads) result(loop_method)
        logical, intent(in), optional :: use_do_classic, use_openmp_simd, use_array_syntax, &
                                         use_do_concurrent, use_openmp, use_mpi
        integer(isp), intent(in), optional :: num_threads
        type(LoopMethod) :: loop_method
        logical :: method_used

        loop_method = empty_loop_method
        method_used = .false.

        if (present_arg(use_do_classic, .false.)) call set_flag(loop_method%use_do_classic, method_used)
        if (present_arg(use_openmp_simd, .false.)) call set_flag(loop_method%vectorization%use_simd_openmp, method_used)
        if (present_arg(use_array_syntax, .false.)) call set_flag(loop_method%vectorization%use_array_syntax, method_used)
        if (present_arg(use_do_concurrent, .false.)) call set_flag(loop_method%use_do_concurrent, method_used)

        if (present_arg(use_openmp, .false.)) then
            call set_flag(loop_method%parallel%use_openmp, method_used)
            call apply_num_threads(loop_method%parallel%num_threads, num_threads)
        end if

        if (present_arg(use_mpi, .false.)) then
            call set_flag(loop_method%parallel%use_mpi, method_used)
            call apply_num_threads(loop_method%parallel%num_threads, num_threads)
        end if

        if (.not. method_used) loop_method = default_loop_method

    end function init_loop_method

    pure subroutine set_flag(flag, method_used)
        logical, intent(out) :: flag
        logical, intent(inout) :: method_used

        if (method_used) error stop "Multiple loop methods cannot be used simultaneously"
        flag = .true.
        method_used = .true.
    end subroutine set_flag

    pure subroutine apply_num_threads(num_threads_field, num_threads)
        integer(isp), intent(out) :: num_threads_field
        integer(isp), intent(in), optional :: num_threads

        if (present(num_threads)) then
            if (num_threads > 0) then
                num_threads_field = num_threads
            else
                error stop "num_threads must be a positive integer"
            end if
        end if
    end subroutine apply_num_threads

    pure function count_true_methods(loop_method) result(count_true)
        type(LoopMethod), intent(in) :: loop_method
        integer(isp) :: count_true

        count_true = merge(1, 0, loop_method%use_do_classic) &
                     + merge(1, 0, loop_method%vectorization%use_simd_openmp) &
                     + merge(1, 0, loop_method%vectorization%use_array_syntax) &
                     + merge(1, 0, loop_method%use_do_concurrent) &
                     + merge(1, 0, loop_method%parallel%use_openmp) &
                     + merge(1, 0, loop_method%parallel%use_mpi)
    end function count_true_methods

    pure function check_loop_method(loop_method) result(loop_method_used)
        type(LoopMethod), intent(in) :: loop_method
        type(LoopMethod) :: loop_method_used

        if (count_true_methods(loop_method) == 1) then
            loop_method_used = loop_method
        else
            loop_method_used = default_loop_method
        end if
    end function check_loop_method

end module NAFPack_loop_method
