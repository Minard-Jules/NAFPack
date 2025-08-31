
module NAFPack_loop_method

    implicit none(type, external)

    private
    public :: LoopMethod, init_loop_method, count_true_methods
    public :: default_loop_method

    type :: ParallelMethod
        logical :: use_openmp = .false.
        logical :: use_mpi = .false.
        integer :: num_threads = 1
    end type ParallelMethod

    type :: LoopMethod
        logical :: use_do_classic = .false.
        logical :: use_vectorized = .false.
        logical :: use_do_concurrent = .false.
        type(ParallelMethod) :: parallel
    end type LoopMethod

    type(LoopMethod), parameter :: default_loop_method = LoopMethod(use_do_classic=.true.), &
                                   empty_loop_method = LoopMethod()

contains

    pure function init_loop_method( &
        use_do_classic, use_vectorized, use_do_concurrent, use_openmp, use_mpi, num_threads) result(loop_method)
        logical, intent(in), optional :: use_do_classic, use_vectorized, use_do_concurrent, use_openmp, use_mpi
        integer, intent(in), optional :: num_threads
        logical :: method_used
        type(LoopMethod) :: loop_method

        loop_method = empty_loop_method

        method_used = .false.
        if (present(use_do_classic)) then
            if (use_do_classic) loop_method%use_do_classic = .true.
                
        end if

        if (present(use_vectorized)) then
            if (use_vectorized) loop_method%use_vectorized = .true.
            call check_method_used(method_used)
        end if

        if (present(use_do_concurrent)) then
            if (use_do_concurrent) loop_method%use_do_concurrent = .true.
            call check_method_used(method_used)
        end if

        if (present(use_openmp)) then
            if (use_openmp) loop_method%parallel%use_openmp = .true.
            if (present(num_threads)) then
                if (num_threads > 0) then
                    loop_method%parallel%num_threads = num_threads
                else
                    error stop "num_threads must be a positive integer"
                end if
            end if
            call check_method_used(method_used)
        end if

        if (present(use_mpi)) then
            if (use_mpi) loop_method%parallel%use_mpi = .true.
            if (present(num_threads)) then
                if (num_threads > 0) then
                    loop_method%parallel%num_threads = num_threads
                else
                    error stop "num_threads must be a positive integer"
                end if
            end if
            call check_method_used(method_used)
        end if

        if (.not. method_used) then
            loop_method = default_loop_method
        end if

    end function init_loop_method

    pure subroutine check_method_used(method_used)
        logical, intent(inout) :: method_used

        if (.not. method_used) then
            method_used = .true.
        else
            error stop "Multiple loop methods cannot be used simultaneously"
        end if
    end subroutine

    pure function count_true_methods(loop_method) result(count_true)
        type(LoopMethod), intent(in) :: loop_method
        integer :: count_true

        count_true = 0
        if (loop_method%use_do_classic) count_true = count_true + 1
        if (loop_method%use_vectorized) count_true = count_true + 1
        if (loop_method%use_do_concurrent) count_true = count_true + 1
        if (loop_method%parallel%use_openmp) count_true = count_true + 1
        if (loop_method%parallel%use_mpi) count_true = count_true + 1
    end function count_true_methods

end module NAFPack_loop_method
