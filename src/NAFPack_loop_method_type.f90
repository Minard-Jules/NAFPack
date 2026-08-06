module NAFPack_loop_method_type

    implicit none(type, external)

    private
    public :: LoopMethod, default_loop_method, empty_loop_method

    type :: ParallelMethod
        logical :: use_openmp = .false.
        logical :: use_mpi = .false.
        integer :: num_threads = 1
    end type ParallelMethod

    type :: VectorizationMethod
        logical :: use_simd_openmp = .false.
        logical :: use_array_syntax = .false.
    end type VectorizationMethod

    type :: LoopMethod
        logical :: use_do_classic = .false.
        logical :: use_array_syntax = .false.
        logical :: use_do_concurrent = .false.
        type(ParallelMethod) :: parallel
        type(VectorizationMethod) :: vectorization
    end type LoopMethod

    type(LoopMethod), parameter :: default_loop_method = LoopMethod(use_do_classic=.true.), &
                                   empty_loop_method = LoopMethod()

end module NAFPack_loop_method_type
