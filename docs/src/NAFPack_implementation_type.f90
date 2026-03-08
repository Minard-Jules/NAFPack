module NAFPack_implementation_type
    use NAFPack_kinds, only: isp
    implicit none
    type :: ImplementationType
        integer(isp) :: id
        character(len=20) :: name
    end type ImplementationType

    type(ImplementationType), parameter :: RECURSIVE = ImplementationType(1, "Recursive"), &
                                           ITERATIVE = ImplementationType(2, "Iterative")
end module NAFPack_implementation_type