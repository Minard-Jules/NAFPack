module NAFPack_memory_monitor

    implicit none(type, external)

    private

    public :: get_memory_kb

    interface
        function get_memory_usage() bind(C, name="get_memory_usage") result(mem_usage)
            use iso_c_binding, only: c_int
            implicit none(type, external)
            integer(c_int) :: mem_usage
        end function get_memory_usage
    end interface

contains

    function get_memory_kb() result(memory_kb)
        integer :: memory_kb

        memory_kb = get_memory_usage()

    end function get_memory_kb

end module NAFPack_memory_monitor
