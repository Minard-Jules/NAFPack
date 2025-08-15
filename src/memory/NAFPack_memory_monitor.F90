MODULE NAFPack_memory_monitor

    USE iso_c_binding

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE

    PUBLIC :: get_memory_kb

#ifdef _WIN32
    INTERFACE
        SUBROUTINE get_windows_memory(mem_kb) BIND(C, name="get_windows_memory")
            IMPORT c_int
            INTEGER(c_int), INTENT(OUT) :: mem_kb
        END SUBROUTINE
    END INTERFACE
#else
    INTERFACE
        SUBROUTINE get_linux_memory(mem_kb) BIND(C, NAME="get_linux_memory")
            IMPORT c_int
            INTEGER(c_int), INTENT(OUT) :: mem_kb
        END SUBROUTINE
    END INTERFACE
#endif

CONTAINS

    FUNCTION get_memory_kb() RESULT(memory_kb)
        INTEGER :: memory_kb

        memory_kb = 0
        
#ifdef _WIN32
        CALL get_windows_memory(memory_kb)
#else
        CALL get_linux_memory(memory_kb)
#endif

    END FUNCTION get_memory_kb

END MODULE NAFPack_memory_monitor
