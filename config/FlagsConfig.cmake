if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(
        CMAKE_Fortran_FLAGS_INIT
        "-fimplicit-none"
        "-ffree-line-length-132"
    )
    set(
        CMAKE_Fortran_FLAGS_RELEASE_INIT
        "-O2"
        "-fimplicit-none"
        "-ffree-line-length-132"
    )
    set(
        CMAKE_Fortran_FLAGS_DEBUG_INIT
        "-Wall"
        "-Wextra"
        "-Wimplicit-procedure"
        "-Werror=implicit-procedure"
        "-std=f2018"
        "-g"
        "-fbacktrace"
    )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran")
endif()
string(REPLACE ";" " " CMAKE_Fortran_FLAGS_INIT "${CMAKE_Fortran_FLAGS_INIT}")
string(REPLACE ";" " " CMAKE_Fortran_FLAGS_RELEASE_INIT "${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
string(REPLACE ";" " " CMAKE_Fortran_FLAGS_DEBUG_INIT "${CMAKE_Fortran_FLAGS_DEBUG_INIT}")