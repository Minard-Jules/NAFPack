cmake_minimum_required(VERSION 3.15)
set(CMAKE_Fortran_COMPILER "gfortran" CACHE FILEPATH "")
project(NAFPack LANGUAGES Fortran)

# Option
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -O2")

# FFTW
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

# source directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)

# main source files
set(NAFPACK_SOURCES
    ${SRC_DIR}/NAFPack_constant.f90
    ${SRC_DIR}/Array_creation/NAFPack_meshgrid.f90
    ${SRC_DIR}/algebra/NAFPack_matricielle.f90
    ${SRC_DIR}/algebra/NAFPack_matrix_decomposition.f90
    ${SRC_DIR}/algebra/NAFPack_linear_algebra.f90
    ${SRC_DIR}/FFT/FFTW3.f90
    ${SRC_DIR}/FFT/NAFPack_fft.f90
)

# Test files
set(TEST_SOURCES
    ${TEST_DIR}/test_NAFPack_linear_algebra.f90
    ${TEST_DIR}/test_NAFPack_fft.f90
    ${TEST_DIR}/main_test.f90
)

# Main library
add_library(NAFPack STATIC ${NAFPACK_SOURCES})
target_include_directories(NAFPack PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_libraries(NAFPack PRIVATE PkgConfig::FFTW3 fftw3_threads)


# Test executable
add_executable(main_test ${TEST_SOURCES})
target_link_libraries(main_test PRIVATE NAFPack PkgConfig::FFTW3 fftw3_threads)
target_include_directories(main_test PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Installation (optional)
install(TARGETS NAFPack main_test
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
