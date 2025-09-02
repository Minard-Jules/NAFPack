cmake_minimum_required(VERSION 3.15)

set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_CURRENT_SOURCE_DIR}/config/FlagsConfig.cmake)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" PROJECT_VERSION)
project(NAFPack
        LANGUAGES Fortran C
        VERSION ${PROJECT_VERSION}
        DESCRIPTION "NAFPack: Numerical Analysis in Fortran Package"
        HOMEPAGE_URL "https://github.com/Minard-Jules/NAFPack"
)

message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} ${CMAKE_Fortran_COMPILER}")
message(STATUS "C Compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION} ${CMAKE_C_COMPILER}")
message(STATUS "Linker: ${CMAKE_Fortran_COMPILER_LINKER_ID} ${CMAKE_Fortran_COMPILER_LINKER_VERSION} ${CMAKE_Fortran_COMPILER_LINKER}")
message(STATUS "System: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION} ${CMAKE_SYSTEM_PROCESSOR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# FFTW and OpenMP
message(STATUS "Searching for external OpenMP and FFTW3 libraries...")
find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

# GNU conventions for installation directories
include(GNUInstallDirs)

# source directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)

add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)
  # add_subdirectory(example)
endif()

install(EXPORT ${PROJECT_NAME}-targets
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)