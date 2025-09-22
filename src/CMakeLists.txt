set(NAFPACK_SOURCES
    ${SRC_DIR}/NAFPack_kinds.f90
    ${SRC_DIR}/NAFPack_constant.f90
    ${SRC_DIR}/NAFPack_loop_method.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_real_2D_compute.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_real_2D.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_real_3D_compute.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_real_3D.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_integer_2D_compute.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_integer_2D.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_integer_3D_compute.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_integer_3D.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_complex_2D_compute.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_complex_2D.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_complex_3D_compute.f90
    ${SRC_DIR}/Array_creation/meshgrid/NAFPack_meshgrid_complex_3D.f90
    ${SRC_DIR}/Fourier_Transform/FFTW3.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_fft.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_Fourier_Transform_dft_compute.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_Fourier_Transform_dft.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_Fourier_Transform_dft2.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_Fourier_Transform_dft3.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_Fourier_Transform_fft_compute_mixed_radix.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_Fourier_Transform_fft_compute_radix2.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_Fourier_Transform_fft_compute_split_radix.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_Fourier_Transform_fft.f90
    ${SRC_DIR}/Fourier_Transform/NAFPack_Fourier_Transform.f90
    ${SRC_DIR}/io/NAFPack_ANSI_ASCII.f90
    ${SRC_DIR}/io/NAFPack_ANSI_ucs4.f90
    ${SRC_DIR}/io/NAFPack_ANSI.f90
    ${SRC_DIR}/io/NAFPack_io_utils.f90
    ${SRC_DIR}/io/NAFPack_Logger_mod.f90
    ${SRC_DIR}/io/NAFPack_matrix_market.f90
    ${SRC_DIR}/io/NAFPack_terminal.f90
    ${SRC_DIR}/linear_algebra/NAFPack_linalg.f90
    ${SRC_DIR}/linear_algebra/NAFPack_Preconditioners.f90
    ${SRC_DIR}/linear_algebra/Direct_method/NAFPack_Direct_methode.f90
    ${SRC_DIR}/linear_algebra/Direct_method/NAFPack_Direct_types.f90
    ${SRC_DIR}/linear_algebra/Eigen/NAFPack_Eigen.f90
    ${SRC_DIR}/linear_algebra/Iterative_method/NAFPack_Iterative_methods.f90
    ${SRC_DIR}/linear_algebra/Iterative_method/NAFPack_Iterative_Params.f90
    ${SRC_DIR}/linear_algebra/Iterative_method/NAFPack_Iterative_types.f90
    ${SRC_DIR}/linear_algebra/Iterative_method/NAFPack_Krylov_method.f90
    ${SRC_DIR}/linear_algebra/Matrix/NAFPack_matricielle.f90
    ${SRC_DIR}/linear_algebra/Matrix/NAFPack_matrix_decomposition.f90
    ${SRC_DIR}/linear_algebra/Matrix/NAFPack_matrix_properties.f90
    ${SRC_DIR}/linear_algebra/Matrix/NAFPack_matrix_tools.f90
    ${SRC_DIR}/math_utils/NAFPack_math_utils.f90
    ${SRC_DIR}/memory/memory_monitor.c
    ${SRC_DIR}/memory/NAFPack_memory_management_complex.f90
    ${SRC_DIR}/memory/NAFPack_memory_management_integer.f90
    ${SRC_DIR}/memory/NAFPack_memory_management_real.f90
    ${SRC_DIR}/memory/NAFPack_memory_management.f90
    ${SRC_DIR}/memory/NAFPack_memory_monitor.f90
)

add_library(${PROJECT_NAME} ${NAFPACK_SOURCES})

set_target_properties(${PROJECT_NAME}
                      PROPERTIES
                      POSITION_INDEPENDENT_CODE ON
                      WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set_target_properties(${PROJECT_NAME}
                      PROPERTIES
                      Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_include_directories(${PROJECT_NAME} PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_Fortran)
target_link_libraries(${PROJECT_NAME} PkgConfig::FFTW3 fftw3_threads)

install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}-targets
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION "${CMAKE_INSTALL_MODULEDIR}")