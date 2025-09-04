set(TEST_SOURCES_old
    ${TEST_DIR}/test_old/test_old.f90
    ${TEST_DIR}/test_old/test_NAFPack_fft.f90
    ${TEST_DIR}/test_old/test_NAFPack_linear_algebra.f90
)

include(FetchContent)

# Downloader et compiler test-drive
set(TEST_DRIVE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  test-drive
  GIT_REPOSITORY https://github.com/fortran-lang/test-drive.git
  GIT_TAG        main
)
FetchContent_MakeAvailable(test-drive)

macro(ADDTEST name)
    add_executable(test_${name} test_${name}.f90)
    target_link_libraries(test_${name} "${PROJECT_NAME}" "test-drive")
    add_test(NAME ${name}
             COMMAND $<TARGET_FILE:test_${name}> ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADDTEST)

add_executable(test_old ${TEST_SOURCES_old})
target_link_libraries(test_old "${PROJECT_NAME}" "test-drive")
# target_include_directories(test_old PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

add_subdirectory(Array_creation)
add_subdirectory(Fourier_Transform)