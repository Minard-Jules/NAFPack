module test_realloc

    use NAFPack_kinds, only: sp, dp, qp, i8, i16, isp, idp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_memory_management

    implicit none(type, external)

    private
    public :: collect_realloc_complex, collect_realloc_integer, collect_realloc_real

    interface
        module subroutine test_realloc_complex_sp_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_sp_1D_grow
        module subroutine test_realloc_complex_sp_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_sp_1D_shrink
        module subroutine test_realloc_complex_dp_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_dp_1D_grow
        module subroutine test_realloc_complex_dp_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_dp_1D_shrink
        module subroutine test_realloc_complex_qp_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_qp_1D_grow
        module subroutine test_realloc_complex_qp_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_qp_1D_shrink
    end interface
    interface
        module subroutine test_realloc_complex_sp_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_sp_2D_grow
        module subroutine test_realloc_complex_sp_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_sp_2D_shrink
        module subroutine test_realloc_complex_dp_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_dp_2D_grow
        module subroutine test_realloc_complex_dp_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_dp_2D_shrink
        module subroutine test_realloc_complex_qp_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_qp_2D_grow
        module subroutine test_realloc_complex_qp_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_qp_2D_shrink
    end interface
    interface
        module subroutine test_realloc_complex_sp_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_sp_3D_grow
        module subroutine test_realloc_complex_sp_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_sp_3D_shrink
        module subroutine test_realloc_complex_dp_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_dp_3D_grow
        module subroutine test_realloc_complex_dp_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_dp_3D_shrink
        module subroutine test_realloc_complex_qp_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_qp_3D_grow
        module subroutine test_realloc_complex_qp_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_complex_qp_3D_shrink
    end interface

    interface
        module subroutine test_realloc_integer_i8_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i8_1D_grow
        module subroutine test_realloc_integer_i8_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i8_1D_shrink
        module subroutine test_realloc_integer_i16_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i16_1D_grow
        module subroutine test_realloc_integer_i16_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i16_1D_shrink
        module subroutine test_realloc_integer_isp_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_isp_1D_grow
        module subroutine test_realloc_integer_isp_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_isp_1D_shrink
        module subroutine test_realloc_integer_idp_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_idp_1D_grow
        module subroutine test_realloc_integer_idp_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_idp_1D_shrink
    end interface
    interface
        module subroutine test_realloc_integer_i8_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i8_2D_grow
        module subroutine test_realloc_integer_i8_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i8_2D_shrink
        module subroutine test_realloc_integer_i16_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i16_2D_grow
        module subroutine test_realloc_integer_i16_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i16_2D_shrink
        module subroutine test_realloc_integer_isp_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_isp_2D_grow
        module subroutine test_realloc_integer_isp_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_isp_2D_shrink
        module subroutine test_realloc_integer_idp_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_idp_2D_grow
        module subroutine test_realloc_integer_idp_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_idp_2D_shrink
    end interface
    interface
        module subroutine test_realloc_integer_i8_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i8_3D_grow
        module subroutine test_realloc_integer_i8_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i8_3D_shrink
        module subroutine test_realloc_integer_i16_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i16_3D_grow
        module subroutine test_realloc_integer_i16_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_i16_3D_shrink
        module subroutine test_realloc_integer_isp_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_isp_3D_grow
        module subroutine test_realloc_integer_isp_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_isp_3D_shrink
        module subroutine test_realloc_integer_idp_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_idp_3D_grow
        module subroutine test_realloc_integer_idp_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_integer_idp_3D_shrink
    end interface

    interface
        module subroutine test_realloc_real_sp_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_sp_1D_grow
        module subroutine test_realloc_real_sp_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_sp_1D_shrink
        module subroutine test_realloc_real_dp_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_dp_1D_grow
        module subroutine test_realloc_real_dp_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_dp_1D_shrink
        module subroutine test_realloc_real_qp_1D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_qp_1D_grow
        module subroutine test_realloc_real_qp_1D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_qp_1D_shrink
    end interface
    interface
        module subroutine test_realloc_real_sp_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_sp_2D_grow
        module subroutine test_realloc_real_sp_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_sp_2D_shrink
        module subroutine test_realloc_real_dp_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_dp_2D_grow
        module subroutine test_realloc_real_dp_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_dp_2D_shrink
        module subroutine test_realloc_real_qp_2D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_qp_2D_grow
        module subroutine test_realloc_real_qp_2D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_qp_2D_shrink
    end interface
    interface
        module subroutine test_realloc_real_sp_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_sp_3D_grow
        module subroutine test_realloc_real_sp_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_sp_3D_shrink
        module subroutine test_realloc_real_dp_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_dp_3D_grow
        module subroutine test_realloc_real_dp_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_dp_3D_shrink
        module subroutine test_realloc_real_qp_3D_grow(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_qp_3D_grow
        module subroutine test_realloc_real_qp_3D_shrink(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_realloc_real_qp_3D_shrink
    end interface

contains

    subroutine collect_realloc_complex(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("realloc complex(sp) 1D grow", test_realloc_complex_sp_1D_grow), &
                    new_unittest("realloc complex(sp) 1D shrink", test_realloc_complex_sp_1D_shrink), &
                    new_unittest("realloc complex(dp) 1D grow", test_realloc_complex_dp_1D_grow), &
                    new_unittest("realloc complex(dp) 1D shrink", test_realloc_complex_dp_1D_shrink), &
                    new_unittest("realloc complex(qp) 1D grow", test_realloc_complex_qp_1D_grow), &
                    new_unittest("realloc complex(qp) 1D shrink", test_realloc_complex_qp_1D_shrink), &
                    new_unittest("realloc complex(sp) 2D grow", test_realloc_complex_sp_2D_grow), &
                    new_unittest("realloc complex(sp) 2D shrink", test_realloc_complex_sp_2D_shrink), &
                    new_unittest("realloc complex(dp) 2D grow", test_realloc_complex_dp_2D_grow), &
                    new_unittest("realloc complex(dp) 2D shrink", test_realloc_complex_dp_2D_shrink), &
                    new_unittest("realloc complex(qp) 2D grow", test_realloc_complex_qp_2D_grow), &
                    new_unittest("realloc complex(qp) 2D shrink", test_realloc_complex_qp_2D_shrink), &
                    new_unittest("realloc complex(sp) 3D grow", test_realloc_complex_sp_3D_grow), &
                    new_unittest("realloc complex(sp) 3D shrink", test_realloc_complex_sp_3D_shrink), &
                    new_unittest("realloc complex(dp) 3D grow", test_realloc_complex_dp_3D_grow), &
                    new_unittest("realloc complex(dp) 3D shrink", test_realloc_complex_dp_3D_shrink), &
                    new_unittest("realloc complex(qp) 3D grow", test_realloc_complex_qp_3D_grow), &
                    new_unittest("realloc complex(qp) 3D shrink", test_realloc_complex_qp_3D_shrink) &
                    ]
    end subroutine collect_realloc_complex

    subroutine collect_realloc_integer(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("realloc integer(i8)  1D grow", test_realloc_integer_i8_1D_grow), &
                    new_unittest("realloc integer(i8)  1D shrink", test_realloc_integer_i8_1D_shrink), &
                    new_unittest("realloc integer(i16) 1D grow", test_realloc_integer_i16_1D_grow), &
                    new_unittest("realloc integer(i16) 1D shrink", test_realloc_integer_i16_1D_shrink), &
                    new_unittest("realloc integer(isp) 1D grow", test_realloc_integer_isp_1D_grow), &
                    new_unittest("realloc integer(isp) 1D shrink", test_realloc_integer_isp_1D_shrink), &
                    new_unittest("realloc integer(idp) 1D grow", test_realloc_integer_idp_1D_grow), &
                    new_unittest("realloc integer(idp) 1D shrink", test_realloc_integer_idp_1D_shrink), &
                    new_unittest("realloc integer(i8)  2D grow", test_realloc_integer_i8_2D_grow), &
                    new_unittest("realloc integer(i8)  2D shrink", test_realloc_integer_i8_2D_shrink), &
                    new_unittest("realloc integer(i16) 2D grow", test_realloc_integer_i16_2D_grow), &
                    new_unittest("realloc integer(i16) 2D shrink", test_realloc_integer_i16_2D_shrink), &
                    new_unittest("realloc integer(isp) 2D grow", test_realloc_integer_isp_2D_grow), &
                    new_unittest("realloc integer(isp) 2D shrink", test_realloc_integer_isp_2D_shrink), &
                    new_unittest("realloc integer(idp) 2D grow", test_realloc_integer_idp_2D_grow), &
                    new_unittest("realloc integer(idp) 2D shrink", test_realloc_integer_idp_2D_shrink), &
                    new_unittest("realloc integer(i8)  3D grow", test_realloc_integer_i8_3D_grow), &
                    new_unittest("realloc integer(i8)  3D shrink", test_realloc_integer_i8_3D_shrink), &
                    new_unittest("realloc integer(i16) 3D grow", test_realloc_integer_i16_3D_grow), &
                    new_unittest("realloc integer(i16) 3D shrink", test_realloc_integer_i16_3D_shrink), &
                    new_unittest("realloc integer(isp) 3D grow", test_realloc_integer_isp_3D_grow), &
                    new_unittest("realloc integer(isp) 3D shrink", test_realloc_integer_isp_3D_shrink), &
                    new_unittest("realloc integer(idp) 3D grow", test_realloc_integer_idp_3D_grow), &
                    new_unittest("realloc integer(idp) 3D shrink", test_realloc_integer_idp_3D_shrink) &
                    ]
    end subroutine collect_realloc_integer

    subroutine collect_realloc_real(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("realloc real(sp) 1D grow", test_realloc_real_sp_1D_grow), &
                    new_unittest("realloc real(sp) 1D shrink", test_realloc_real_sp_1D_shrink), &
                    new_unittest("realloc real(dp) 1D grow", test_realloc_real_dp_1D_grow), &
                    new_unittest("realloc real(dp) 1D shrink", test_realloc_real_dp_1D_shrink), &
                    new_unittest("realloc real(qp) 1D grow", test_realloc_real_qp_1D_grow), &
                    new_unittest("realloc real(qp) 1D shrink", test_realloc_real_qp_1D_shrink), &
                    new_unittest("realloc real(sp) 2D grow", test_realloc_real_sp_2D_grow), &
                    new_unittest("realloc real(sp) 2D shrink", test_realloc_real_sp_2D_shrink), &
                    new_unittest("realloc real(dp) 2D grow", test_realloc_real_dp_2D_grow), &
                    new_unittest("realloc real(dp) 2D shrink", test_realloc_real_dp_2D_shrink), &
                    new_unittest("realloc real(qp) 2D grow", test_realloc_real_qp_2D_grow), &
                    new_unittest("realloc real(qp) 2D shrink", test_realloc_real_qp_2D_shrink), &
                    new_unittest("realloc real(sp) 3D grow", test_realloc_real_sp_3D_grow), &
                    new_unittest("realloc real(sp) 3D shrink", test_realloc_real_sp_3D_shrink), &
                    new_unittest("realloc real(dp) 3D grow", test_realloc_real_dp_3D_grow), &
                    new_unittest("realloc real(dp) 3D shrink", test_realloc_real_dp_3D_shrink), &
                    new_unittest("realloc real(qp) 3D grow", test_realloc_real_qp_3D_grow), &
                    new_unittest("realloc real(qp) 3D shrink", test_realloc_real_qp_3D_shrink) &
                    ]
    end subroutine collect_realloc_real

end module test_realloc

program test
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_realloc

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("realloc complex", collect_realloc_complex), &
                 new_testsuite("realloc integer", collect_realloc_integer), &
                 new_testsuite("realloc real", collect_realloc_real) &
                 ]

    call init_color_output(.true.)

    do is = 1, size(testsuites)
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(realloc) failed!"
        error stop
    end if

    deallocate(testsuites)

end program test
