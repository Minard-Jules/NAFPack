module test_meshgrid

    use NAFPack_kinds, only: i8, i16, isp, idp, sp, dp, qp
    use NAFPack_constant, only: TOL_TEST_sp, TOL_TEST_dp, TOL_TEST_qp
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NAFPack_meshgrid, only: meshgrid, INDEXING_IJ, INDEXING_XY
    use NAFPack_loop_method, only: LoopMethod, init_loop_method

    implicit none(type, external)

    private
    public :: collect_meshgrid_integer, collect_meshgrid_complex, collect_meshgrid_real

    interface
        module subroutine test_meshgrid_i8_2D(error)
        type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_i8_2D

        module subroutine test_meshgrid_i16_2D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_i16_2D

        module subroutine test_meshgrid_isp_2D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_isp_2D

        module subroutine test_meshgrid_idp_2D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_idp_2D
    end interface

    interface
        module subroutine test_meshgrid_i8_3D(error)
        type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_i8_3D

        module subroutine test_meshgrid_i16_3D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_i16_3D

        module subroutine test_meshgrid_isp_3D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_isp_3D

        module subroutine test_meshgrid_idp_3D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_idp_3D
    end interface

    interface
        module subroutine test_meshgrid_cmplx_sp_2D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_cmplx_sp_2D

        module subroutine test_meshgrid_cmplx_dp_2D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_cmplx_dp_2D

        module subroutine test_meshgrid_cmplx_qp_2D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_cmplx_qp_2D
    end interface

    interface
        module subroutine test_meshgrid_cmplx_sp_3D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_cmplx_sp_3D

        module subroutine test_meshgrid_cmplx_dp_3D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_cmplx_dp_3D

        module subroutine test_meshgrid_cmplx_qp_3D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_cmplx_qp_3D
    end interface

    interface 
        module subroutine test_meshgrid_real_sp_2D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_real_sp_2D

        module subroutine test_meshgrid_real_dp_2D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_real_dp_2D

        module subroutine test_meshgrid_real_qp_2D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_real_qp_2D
    end interface

    interface 
        module subroutine test_meshgrid_real_sp_3D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_real_sp_3D

        module subroutine test_meshgrid_real_dp_3D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_real_dp_3D

        module subroutine test_meshgrid_real_qp_3D(error)
            type(error_type), allocatable, intent(out) :: error
        end subroutine test_meshgrid_real_qp_3D
    end interface


contains

    subroutine collect_meshgrid_integer(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("meshgrid integer(i8)  2D", test_meshgrid_i8_2D), &
                    new_unittest("meshgrid integer(i16) 2D", test_meshgrid_i16_2D), &
                    new_unittest("meshgrid integer(isp) 2D", test_meshgrid_isp_2D), &
                    new_unittest("meshgrid integer(idp) 2D", test_meshgrid_idp_2D), &
                    new_unittest("meshgrid integer(i8)  3D", test_meshgrid_i8_3D), &
                    new_unittest("meshgrid integer(i16) 3D", test_meshgrid_i16_3D), &
                    new_unittest("meshgrid integer(isp) 3D", test_meshgrid_isp_3D), &
                    new_unittest("meshgrid integer(idp) 3D", test_meshgrid_idp_3D) &
                    ]

    end subroutine collect_meshgrid_integer

    subroutine collect_meshgrid_complex(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("meshgrid complex(sp) 2D", test_meshgrid_cmplx_sp_2D), &
                    new_unittest("meshgrid complex(dp) 2D", test_meshgrid_cmplx_dp_2D), &
                    new_unittest("meshgrid complex(qp) 2D", test_meshgrid_cmplx_qp_2D), &
                    new_unittest("meshgrid complex(sp) 3D", test_meshgrid_cmplx_sp_3D), &
                    new_unittest("meshgrid complex(dp) 3D", test_meshgrid_cmplx_dp_3D), &
                    new_unittest("meshgrid complex(qp) 3D", test_meshgrid_cmplx_qp_3D) &
                    ]

    end subroutine collect_meshgrid_complex

    subroutine collect_meshgrid_real(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("meshgrid real(sp) 2D", test_meshgrid_real_sp_2D), &
                    new_unittest("meshgrid real(dp) 2D", test_meshgrid_real_dp_2D), &
                    new_unittest("meshgrid real(qp) 2D", test_meshgrid_real_qp_2D), &
                    new_unittest("meshgrid real(sp) 3D", test_meshgrid_real_sp_3D), &
                    new_unittest("meshgrid real(dp) 3D", test_meshgrid_real_dp_3D), &
                    new_unittest("meshgrid real(qp) 3D", test_meshgrid_real_qp_3D) &
                    ]

    end subroutine collect_meshgrid_real

end module test_meshgrid

program test
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, init_color_output
    use test_meshgrid

    implicit none(type, external)

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("meshgrid integer", collect_meshgrid_integer), &
                 new_testsuite("meshgrid complex", collect_meshgrid_complex), &
                 new_testsuite("meshgrid real", collect_meshgrid_real) &
                 ]

    call init_color_output(.true.)

    do is = 1, size(testsuites)
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(s) meshgrid failed!"
        error stop
    end if

    deallocate(testsuites)

end program test
