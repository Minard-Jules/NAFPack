!> Module for advanced configuration and tuning parameters
module NAFPack_config

    use NAFPack_constant, only: dp, TOL_CONVERGENCE, TOL_PIVOT, epsi, kmax

    implicit none(type, external)

    private
    public :: config_type, get_default_config, validate_config

    !> Configuration type for NAFPack
    type :: config_type
        ! Numerical tolerances
        real(dp) :: pivot_tolerance = 1.0e-14_dp
        real(dp) :: convergence_tolerance = 1.0e-12_dp
        real(dp) :: residual_tolerance = 1.0e-10_dp

        ! Performance settings
        integer :: max_iterations = 10000
        integer :: block_size = 64
        logical :: use_openmp = .true.
        logical :: use_blas = .true.

        ! Memory management
        logical :: preallocate_workspace = .true.
        integer :: workspace_size = 1000

        ! Debugging and logging
        logical :: enable_debug = .false.
        logical :: enable_timing = .false.
        character(LEN=100) :: log_file = "nafpack.log"

        ! Method selection
        character(LEN=50) :: default_direct_method = "A_LU"
        character(LEN=50) :: default_iterative_method = "Gauss_Seidel"
        character(LEN=50) :: default_preconditioner = "ILU"
    end type config_type

contains

    !> Get default configuration
    function get_default_config() result(config)
        type(config_type) :: config
        ! Default values are already set in the type definition
        ! No need to explicitly set them as they're defined in the type
        config = config_type(pivot_tolerance=TOL_PIVOT, &
                             convergence_tolerance=TOL_CONVERGENCE, &
                             residual_tolerance=epsi, &
                             max_iterations=kmax, &
                             block_size=64, &
                             use_openmp=.true., &
                             use_blas=.true., &
                             preallocate_workspace=.true., &
                             workspace_size=1000, &
                             enable_debug=.false., &
                             enable_timing=.false., &
                             log_file="nafpack.log", &
                             default_direct_method="A_LU", &
                             default_iterative_method="Gauss_Seidel", &
                             default_preconditioner="ILU")
    end function get_default_config

    !> Validate configuration parameters
    subroutine validate_config(config, is_valid, error_msg)
        type(config_type), intent(in) :: config
        logical, intent(out) :: is_valid
        character(LEN=40), intent(out) :: error_msg

        is_valid = .true.
        error_msg = ""

        ! Check tolerances
        if (config%pivot_tolerance <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Pivot tolerance must be positive"
            return
        end if

        if (config%convergence_tolerance <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Convergence tolerance must be positive"
            return
        end if

        ! Check iteration limits
        if (config%max_iterations <= 0) then
            is_valid = .false.
            error_msg = "Maximum iterations must be positive"
            return
        end if

        ! Check block size
        if (config%block_size <= 0) then
            is_valid = .false.
            error_msg = "Block size must be positive"
            return
        end if

    end subroutine validate_config

end module NAFPack_config
