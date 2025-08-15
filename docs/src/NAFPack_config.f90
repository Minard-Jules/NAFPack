!> Module for advanced configuration and tuning parameters
MODULE NAFPack_config

    USE NAFPack_constant

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE
    PUBLIC :: config_type, get_default_config, validate_config

    !> Configuration type for NAFPack
    TYPE :: config_type
        ! Numerical tolerances
        REAL(dp) :: pivot_tolerance = 1.0e-14_dp
        REAL(dp) :: convergence_tolerance = 1.0e-12_dp
        REAL(dp) :: residual_tolerance = 1.0e-10_dp

        ! Performance settings
        INTEGER :: max_iterations = 10000
        INTEGER :: block_size = 64
        LOGICAL :: use_openmp = .TRUE.
        LOGICAL :: use_blas = .TRUE.

        ! Memory management
        LOGICAL :: preallocate_workspace = .TRUE.
        INTEGER :: workspace_size = 1000

        ! Debugging and logging
        LOGICAL :: enable_debug = .FALSE.
        LOGICAL :: enable_timing = .FALSE.
        CHARACTER(LEN=100) :: log_file = "nafpack.log"

        ! Method selection
        CHARACTER(LEN=50) :: default_direct_method = "A_LU"
        CHARACTER(LEN=50) :: default_iterative_method = "Gauss_Seidel"
        CHARACTER(LEN=50) :: default_preconditioner = "ILU"
    END TYPE config_type

CONTAINS

    !> Get default configuration
    FUNCTION get_default_config() RESULT(config)
        TYPE(config_type) :: config
        ! Default values are already set in the type definition
        ! No need to explicitly set them as they're defined in the type
        config = config_type(pivot_tolerance=TOL_PIVOT, &
                             convergence_tolerance=TOL_CONVERGENCE, &
                             residual_tolerance=epsi, &
                             max_iterations=kmax, &
                             block_size=64, &
                             use_openmp=.TRUE., &
                             use_blas=.TRUE., &
                             preallocate_workspace=.TRUE., &
                             workspace_size=1000, &
                             enable_debug=.FALSE., &
                             enable_timing=.FALSE., &
                             log_file="nafpack.log", &
                             default_direct_method="A_LU", &
                             default_iterative_method="Gauss_Seidel", &
                             default_preconditioner="ILU")
    END FUNCTION get_default_config

    !> Validate configuration parameters
    SUBROUTINE validate_config(config, is_valid, error_msg)
        TYPE(config_type), INTENT(IN) :: config
        LOGICAL, INTENT(OUT) :: is_valid
        CHARACTER(LEN=*), INTENT(OUT) :: error_msg

        is_valid = .TRUE.
        error_msg = ""

        ! Check tolerances
        IF (config%pivot_tolerance <= 0.0_dp) THEN
            is_valid = .FALSE.
            error_msg = "Pivot tolerance must be positive"
            RETURN
        END IF

        IF (config%convergence_tolerance <= 0.0_dp) THEN
            is_valid = .FALSE.
            error_msg = "Convergence tolerance must be positive"
            RETURN
        END IF

        ! Check iteration limits
        IF (config%max_iterations <= 0) THEN
            is_valid = .FALSE.
            error_msg = "Maximum iterations must be positive"
            RETURN
        END IF

        ! Check block size
        IF (config%block_size <= 0) THEN
            is_valid = .FALSE.
            error_msg = "Block size must be positive"
            RETURN
        END IF

    END SUBROUTINE validate_config

END MODULE NAFPack_config
