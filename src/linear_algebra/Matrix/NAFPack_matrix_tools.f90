MODULE NAFPack_matrix_tools 

    USE NAFPack_constant
    USE NAFPack_matricielle

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE Faddeev_Leverrier(A, c, Ainv, success, check)
        INTEGER, PARAMETER :: dp = KIND(1.0d0)
        REAL(dp), DIMENSION(:, :), INTENT(IN)  :: A
        LOGICAL, OPTIONAL, INTENT(IN) :: check
        REAL(dp), DIMENSION(:),     INTENT(OUT) :: c
        REAL(dp), DIMENSION(SIZE(A,1), SIZE(A,1)), OPTIONAL, INTENT(OUT) :: Ainv
        LOGICAL, OPTIONAL, INTENT(OUT) :: success
        REAL(dp), DIMENSION(SIZE(A,1), SIZE(A,1)) :: Bk, I, B_Nm1, AB
        LOGICAL :: do_check = .TRUE.
        INTEGER :: N, k

        N = SIZE(A,1)

        IF (PRESENT(check)) do_check = check

        IF (do_check) THEN
            PRINT*, "Checking if the matrix A is square and size of c is correct"
            IF (SIZE(A,2) /= N .OR. SIZE(c) < N+1) THEN
                PRINT *, "Error : Matrix A must be square and size of c must be at least N+1"
                STOP
            END IF
        END IF

        ! Initialization
        I = Identity_n(N)
        c = 0.0_dp
        c(1) = 1.0_dp
        c(2) = -Trace(A)
        Bk = A + c(2)*I

        DO k = 2, N
            AB = MATMUL(A, Bk)
            c(k+1) = -Trace(AB) / REAL(k, dp)
            Bk = AB + c(k+1)*I
            IF (k == N-1 .AND. PRESENT(Ainv)) B_Nm1 = -Bk
        END DO

        IF (PRESENT(Ainv) .AND. PRESENT(success)) THEN
            IF (ABS(c(N+1)) < 1.0e-12_dp) THEN
                success = .FALSE.
                Ainv = 0.0_dp
            ELSE
                success = .TRUE.
                Ainv = B_Nm1 / c(N+1)
            END IF
        ELSE IF (PRESENT(Ainv)) THEN
            IF (ABS(c(N+1)) < 1.0e-12_dp) THEN
                Ainv = 0.0_dp
            ELSE
                Ainv = B_Nm1 / c(N+1)
            END IF
        END IF

    END SUBROUTINE Faddeev_Leverrier

END MODULE NAFPack_matrix_tools