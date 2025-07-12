MODULE NAFPack_matrix_properties

    USE NAFPack_constant
    USE NAFPack_matricielle
    USE NAFPack_Eigen

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: is_square_matrix, is_symmetric, is_orthogonal, is_SPD, is_tridiagonal, &
              is_diagonally_dominant, is_non_zero_diagonal

    CONTAINS

    FUNCTION is_square_matrix(A) RESULT(is_square)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL :: is_square

        PRINT *, "Checking if the matrix is square..."

        is_square = (SIZE(A, 1) == SIZE(A, 2))

    END FUNCTION is_square_matrix

    FUNCTION is_symmetric(A, is_square) RESULT(is_sym)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL, OPTIONAL, INTENT(IN) :: is_square
        LOGICAL :: is_sym

        IF (PRESENT(is_square)) THEN
            IF (.NOT. is_square) THEN
                is_sym = .FALSE.
                RETURN
            END IF
        ELSE IF (.NOT. is_square_matrix(A)) THEN
            is_sym = .FALSE.
            RETURN
        END IF

        PRINT *, "Checking if the matrix is symmetric..."

        is_sym = ALL(A == TRANSPOSE(A))

    END FUNCTION is_symmetric

    FUNCTION is_orthogonal(A, is_square) RESULT(is_orth)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL, OPTIONAL, INTENT(IN) :: is_square
        LOGICAL :: is_orth

        IF (PRESENT(is_square)) THEN
            IF (.NOT. is_square) THEN
                is_orth = .FALSE.
                RETURN
            END IF
        ELSE IF (.NOT. is_square_matrix(A)) THEN
            is_orth = .FALSE.
            RETURN
        END IF

        PRINT *, "Checking if the matrix is orthogonal..."

        is_orth = ALL(ABS(MATMUL(A, TRANSPOSE(A)) - Identity_n(SIZE(A, 1))) < epsi_test)

    END FUNCTION is_orthogonal

    FUNCTION is_SPD(A, is_sym) RESULT(is_spd_matrix)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL, OPTIONAL, INTENT(IN) :: is_sym
        REAL(dp), DIMENSION(SIZE(A, 1)) :: lambda
        LOGICAL :: is_spd_matrix

        IF (PRESENT(is_sym)) THEN
            IF (.NOT. is_sym) THEN
                is_spd_matrix = .FALSE.
                RETURN
            END IF
        ELSE IF (.NOT. is_symmetric(A)) THEN
            is_spd_matrix = .FALSE.
            RETURN
        END IF

        PRINT *, "Checking if the matrix is SPD..."

        CALL Eigen(A, lambda, method = "Power_iteration")
        IF(MINVAL(lambda) < 0) THEN
            is_spd_matrix = .FALSE.
        ELSE
            is_spd_matrix = .TRUE.
        END IF 

    END FUNCTION is_SPD

    FUNCTION is_tridiagonal(A, is_square) RESULT(is_tridiag)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL, OPTIONAL, INTENT(IN) :: is_square
        LOGICAL :: is_tridiag
        INTEGER :: i, j, N

        N = SIZE(A, 1)
        is_tridiag = .TRUE.

        IF (PRESENT(is_square)) THEN
            IF (.NOT. is_square) THEN
                is_tridiag = .FALSE.
                RETURN
            END IF
        ELSE IF (.NOT. is_square_matrix(A)) THEN
            is_tridiag = .FALSE.
            RETURN
        END IF

        PRINT *, "Checking if the matrix is tridiagonal..."

        DO i = 1, N
            DO j = 1, N
                IF (ABS(i-j) > 1) THEN
                    IF (ABS(A(i,j)) > epsi) THEN
                        is_tridiag = .FALSE.
                        RETURN
                    END IF
                END IF
            END DO
        END DO

    END FUNCTION is_tridiagonal

    FUNCTION is_diagonally_dominant(A, is_square) RESULT(is_diag_dom)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL, OPTIONAL, INTENT(IN) :: is_square
        LOGICAL :: is_diag_dom
        INTEGER :: i, N
        REAL(dp) :: row_sum

        print*, "ici"

        N = SIZE(A, 1)
        is_diag_dom = .TRUE.

        IF (PRESENT(is_square)) THEN
            IF (.NOT. is_square) THEN
                is_diag_dom = .FALSE.
                RETURN
            END IF
        ELSE IF (.NOT. is_square_matrix(A)) THEN
            is_diag_dom = .FALSE.
            RETURN
        END IF

        PRINT *, "Checking if the matrix is diagonally dominant..."

        DO i = 1, N
            row_sum = SUM(ABS(A(i, :))) - ABS(A(i, i))
            IF (ABS(A(i, i)) < row_sum) THEN
                is_diag_dom = .FALSE.
                RETURN
            END IF
        END DO

    END FUNCTION is_diagonally_dominant

    FUNCTION is_non_zero_diagonal(A, is_square) RESULT(is_non_zero_diag)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL, OPTIONAL, INTENT(IN) :: is_square
        LOGICAL :: is_non_zero_diag

        is_non_zero_diag = .TRUE.

        IF (PRESENT(is_square)) THEN
            IF (.NOT. is_square) THEN
                is_non_zero_diag = .FALSE.
                RETURN
            END IF
        ELSE IF (.NOT. is_square_matrix(A)) THEN
            is_non_zero_diag = .FALSE.
            RETURN
        END IF

        PRINT *, "Checking if the matrix has non-zero diagonal elements..."

        IF (ANY(ABS(Diag(A)) < epsi)) is_non_zero_diag = .FALSE.

    END FUNCTION is_non_zero_diagonal

END MODULE NAFPack_matrix_properties