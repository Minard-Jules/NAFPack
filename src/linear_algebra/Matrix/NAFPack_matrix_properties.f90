MODULE NAFPack_matrix_properties

    USE NAFPack_constant
    USE NAFPack_matricielle
    USE NAFPack_Eigen

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE
    PUBLIC :: is_square_matrix, is_symmetric, is_orthogonal, is_SPD, is_tridiagonal, &
              is_diagonally_dominant, is_non_zero_diagonal

CONTAINS

    FUNCTION is_square_matrix(A) RESULT(is_square)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL :: is_square

        is_square = (size(A, 1) == size(A, 2))

    END FUNCTION is_square_matrix

    FUNCTION is_symmetric(A) RESULT(is_sym)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL :: is_sym

        is_sym = all(A == transpose(A))

    END FUNCTION is_symmetric

    FUNCTION is_orthogonal(A) RESULT(is_orth)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL :: is_orth

        is_orth = all(abs(matmul(A, transpose(A)) - Identity_n(size(A, 1))) < epsi_test)

    END FUNCTION is_orthogonal

    FUNCTION is_SPD(A, is_sym) RESULT(is_spd_matrix)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL, OPTIONAL, INTENT(IN) :: is_sym
        REAL(dp), DIMENSION(size(A, 1)) :: lambda
        LOGICAL :: is_spd_matrix

        IF (present(is_sym)) THEN
            IF (.NOT. is_sym) THEN
                is_spd_matrix = .FALSE.
                RETURN
            END IF
        ELSE IF (.NOT. is_symmetric(A)) THEN
            is_spd_matrix = .FALSE.
            RETURN
        END IF

        CALL Eigen(A, lambda, method="Power_iteration")
        IF (minval(lambda) < 0) THEN
            is_spd_matrix = .FALSE.
        ELSE
            is_spd_matrix = .TRUE.
        END IF

    END FUNCTION is_SPD

    FUNCTION is_tridiagonal(A) RESULT(is_tridiag)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL :: is_tridiag
        INTEGER :: i, j, N

        N = size(A, 1)
        is_tridiag = .TRUE.

        DO i = 1, N
            DO j = 1, N
                IF (abs(i - j) > 1) THEN
                    IF (abs(A(i, j)) > epsi) THEN
                        is_tridiag = .FALSE.
                        RETURN
                    END IF
                END IF
            END DO
        END DO

    END FUNCTION is_tridiagonal

    FUNCTION is_diagonally_dominant(A) RESULT(is_diag_dom)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL :: is_diag_dom
        INTEGER :: i, N
        REAL(dp) :: row_sum

        N = size(A, 1)
        is_diag_dom = .TRUE.

        DO i = 1, N
            row_sum = sum(abs(A(i, :))) - abs(A(i, i))
            IF (abs(A(i, i)) < row_sum) THEN
                is_diag_dom = .FALSE.
                RETURN
            END IF
        END DO

    END FUNCTION is_diagonally_dominant

    FUNCTION is_non_zero_diagonal(A) RESULT(is_non_zero_diag)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL :: is_non_zero_diag

        is_non_zero_diag = .TRUE.

        IF (any(abs(Diag(A)) < epsi)) is_non_zero_diag = .FALSE.

    END FUNCTION is_non_zero_diagonal

END MODULE NAFPack_matrix_properties
