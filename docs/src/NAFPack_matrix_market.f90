MODULE NAFPack_matrix_market

    USE NAFPack_constant

    IMPLICIT NONE(TYPE, EXTERNAL)

CONTAINS

    SUBROUTINE readMatrixMarket(filename, A)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        REAL(dp), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: A

        INTEGER :: i, j, nrows, ncols, nnz, ios, k, row_idx, col_idx
        CHARACTER(LEN=256) :: line, header
        CHARACTER(LEN=50) :: object, format_type, field, symmetry
        INTEGER :: unit
        REAL(dp) :: val

        unit = 10

        ! Open the file
        OPEN (NEWUNIT=unit, FILE=filename, STATUS='old', ACTION='read', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT*,'Error: opening file: ', filename
            STOP
        END IF

        ! read the header
        READ (unit, '(A)', IOSTAT=ios) header
        IF (ios /= 0) THEN
            PRINT*,'Error: reading header Matrix Market'
            CLOSE (unit)
            STOP
        END IF

        ! Check Matrix Market format
        IF (index(header, '%%MatrixMarket') == 0) THEN
            PRINT*,'Error: file is not in Matrix Market format'
            CLOSE (unit)
            RETURN
        END IF

        ! Parse the header to extract information
        READ (header, *, iostat=ios) object, object, format_type, field, symmetry
        IF (ios /= 0) THEN
            PRINT*,'Error reading Matrix Market header'
            CLOSE (unit)
            STOP
        END IF

        ! Skip comment lines and find the size line
        DO
            READ (unit, '(A)', IOSTAT=ios) line
            IF (ios /= 0) THEN
                PRINT*,'Error: malformed file'
                CLOSE (unit)
                STOP
            END IF
            IF (line(1:1) /= '%') THEN
                EXIT
            END IF
        END DO

        ! Read matrix dimensions
        READ (line, *, iostat=ios) nrows, ncols, nnz
        IF (ios /= 0) THEN
            PRINT*,'Error: unable to read dimensions'
            CLOSE (unit)
            RETURN
        END IF

        ALLOCATE (A(nrows, ncols))
        A = 0.0_dp

        ! Read matrix entries
        IF (trim(format_type) == 'coordinate') THEN
            ! Coordinate format (sparse)
            DO k = 1, nnz
                READ (unit, *, iostat=ios) row_idx, col_idx, val
                IF (ios /= 0) THEN
                    PRINT*,'Error: reading entry', k
                    EXIT
                END IF

                A(row_idx, col_idx) = val

                ! If the matrix is symmetric, also fill A(j,i)
                IF (trim(symmetry) == 'symmetric' .AND. row_idx /= col_idx) THEN
                    A(col_idx, row_idx) = val
                END IF
            END DO

        ELSE IF (trim(format_type) == 'array') THEN
            ! Dense format (array)
            DO j = 1, ncols
                DO i = 1, nrows
                    READ (unit, *, iostat=ios) val
                    IF (ios /= 0) THEN
                        PRINT*,'Error: reading element (', i, ',', j, ')'
                        CLOSE (unit)
                        RETURN
                    END IF
                    A(i, j) = val
                END DO
            END DO

        ELSE
            PRINT*,'Error: unsupported format: ', trim(format_type)
            CLOSE (unit)
            RETURN
        END IF

        CLOSE (unit)
    END SUBROUTINE readMatrixMarket

END MODULE NAFPack_matrix_market
