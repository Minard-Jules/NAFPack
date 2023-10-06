MODULE NAFPack_matricielle

    USE NAFPack_constant
    
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: dot, cross
    PUBLIC :: norm_2, norm_2_complex
    PUBLIC :: normalise, normalise_complexe
    PUBLIC :: Diagonally_Dominant_Matrix
    PUBLIC :: Identity_n

    CONTAINS

    !function that calculates the dot product between a and b
    FUNCTION dot(a,b) RESULT(result)

        REAL(dp), DIMENSION(3) :: a, b
        REAL(dp) :: result

        result = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)

    END FUNCTION dot

    !function that calculates the cross product between a and b
    FUNCTION cross(a,b) RESULT(result)

        REAL(dp), DIMENSION(3) :: a, b
        REAL(dp), DIMENSION(3) :: result

        result(1) = a(2) * b(3) - b(2) * a(3)
        result(2) = - (a(1) * b(3) - b(1) * a(3))
        result(3) = a(1) * b(2) - b(1) * a(2)

    END FUNCTION cross

    !function that calculates the squared norm 2 of a
    FUNCTION norm_2(a) RESULT(result)

        REAL(dp), DIMENSION(:) :: a
        REAL(dp) :: result

        result = SQRT(SUM(a**2))
    
    END FUNCTION norm_2

    !function that calculates the squared norm 2 of a
    FUNCTION norm_2_complex(a) RESULT(result)

        COMPLEX(dp), DIMENSION(:) :: a
        COMPLEX(dp) :: result

        result = SQRT(SUM(a**2))
    
    END FUNCTION norm_2_complex

    !function that transforms a into a unit vector
    FUNCTION normalise(a) RESULT(result)

        REAL(dp), DIMENSION(:) :: a
        REAL(dp), DIMENSION(SIZE(a)) :: result

        result = a / NORM2(a)
    
    END FUNCTION normalise

    !function that transforms a into a unit vector
    FUNCTION normalise_complexe(a) RESULT(result)
        COMPLEX(dp), DIMENSION(:) :: a
        COMPLEX(dp), DIMENSION(SIZE(a)) :: result

        result = a / norm_2_complex(a)
    
    END FUNCTION normalise_complexe

    !function which checks if A is diagonally dominant
    FUNCTION Diagonally_Dominant_Matrix(A) RESULT(diagonally_dominant)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL :: diagonally_dominant
        REAL(dp) :: summation
        INTEGER :: i, N

        N = SIZE(A, 1)

        diagonally_dominant = .TRUE.

        DO i = 1, N
            summation = SUM(ABS(A(i, :) - A(i, i)))
            if (ABS(A(i, i)) < summation)THEN
                diagonally_dominant = .FALSE.
                EXIT
            END IF
        END DO

    END FUNCTION Diagonally_Dominant_Matrix

    !function that returns the identity matrix
    FUNCTION Identity_n(N) RESULT(Identity)
        INTEGER, INTENT(IN) :: N
        REAL(dp), DIMENSION(N, N) :: Identity
        INTEGER :: i

        Identity = 0.d0

        DO i = 1, N
            Identity(i, i) = 1.d0
        END DO

    END FUNCTION Identity_n

END MODULE NAFPack_matricielle