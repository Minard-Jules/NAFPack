!> Module for Tensor operations in NAFPack
MODULE NAFPack_matricielle

    USE NAFPack_constant
    
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: dot, cross
    PUBLIC :: norm_2, norm_2_complex
    PUBLIC :: normalise, normalise_complexe
    PUBLIC :: Diagonally_Dominant_Matrix
    PUBLIC :: Identity_n
    PUBLIC :: rotation_matrix

    CONTAINS

    !> function that calculates the dot product of two real 3-dimensional vectors a and b
    FUNCTION dot(a,b) RESULT(result)

        REAL(dp), DIMENSION(3) :: a, b
        REAL(dp) :: result

        result = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)

    END FUNCTION dot

    !> function that calculates the cross product between two real 3-dimensional vectors a and b
    FUNCTION cross(a,b) RESULT(result)

        REAL(dp), DIMENSION(3) :: a, b
        REAL(dp), DIMENSION(3) :: result

        result(1) = a(2) * b(3) - b(2) * a(3)
        result(2) = - (a(1) * b(3) - b(1) * a(3))
        result(3) = a(1) * b(2) - b(1) * a(2)

    END FUNCTION cross

    !> function that calculates the Euclidean norm (L2 norm) of a real vector a
    FUNCTION norm_2(a) RESULT(result)

        REAL(dp), DIMENSION(:) :: a
        REAL(dp) :: result

        result = SQRT(SUM(a**2))
    
    END FUNCTION norm_2

    !> function that calculates the Euclidean norm (L2 norm) of a complex vector a
    FUNCTION norm_2_complex(a) RESULT(result)

        COMPLEX(dp), DIMENSION(:) :: a
        COMPLEX(dp) :: result

        result = SQRT(SUM(a**2))
    
    END FUNCTION norm_2_complex

    !> function that normalises a real vector a to make it a unit vector
    FUNCTION normalise(a) RESULT(result)

        REAL(dp), DIMENSION(:) :: a
        REAL(dp), DIMENSION(SIZE(a)) :: result

        result = a / norm_2(a)

    END FUNCTION normalise

    !> function that normalises a complex vector a to make it a unit vector
    FUNCTION normalise_complexe(a) RESULT(result)
        COMPLEX(dp), DIMENSION(:) :: a
        COMPLEX(dp), DIMENSION(SIZE(a)) :: result

        result = a / norm_2_complex(a)
    
    END FUNCTION normalise_complexe

    !> function which checks if A is diagonally dominant
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

    !> function that returns the identity matrix for a given size N
    FUNCTION Identity_n(N, use_concurrent) RESULT(Identity)
        INTEGER, INTENT(IN) :: N
        LOGICAL, INTENT(IN), OPTIONAL :: use_concurrent
        REAL(dp), DIMENSION(N, N) :: Identity
        INTEGER :: i
        LOGICAL :: concurrent_mode

        concurrent_mode = .FALSE.
        IF (PRESENT(use_concurrent)) concurrent_mode = use_concurrent

        Identity = 0.d0

        IF (concurrent_mode) THEN
            DO CONCURRENT (i = 1:N)
                Identity(i, i) = 1.0_dp
            END DO
        ELSE
            FORALL(i = 1:N) Identity(i, i) = 1.0_dp
        END IF
        
    END FUNCTION Identity_n

    !> Function to create a rotation matrix 
    !>
    !> This function generates a rotation matrix G based on the input matrix A and the specified rotation indices.
    FUNCTION rotation_matrix(A,rotation) RESULT(G)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        INTEGER, DIMENSION(2), INTENT(IN) :: rotation
        REAL(dp), DIMENSION(SIZE(A, 1),SIZE(A, 2)) :: G
        REAL(dp) :: frac, val_1, val_2
        INTEGER :: i, j

        i = rotation(1)
        j = rotation(2)

        G = Identity_n(SIZE(A, 1))

        val_1 = A(j, j)
        val_2 = A(i, j)

        frac = SQRT(val_1**2 + val_2**2)

        G(i,i) = val_1 / frac
        G(j,j) = val_1 / frac
        G(i,j) = - val_2 / frac
        G(j,i) = val_2 / frac

    END FUNCTION rotation_matrix

END MODULE NAFPack_matricielle