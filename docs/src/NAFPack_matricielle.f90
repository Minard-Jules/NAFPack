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

    !> function that calculates the dot product of two real 3-dimensional vectors \( \vec{a} \) and \( \vec{b} \)
    !> \[ \vec{a} \cdot \vec{b} \]
    FUNCTION dot(a, b) RESULT(result)
        REAL(dp), DIMENSION(:), INTENT(IN) :: a, b
        REAL(dp) :: result
        INTEGER :: i

        IF (SIZE(a) /= SIZE(b)) STOP "Error: Vectors must be of the same size."

        result = 0.0_dp
        DO i = 1, SIZE(a)
            result = result + a(i) * b(i)
        END DO

    END FUNCTION dot

    !> function that calculates the cross product between two real 3-dimensional vectors \( \vec{a} \) and \( \vec{b} \)
    !> \[ \vec{a} \times \vec{b} \][^1]
    !> [^1]: the wedge notation \( \vec{a} \wedge \vec{b} \) can sometimes be used to denote the vector product.
    FUNCTION cross(a,b) RESULT(result)

        REAL(dp), DIMENSION(3) :: a, b
        REAL(dp), DIMENSION(3) :: result

        result(1) = a(2) * b(3) - b(2) * a(3)
        result(2) = - (a(1) * b(3) - b(1) * a(3))
        result(3) = a(1) * b(2) - b(1) * a(2)

    END FUNCTION cross

    !> function that calculates the Euclidean norm (L2 norm) of a vector \( \vec{a} \),
    !> where \( \vec{a} \in \mathbb{R}^n \)
    !> \[ ||\vec{a}||_2 = \sqrt{\sum_{i=1}^{n} a_i^2} \quad \text{ with } \quad \sum_{i=1}^{n} a_i^2 = \vec{a} \cdot \vec{a} \]
    !> where \( n \) is the dimension of the real vector \( \vec{a} \).
    FUNCTION norm_2(a) RESULT(result)

        REAL(dp), DIMENSION(:) :: a
        REAL(dp) :: result

        result = SQRT(DOT_PRODUCT(a, a))
    
    END FUNCTION norm_2

    !> function that calculates the Euclidean norm (L2 norm or modulus) of a vector \( \vec{a} \),
    !> where \( \vec{a} \in \mathbb{C}^n \)
    !> \[ ||\vec{a}||_2 = \sqrt{\sum_{i=1}^{n} |a_i|^2} \quad \text{ with } \quad \sum_{i=1}^{n} |a_i|^2 = \vec{a} \cdot \overline{\vec{a}} \]
    !> where \( n \) is the dimension of the complex vector \( \vec{a} \).
    FUNCTION norm_2_complex(a) RESULT(result)

        COMPLEX(dp), DIMENSION(:) :: a
        REAL(dp) :: result

        result = SQRT(REAL(DOT_PRODUCT(a, CONJG(a))))
    
    END FUNCTION norm_2_complex

    !> function that normalises a real vector a to make it a unit vector, 
    !> where \( \vec{a} \in \mathbb{R}^n \)
    !> \[ \hat{a} = \frac{\vec{a}}{||\vec{a}||_2} \]
    FUNCTION normalise(a) RESULT(result)

        REAL(dp), DIMENSION(:) :: a
        REAL(dp), DIMENSION(SIZE(a)) :: result

        result = a / norm_2(a)

    END FUNCTION normalise

    !> function that normalises a complex vector a to make it a unit vector,
    !> where \( \vec{a} \in \mathbb{C}^n \)
    !> \[ \hat{a} = \frac{\vec{a}}{||\vec{a}||_2} \]
    FUNCTION normalise_complexe(a) RESULT(result)
        COMPLEX(dp), DIMENSION(:) :: a
        COMPLEX(dp), DIMENSION(SIZE(a)) :: result

        result = a / norm_2_complex(a)

    END FUNCTION normalise_complexe

    !> function which checks if **A** is diagonally dominant
    !> \[ \forall i, |A(i,i)| \geq \sum_{j \neq i} |A(i,j)| \]
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
    !> \[ I_N = \begin{pmatrix} 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & 1 \end{pmatrix} \]
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
    !> This function generates a rotation matrix **G** based on the input matrix **A** and the specified rotation indices.
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