!> Module for Tensor operations in NAFPack
MODULE NAFPack_matricielle

    USE NAFPack_constant

    IMPLICIT NONE(TYPE, EXTERNAL)

    PRIVATE
    PUBLIC :: dot, cross
    PUBLIC :: norm_2_real, norm_2_complex
    PUBLIC :: normalise, normalise_complexe
    PUBLIC :: Diagonally_Dominant_Matrix
    PUBLIC :: Identity_n
    PUBLIC :: rotation_matrix
    PUBLIC :: Trace
    PUBLIC :: Diag, Make_Tridiagonal

CONTAINS

    !> function that calculates the dot product of two real 3-dimensional vectors \( \vec{a} \) and \( \vec{b} \)
    !> \[ \vec{a} \cdot \vec{b} \]
    FUNCTION dot(a, b) RESULT(RESULT)
        REAL(dp), DIMENSION(:), INTENT(IN) :: a, b
        REAL(dp) :: RESULT
        INTEGER :: i

        IF (size(a) /= size(b)) STOP "Error: Vectors must be of the same size."

        RESULT = 0.0_dp
        DO i = 1, size(a)
            RESULT = RESULT + a(i) * b(i)
        END DO

    END FUNCTION dot

    !> function that calculates the cross product between two real 3-dimensional vectors \( \vec{a} \) and \( \vec{b} \)
    !> \[ \vec{a} \times \vec{b} \][^1]
    !> [^1]: the wedge notation \( \vec{a} \wedge \vec{b} \) can sometimes be used to denote the vector product.
    FUNCTION cross(a, b) RESULT(RESULT)

        REAL(dp), DIMENSION(3) :: a, b
        REAL(dp), DIMENSION(3) :: RESULT

        RESULT(1) = a(2) * b(3) - b(2) * a(3)
        RESULT(2) = -(a(1) * b(3) - b(1) * a(3))
        RESULT(3) = a(1) * b(2) - b(1) * a(2)

    END FUNCTION cross

    !> function that calculates the Euclidean norm (L2 norm) of a vector \( \vec{a} \),
    !> where \( \vec{a} \in \mathbb{R}^n \)
    !> \[ ||\vec{a}||_2 = \sqrt{\sum_{i=1}^{n} a_i^2} \quad \text{ with } \quad \sum_{i=1}^{n} a_i^2 = \vec{a} \cdot \vec{a} \]
    !> where \( n \) is the dimension of the real vector \( \vec{a} \).
    FUNCTION norm_2_real(a) RESULT(RESULT)

        REAL(dp), DIMENSION(:) :: a
        REAL(dp) :: RESULT

        RESULT = sqrt(dot_product(a, a))

    END FUNCTION norm_2_real

    !> Optimized norm calculation avoiding overflow/underflow
    PURE FUNCTION norm_2_safe(a) RESULT(RESULT)
        REAL(dp), DIMENSION(:), INTENT(IN) :: a
        REAL(dp) :: RESULT
        REAL(dp) :: scale, sum_of_squares

        scale = maxval(abs(a))
        IF (scale == 0.0_dp) THEN
            RESULT = 0.0_dp
        ELSE
            sum_of_squares = sum((a / scale)**2)
            RESULT = scale * sqrt(sum_of_squares)
        END IF
    END FUNCTION norm_2_safe

    !> function that calculates the Euclidean norm (L2 norm or modulus) of a vector \( \vec{a} \),
    !> where \( \vec{a} \in \mathbb{C}^n \)
    !> \[ ||\vec{a}||_2 = \sqrt{\sum_{i=1}^{n} |a_i|^2} \quad \text{ with } \quad \sum_{i=1}^{n} |a_i|^2 = \vec{a} \cdot \overline{\vec{a}} \]
    !> where \( n \) is the dimension of the complex vector \( \vec{a} \).
    FUNCTION norm_2_complex(a) RESULT(RESULT)

        COMPLEX(dp), DIMENSION(:) :: a
        REAL(dp) :: RESULT

        RESULT = sqrt(REAL(dot_product(a, conjg(a))))

    END FUNCTION norm_2_complex

    !> function that normalises a real vector a to make it a unit vector,
    !> where \( \vec{a} \in \mathbb{R}^n \)
    !> \[ \hat{a} = \frac{\vec{a}}{||\vec{a}||_2} \]
    FUNCTION normalise(a) RESULT(RESULT)

        REAL(dp), DIMENSION(:) :: a
        REAL(dp), DIMENSION(size(a)) :: RESULT

        RESULT = a / norm_2_real(a)

    END FUNCTION normalise

    !> function that normalises a complex vector a to make it a unit vector,
    !> where \( \vec{a} \in \mathbb{C}^n \)
    !> \[ \hat{a} = \frac{\vec{a}}{||\vec{a}||_2} \]
    FUNCTION normalise_complexe(a) RESULT(RESULT)
        COMPLEX(dp), DIMENSION(:) :: a
        COMPLEX(dp), DIMENSION(size(a)) :: RESULT

        RESULT = a / norm_2_complex(a)

    END FUNCTION normalise_complexe

    !> function that calculates the trace of a square matrix \( A \)
    !> \[ \text{Tr}(A) = \sum_{i=1}^{n} A(i,i) \]
    FUNCTION Trace(A) RESULT(RESULT)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp) :: RESULT
        INTEGER :: i, N

        N = size(A, 1)
        IF (size(A, 2) /= N) STOP "Error: Matrix must be square."

        RESULT = sum([(A(i, i), i=1, N)])

    END FUNCTION Trace

    !> function which checks if **A** is diagonally dominant
    !> \[ \forall i, |A(i,i)| \geq \sum_{j \neq i} |A(i,j)| \]
    FUNCTION Diagonally_Dominant_Matrix(A) RESULT(diagonally_dominant)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        LOGICAL :: diagonally_dominant
        REAL(dp) :: summation
        INTEGER :: i, N

        N = size(A, 1)

        diagonally_dominant = .TRUE.

        DO i = 1, N
            summation = sum(abs(A(i, :) - A(i, i)))
            IF (abs(A(i, i)) < summation) THEN
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
        IF (present(use_concurrent)) concurrent_mode = use_concurrent

        Identity = 0.d0

        IF (concurrent_mode) THEN
            DO CONCURRENT(i=1:N)
                Identity(i, i) = 1.0_dp
            END DO
        ELSE
            FORALL (i=1:N) Identity(i, i) = 1.0_dp
        END IF

    END FUNCTION Identity_n

    !> function that extracts the diagonal of a matrix
    !> \[ D = \begin{pmatrix} A(1,1) & 0 & \cdots & 0 \\ 0 & A(2,2) & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & A(n,n) \end{pmatrix} \]
    !> where \( D \) is a vector containing the diagonal elements of the matrix \( A \).
    FUNCTION Diag(A) RESULT(D)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(size(A, 1)) :: D
        INTEGER :: i, N

        N = size(A, 1)

        FORALL (i=1:N) D(i) = A(i, i)

    END FUNCTION Diag

    FUNCTION Make_Tridiagonal(d_minus, d, d_plus) RESULT(T)
        REAL(dp), DIMENSION(:), INTENT(IN) :: d_minus, d, d_plus
        REAL(dp), DIMENSION(size(d, 1), size(d, 1)) :: T
        INTEGER :: i, N

        N = size(d, 1)

        T = 0.d0
        DO i = 1, N
            T(i, i) = d(i)
            IF (i > 1) T(i, i-1) = d_minus(i)
            IF (i < N) T(i, i+1) = d_plus(i)
        END DO

    END FUNCTION Make_Tridiagonal

    !> Function to create a rotation matrix
    !>
    !> This function generates a rotation matrix **G** based on the input matrix **A** and the specified rotation indices.
    FUNCTION rotation_matrix(A, rotation) RESULT(G)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        INTEGER, DIMENSION(2), INTENT(IN) :: rotation
        REAL(dp), DIMENSION(size(A, 1), size(A, 2)) :: G
        REAL(dp) :: frac, val_1, val_2
        INTEGER :: i, j

        i = rotation(1)
        j = rotation(2)

        G = Identity_n(size(A, 1))

        val_1 = A(j, j)
        val_2 = A(i, j)

        frac = sqrt(val_1**2 + val_2**2)

        G(i, i) = val_1 / frac
        G(j, j) = val_1 / frac
        G(i, j) = -val_2 / frac
        G(j, i) = val_2 / frac

    END FUNCTION rotation_matrix

END MODULE NAFPack_matricielle
