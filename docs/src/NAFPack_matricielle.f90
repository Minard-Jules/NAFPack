!> Module for Tensor operations in NAFPack
module NAFPack_matricielle

    use NAFPack_constant

    implicit none(type, external)

    private
    public :: dot, cross
    public :: norm_2_real, norm_2_complex
    public :: normalise, normalise_complexe
    public :: Diagonally_Dominant_Matrix
    public :: Identity_n
    public :: rotation_matrix
    public :: Trace
    public :: Diag, Make_Tridiagonal

contains

    !> function that calculates the dot product of two real 3-dimensional vectors \( \vec{a} \) and \( \vec{b} \)
    !> \[ \vec{a} \cdot \vec{b} \]
    function dot(a, b) result(result)
        real(dp), dimension(:), intent(IN) :: a, b
        real(dp) :: result
        integer :: i

        if (size(a) /= size(b)) stop "Error: Vectors must be of the same size."

        result = 0.0_dp
        do i = 1, size(a)
            result = result + a(i) * b(i)
        end do

    end function dot

    !> function that calculates the cross product between two real 3-dimensional vectors \( \vec{a} \) and \( \vec{b} \)
    !> \[ \vec{a} \times \vec{b} \][^1]
    !> [^1]: the wedge notation \( \vec{a} \wedge \vec{b} \) can sometimes be used to denote the vector product.
    function cross(a, b) result(result)

        real(dp), dimension(3) :: a, b
        real(dp), dimension(3) :: result

        result(1) = a(2) * b(3) - b(2) * a(3)
        result(2) = -(a(1) * b(3) - b(1) * a(3))
        result(3) = a(1) * b(2) - b(1) * a(2)

    end function cross

    !> function that calculates the Euclidean norm (L2 norm) of a vector \( \vec{a} \),
    !> where \( \vec{a} \in \mathbb{R}^n \)
    !> \[ ||\vec{a}||_2 = \sqrt{\sum_{i=1}^{n} a_i^2} \quad \text{ with } \quad \sum_{i=1}^{n} a_i^2 = \vec{a} \cdot \vec{a} \]
    !> where \( n \) is the dimension of the real vector \( \vec{a} \).
    function norm_2_real(a) result(result)

        real(dp), dimension(:) :: a
        real(dp) :: result

        result = sqrt(dot_product(a, a))

    end function norm_2_real

    !> Optimized norm calculation avoiding overflow/underflow
    pure function norm_2_safe(a) result(result)
        real(dp), dimension(:), intent(IN) :: a
        real(dp) :: result
        real(dp) :: scale, sum_of_squares

        scale = maxval(abs(a))
        if (scale == 0.0_dp) then
            result = 0.0_dp
        else
            sum_of_squares = sum((a / scale)**2)
            result = scale * sqrt(sum_of_squares)
        end if
    end function norm_2_safe

    !> function that calculates the Euclidean norm (L2 norm or modulus) of a vector \( \vec{a} \),
    !> where \( \vec{a} \in \mathbb{C}^n \)
    !> \[ ||\vec{a}||_2 = \sqrt{\sum_{i=1}^{n} |a_i|^2} \quad \text{ with } \quad \sum_{i=1}^{n} |a_i|^2 = \vec{a} \cdot \overline{\vec{a}} \]
    !> where \( n \) is the dimension of the complex vector \( \vec{a} \).
    function norm_2_complex(a) result(result)

        complex(dp), dimension(:) :: a
        real(dp) :: result

        result = sqrt(real(dot_product(a, conjg(a))))

    end function norm_2_complex

    !> function that normalises a real vector a to make it a unit vector,
    !> where \( \vec{a} \in \mathbb{R}^n \)
    !> \[ \hat{a} = \frac{\vec{a}}{||\vec{a}||_2} \]
    function normalise(a) result(result)

        real(dp), dimension(:) :: a
        real(dp), dimension(size(a)) :: result

        result = a / norm_2_real(a)

    end function normalise

    !> function that normalises a complex vector a to make it a unit vector,
    !> where \( \vec{a} \in \mathbb{C}^n \)
    !> \[ \hat{a} = \frac{\vec{a}}{||\vec{a}||_2} \]
    function normalise_complexe(a) result(result)
        complex(dp), dimension(:) :: a
        complex(dp), dimension(size(a)) :: result

        result = a / norm_2_complex(a)

    end function normalise_complexe

    !> function that calculates the trace of a square matrix \( A \)
    !> \[ \text{Tr}(A) = \sum_{i=1}^{n} A(i,i) \]
    function Trace(A) result(result)
        real(dp), dimension(:, :), intent(IN) :: A
        real(dp) :: result
        integer :: i, N

        N = size(A, 1)
        if (size(A, 2) /= N) stop "Error: Matrix must be square."

        result = sum([(A(i, i), i=1, N)])

    end function Trace

    !> function which checks if **A** is diagonally dominant
    !> \[ \forall i, |A(i,i)| \geq \sum_{j \neq i} |A(i,j)| \]
    function Diagonally_Dominant_Matrix(A) result(diagonally_dominant)
        real(dp), dimension(:, :), intent(IN) :: A
        logical :: diagonally_dominant
        real(dp) :: summation
        integer :: i, N

        N = size(A, 1)

        diagonally_dominant = .true.

        do i = 1, N
            summation = sum(abs(A(i, :) - A(i, i)))
            if (abs(A(i, i)) < summation) then
                diagonally_dominant = .false.
                exit
            end if
        end do

    end function Diagonally_Dominant_Matrix

    !> function that returns the identity matrix for a given size N
    !> \[ I_N = \begin{pmatrix} 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & 1 \end{pmatrix} \]
    function Identity_n(N, use_concurrent) result(Identity)
        integer, intent(IN) :: N
        logical, intent(IN), optional :: use_concurrent
        real(dp), dimension(N, N) :: Identity
        integer :: i
        logical :: concurrent_mode

        concurrent_mode = .false.
        if (present(use_concurrent)) concurrent_mode = use_concurrent

        Identity = 0.d0

        if (concurrent_mode) then
            do concurrent(i=1:N)
                Identity(i, i) = 1.0_dp
            end do
        else
            forall (i=1:N) Identity(i, i) = 1.0_dp
        end if

    end function Identity_n

    !> function that extracts the diagonal of a matrix
    !> \[ D = \begin{pmatrix} A(1,1) & 0 & \cdots & 0 \\ 0 & A(2,2) & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & A(n,n) \end{pmatrix} \]
    !> where \( D \) is a vector containing the diagonal elements of the matrix \( A \).
    function Diag(A) result(D)
        real(dp), dimension(:, :), intent(IN) :: A
        real(dp), dimension(size(A, 1)) :: D
        integer :: i, N

        N = size(A, 1)

        forall (i=1:N) D(i) = A(i, i)

    end function Diag

    function Make_Tridiagonal(d_minus, d, d_plus) result(T)
        real(dp), dimension(:), intent(IN) :: d_minus, d, d_plus
        real(dp), dimension(size(d, 1), size(d, 1)) :: T
        integer :: i, N

        N = size(d, 1)

        T = 0.d0
        do i = 1, N
            T(i, i) = d(i)
            if (i > 1) T(i, i - 1) = d_minus(i)
            if (i < N) T(i, i + 1) = d_plus(i)
        end do

    end function Make_Tridiagonal

    !> Function to create a rotation matrix
    !>
    !> This function generates a rotation matrix **G** based on the input matrix **A** and the specified rotation indices.
    function rotation_matrix(A, rotation) result(G)

        real(dp), dimension(:, :), intent(IN) :: A
        integer, dimension(2), intent(IN) :: rotation
        real(dp), dimension(size(A, 1), size(A, 2)) :: G
        real(dp) :: frac, val_1, val_2
        integer :: i, j

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

    end function rotation_matrix

end module NAFPack_matricielle
