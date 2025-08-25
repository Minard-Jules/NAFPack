!> Module for matrix decomposition methods
!>
!> This module provides subroutines for various matrix decomposition methods including LU, LDU, Cholesky, and QR decompositions.
module NAFPack_matrix_decomposition

    use NAFPack_constant
    use NAFPack_matricielle

    implicit none(type, external)

    private

    public :: forward, backward
    public :: LU_decomposition, LDU_decomposition, ILU_decomposition
    public :: Cholesky_decomposition, LDL_Cholesky_decomposition, Incomplete_Cholesky_decomposition
    public :: QR_decomposition
    public :: QR_Householder_decomposition, QR_Givens_decomposition, &
              QR_Gram_Schmidt_Classical_decomposition, QR_Gram_Schmidt_Modified_decomposition
    public :: pivot_partial, pivot_total

contains

    !> forward algorithm,
    !> solves the system
    !> \[ L * y = b \]
    !> where **L** is a lower triangular matrix and **b** is a vector
    function forward(L, b) result(y)
        real(dp), dimension(:, :), intent(in) :: L
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(L, 1)) :: y
        integer :: i, N

        N = size(L, 1)

        y(1) = b(1) / L(1, 1)

        do i = 2, N
            y(i) = (b(i) - dot_product(L(i, 1:i - 1), y(1:i - 1))) / L(i, i)
        end do

    end function forward

    !> backward algorithm,
    !> solves the system
    !> \[ U * x = y \]
    !> where **U** is an upper triangular matrix and **y** is a vector
    function backward(U, y) result(x)
        real(dp), dimension(:, :), intent(in) :: U
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(size(U, 1)) :: x
        integer :: i, N

        N = size(U, 1)

        x(N) = y(N) / U(N, N)

        do i = N - 1, 1, -1
            x(i) = (y(i) - dot_product(U(i, i + 1:N), x(i + 1:N))) / U(i, i)
        end do

    end function backward

    !> LU decomposition of a matrix A
    !> \[ A = LU \]
    !> This subroutine performs LU decomposition of a given matrix **A**, where **L** is a lower triangular matrix and **U** is an upper triangular matrix.
    subroutine LU_decomposition(A, L, U)

        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 1)), intent(out) :: L, U
        integer :: i, j, N

        N = size(A, 1)

        L = 0.d0
        U = 0.d0

        do j = 1, N
            L(j, j) = 1.d0

            do i = 1, j
                U(i, j) = A(i, j) - dot_product(L(i, 1:i - 1), U(1:i - 1, j))
            end do

            do i = j + 1, N
                L(i, j) = (A(i, j) - dot_product(L(i, 1:j - 1), U(1:j - 1, j))) / U(j, j)
            end do
        end do

    end subroutine LU_decomposition

    !> LDU decomposition of a matrix A
    !> \[ A = LDU \]
    !> This subroutine performs LDU decomposition of a given matrix **A**, where **L** is a lower triangular matrix, **D** is a diagonal matrix, and **U** is an upper triangular matrix.
    subroutine LDU_decomposition(A, L, D, U)

        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 1)), intent(out) :: L, U, D
        integer :: i, j, k, N

        N = size(A, 1)

        L = 0.d0
        D = 0.d0
        U = 0.d0

        do j = 1, N
            L(j, j) = 1.d0
            U(j, j) = 1.d0

            do i = 1, j - 1
                U(i, j) = (A(i, j) - dot_product(L(i, 1:i - 1), U(1:i - 1, j)*[(D(k, k), k=1, i - 1)])) / D(i, i)
            end do

            i = j
            D(j, j) = A(j, j) - dot_product(L(j, 1:j - 1), U(1:j - 1, j)*[(D(k, k), k=1, j - 1)])

            do i = j + 1, N
                L(i, j) = (A(i, j) - dot_product(L(i, 1:j - 1), U(1:j - 1, j)*[(D(k, k), k=1, j - 1)])) / D(j, j)
            end do
        end do

    end subroutine LDU_decomposition

    !> Incomplete LU decomposition of a matrix A
    !> \[ A \approx LU \]
    !> This subroutine performs incomplete LU decomposition of a given matrix **A**, where **L** is a lower triangular matrix and **U** is an upper triangular matrix.
    subroutine ILU_decomposition(A, L, U, level)

        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 1)), intent(out) :: L, U
        integer, optional, intent(in) :: level
        integer :: N, i, j
        integer, dimension(size(A, 1), size(A, 1)) :: fill_level
        logical, dimension(size(A, 1), size(A, 1)) :: S

        N = size(A, 1)

        L = 0.d0
        U = 0.d0

        if (present(level)) then
            call compute_fill_pattern_ILU(A, fill_level, level, N)
            S = (fill_level <= level)
        else
            S = A /= 0
        end if

        do j = 1, N
            L(j, j) = 1.d0

            do i = 1, j
                if (S(i, j)) U(i, j) = A(i, j) - dot_product(L(i, 1:i - 1), U(1:i - 1, j))
            end do

            if (abs(U(j, j)) < 1.0e-12_dp) then
                print*,"Warning: Near-zero pivot at row ", j, ", value =", U(j, j)
                print*,"Replacing with small value, value =", sign(1.0e-12_dp, U(j, j))
                U(j, j) = sign(1.0e-12_dp, U(j, j))
            end if

            do i = j + 1, N
                if (S(i, j)) L(i, j) = (A(i, j) - dot_product(L(i, 1:j - 1), U(1:j - 1, j))) / U(j, j)
            end do
        end do

    end subroutine ILU_decomposition

    !> Cholesky decomposition of a matrix A
    !> \[ A = LL^T \]
    !> This subroutine performs Cholesky decomposition of a given symmetric positive definite matrix **A**, where **L** is a lower triangular matrix.
    subroutine Cholesky_decomposition(A, L)

        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 1)), intent(out) :: L
        integer :: i, j, N

        N = size(A, 1)

        do j = 1, N
            L(j, j) = sqrt(A(j, j) - dot_product(L(j, 1:j - 1), L(j, 1:j - 1)))

            do i = j + 1, N
                L(i, j) = (A(i, j) - dot_product(L(i, 1:j - 1), L(j, 1:j - 1))) / L(j, j)
            end do
        end do

    end subroutine Cholesky_decomposition

    !> Alternative Cholesky decomposition of a matrix A
    !> \[ A = LDL^T \]
    !> This subroutine performs alternative Cholesky decomposition of a given symmetric positive definite matrix **A**, where **L** is a lower triangular matrix and **D** is a diagonal matrix.
    subroutine LDL_Cholesky_decomposition(A, L, D)

        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 1)), intent(out) :: L, D
        integer :: i, j, N, k

        N = size(A, 1)

        L = Identity_n(N)
        D = 0.d0

        do j = 1, N
            D(j, j) = A(j, j) - dot_product(L(j, 1:j - 1), L(j, 1:j - 1)*[(D(k, k), k=1, j - 1)])

            do i = j + 1, N
                L(i, j) = (A(i, j) - dot_product(L(i, 1:j - 1), L(j, 1:j - 1)*[(D(k, k), k=1, j - 1)])) / D(j, j)
            end do
        end do

    end subroutine LDL_Cholesky_decomposition

    !> Incomplete Cholesky decomposition of a matrix A
    !> \[ A \approx LL^T \]
    !> This subroutine performs incomplete Cholesky decomposition of a given matrix **A**, where **L** is a lower triangular matrix and **U** is an upper triangular matrix.
    subroutine Incomplete_Cholesky_decomposition(A, L, level)

        real(dp), dimension(:, :), intent(in) :: A
        integer, optional, intent(in) :: level
        real(dp), dimension(size(A, 1), size(A, 1)), intent(out) :: L
        logical, dimension(size(A, 1), size(A, 1)) :: S
        integer :: N, i, j
        integer, dimension(size(A, 1), size(A, 1)) :: fill_level

        N = size(A, 1)

        L = 0.d0

        if (present(level)) then
            call compute_fill_pattern_IC(A, fill_level, level, N)
            S = (fill_level <= level)
        else
            S = A /= 0
        end if

        do i = 1, N
            do j = 1, i - 1
                if (S(i, j)) L(i, j) = (A(i, j) - dot_product(L(i, 1:j - 1), L(j, 1:j - 1))) / L(j, j)
            end do

            if (S(i, i)) L(i, i) = sqrt(A(i, i) - dot_product(L(i, 1:i - 1), L(i, 1:i - 1)))
        end do

    end subroutine Incomplete_Cholesky_decomposition

    !> QR decomposition of a matrix **A** using various methods
    !> \[ A = QR \]
    !> This subroutine performs QR decomposition of a given matrix **A** using the specified method (Householder, Givens, Classical Gram-Schmidt, or Modified Gram-Schmidt).
    !> The output matrices **Q** is an orthogonal matrix and **R** is an upper triangular matrix.
    subroutine QR_decomposition(A, method, Q, R)

        real(dp), dimension(:, :), intent(in) :: A
        character(LEN=*), optional, intent(in) :: method
        real(dp), dimension(size(A, 1), size(A, 2)), intent(out) :: Q, R

        if (method == "QR_Householder") then
            call QR_Householder_decomposition(A, Q, R)
        else if (method == "QR_Givens") then
            call QR_Givens_decomposition(A, Q, R)
        else if (method == "QR_Gram_Schmidt_Classical") then
            call QR_Gram_Schmidt_Classical_decomposition(A, Q, R)
        else if (method == "QR_Gram_Schmidt_Modified") then
            call QR_Gram_Schmidt_Modified_decomposition(A, Q, R)
        end if

    end subroutine QR_decomposition

    !> QR decomposition using Householder method
    subroutine QR_Householder_decomposition(A, Q, R)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 2)), intent(out) :: Q, R
        real(dp), dimension(size(A, 1), size(A, 2)) :: Id, H, v_mat_tmp
        real(dp), dimension(size(A, 1)) :: v, u, x
        integer :: N, i, j, k
        real(dp) :: alpha, w, signe, norm_u

        N = size(A, 1)

        R = A

        Id = Identity_n(N)
        Q = Identity_n(N)

        do k = 1, N

            x = 0.d0
            u = 0.d0
            v = 0.d0
            v_mat_tmp = 0.d0
            x(k:N) = R(K:N, K)

            alpha = norm2(R(k:N, k))

            signe = -sign(alpha, x(k))
            u(k:N) = x(k:N) - signe * Id(k:N, k)

            norm_u = norm2(u)
            if (norm_u < epsi) cycle
            v(k:N) = u(k:N) / norm_u

            w = 1.d0
            do i = k, N
                do j = k, N
                    v_mat_tmp(i, j) = v(i) * v(j)
                end do
            end do

            H = Id
            H(k:N, k:N) = Id(k:N, k:N) - (1.d0 + w) * v_mat_tmp(k:N, k:N)

            Q = matmul(Q, H)

            R(k:N, k:N) = matmul(H(k:N, k:N), R(k:N, k:N))

        end do
    end subroutine QR_Householder_decomposition

    !> QR decomposition using Givens rotations
    subroutine QR_Givens_decomposition(A, Q, R)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 2)), intent(out) :: Q, R
        real(dp), dimension(size(A, 1), size(A, 2)) :: G
        integer :: N, i, j

        N = size(A, 1)

        R = A

        Q = Identity_n(N)

        do j = 1, N - 1
            do i = j + 1, N

                G = rotation_matrix(R, [i, j])

                R = matmul(G, R)

                Q = matmul(Q, transpose(G))

            end do
        end do

    end subroutine QR_Givens_decomposition

    !> QR decomposition using Classical Gram-Schmidt method
    subroutine QR_Gram_Schmidt_Classical_decomposition(A, Q, R)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 2)), intent(out) :: Q, R
        real(dp), dimension(size(A, 1)) :: u
        integer :: N, i, j

        N = size(A, 1)
        Q = 0.d0
        R = 0.d0

        do j = 1, N
            u = A(:, j)
            do i = 1, j - 1
                R(i, j) = dot_product(Q(:, i), A(:, j))
                u = u - (R(i, j) * Q(:, i))
            end do
            R(j, j) = norm2(u)
            Q(:, j) = u / R(j, j)
        end do

    end subroutine QR_Gram_Schmidt_Classical_decomposition

    !> QR decomposition using Modified Gram-Schmidt method
    subroutine QR_Gram_Schmidt_Modified_decomposition(A, Q, R)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 2)), intent(out) :: Q, R
        real(dp), dimension(size(A, 1), size(A, 2)) :: u
        integer :: N, i, j

        N = size(A, 1)
        u = A
        Q = 0.d0
        R = 0.d0

        do i = 1, N
            R(i, i) = norm2(u(:, i))
            Q(:, i) = u(:, i) / R(i, i)
            do j = i + 1, N
                R(i, j) = dot_product(Q(:, i), u(:, j))
                u(:, j) = u(:, j) - R(i, j) * Q(:, i)
            end do
        end do

    end subroutine QR_Gram_Schmidt_Modified_decomposition

    subroutine pivot_partial(A, P)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 1)), intent(out) :: P
        integer, dimension(1) :: vlmax
        integer :: N, lmax, k
        real(dp), dimension(size(A, 1), size(A, 1)) :: P_tmp

        N = size(A, 1)
        P = Identity_n(N)

        do k = 1, N - 1

            ! Find the maximum absolute value in the column from row k to N
            vlmax = maxloc(abs(A(k:N, k)))
            lmax = vlmax(1) + k - 1

            !calculate permutation matrix P
            P_tmp = Identity_n(N)
            if (k /= lmax) then
                P_tmp([k, lmax], :) = P_tmp([lmax, k], :)
            end if
            P = matmul(P_tmp, P)

        end do

    end subroutine pivot_partial

    subroutine pivot_total(A, P, Q)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(size(A, 1), size(A, 1)), intent(out) :: P, Q
        real(dp), dimension(size(A, 1), size(A, 1)) :: P_tmp, Q_tmp
        integer, dimension(2) :: vlmax
        integer :: N, lmax, cmax, k

        N = size(A, 1)
        P = Identity_n(N)
        Q = Identity_n(N)

        do k = 1, N - 1
            ! Find max abs element in submatrix
            vlmax = maxloc(abs(A(k:N, k:N)))
            lmax = vlmax(1) + k - 1
            cmax = vlmax(2) + k - 1

            ! permute line if necessary
            P_tmp = Identity_n(N)
            if (k /= lmax) then
                P_tmp([k, lmax], :) = P_tmp([lmax, k], :)
            end if
            P = matmul(P_tmp, P)

            ! permute column if necessary
            Q_tmp = Identity_n(N)
            if (cmax /= k) then
                Q_tmp(:, [k, cmax]) = Q_tmp(:, [cmax, k])
            end if
            Q = matmul(Q_tmp, Q)
        end do

    end subroutine pivot_total

    subroutine compute_fill_pattern_ILU(A, fill_level, max_level, N)
        real(dp), dimension(N, N), intent(in) :: A
        integer, dimension(N, N), intent(out) :: fill_level
        integer, intent(in) :: max_level, N
        logical, dimension(N, N) :: S
        integer :: new_level
        integer :: i, j, k

        ! Niveau initial basé sur A
        fill_level = int_inf
        S = A /= 0.d0
        where (S) fill_level = 0

        ! Calcul symbolique des niveaux de remplissage
        do k = 1, N - 1
            do i = k + 1, N
                if (fill_level(i, k) <= max_level) then
                    do j = k + 1, N
                        if (fill_level(k, j) <= max_level) then
                            new_level = fill_level(i, k) + fill_level(k, j) + 1
                            if (new_level < max_level) then
                                fill_level(i, j) = new_level
                            end if
                        end if
                    end do
                end if
            end do
        end do

    end subroutine compute_fill_pattern_ILU

    subroutine compute_fill_pattern_IC(A, fill_level, max_level, N)
        real(dp), dimension(N, N), intent(in) :: A
        integer, dimension(N, N), intent(out) :: fill_level
        integer, intent(in) :: max_level, N
        logical, dimension(N, N) :: S
        integer :: new_level
        integer :: i, j, k

        ! Niveau initial basé sur A
        fill_level = int_inf
        S = A /= 0.d0
        where (S) fill_level = 0

        ! Calcul symbolique des niveaux de remplissage
        do k = 1, N - 1
            do i = k + 1, N
                if (fill_level(i, k) <= max_level) then
                    do j = k + 1, i
                        if (fill_level(k, j) <= max_level) then
                            new_level = fill_level(i, k) + fill_level(k, j) + 1
                            if (new_level < max_level) then
                                fill_level(i, j) = new_level
                                fill_level(j, i) = new_level
                            end if
                        end if
                    end do
                end if
            end do
        end do

    end subroutine compute_fill_pattern_IC

end module NAFPack_matrix_decomposition
