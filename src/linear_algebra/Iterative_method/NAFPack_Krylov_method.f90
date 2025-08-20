module NAFPack_Krylov_method

    use NAFPack_constant
    use NAFPack_matricielle

    implicit none(type, external)

    private

contains

    subroutine lanczos(A, q1, m, Q, T)
        real(dp), dimension(:, :), intent(IN) :: A
        real(dp), dimension(:), intent(IN) :: q1
        integer, intent(IN) :: m
        real(dp), dimension(:, :), intent(OUT) :: Q
        real(dp), dimension(:, :), intent(OUT) :: T
        real(dp), dimension(size(A, 1)) :: y, z
        real(dp), dimension(m) :: alpha
        real(dp), dimension(m - 1) :: beta
        integer :: N, k

        N = size(A, 1)

        alpha = 0.d0
        beta = 0.d0
        Q = 0.d0

        Q(:, 1) = q1 / norm2(q1)

        do k = 1, m
            if (k == 1) then
                y = matmul(A, Q(:, k))
            else
                y = matmul(A, Q(:, k)) - beta(k - 1) * Q(:, k - 1)
            end if

            alpha(k) = dot_product(Q(:, k), y)
            z = y - alpha(k) * Q(:, k)

            if (k < m) then
                beta(k) = norm2(z)
                if (beta(k) < 1.0d-12) exit
                Q(:, k + 1) = z / beta(k)
            end if
        end do

        T = Make_Tridiagonal(beta, alpha, beta)

    end subroutine lanczos

    subroutine Arnoldi(A, q1, m, Q, H)
        real(dp), dimension(:, :), intent(IN) :: A
        real(dp), dimension(:), intent(IN) :: q1
        integer, intent(IN) :: m
        real(dp), dimension(:, :), intent(OUT) :: Q
        real(dp), dimension(:, :), intent(OUT) :: H
        integer :: k, j, N

        N = size(A, 1)
        Q(:, 1) = q1 / norm2(q1)
        H = 0.0d0

        do k = 2, m
            Q(:, k) = matmul(A, Q(:, k - 1))
            do j = 1, k - 1
                H(j, k - 1) = dot_product(Q(:, j), Q(:, k))
                Q(:, k) = Q(:, k) - H(j, k - 1) * Q(:, j)
            end do
            H(k, k - 1) = norm2(Q(:, k))
            Q(:, k) = Q(:, k) / H(k, k - 1)
        end do

    end subroutine Arnoldi

    subroutine Arnoldi_MGS(A, q1, m, Q, H)
        real(dp), dimension(:, :), intent(IN) :: A
        real(dp), dimension(:), intent(IN) :: q1
        integer, intent(IN) :: m
        real(dp), dimension(:, :), intent(OUT) :: Q ! (n, m+1)
        real(dp), dimension(:, :), intent(OUT) :: H ! (m+1, m)
        integer :: k, j, N
        real(dp), dimension(size(A, 1)) :: w

        N = size(A, 1)
        Q(:, 1) = q1 / norm2(q1)
        H = 0.0d0

        do k = 1, m
            ! w = A * q_k
            w = matmul(A, Q(:, k))

            ! Modified Gram-Schmidt orthonormalization
            do j = 1, k
                H(j, k) = dot_product(Q(:, j), w)
                w = w - H(j, k) * Q(:, j)
            end do

            H(k + 1, k) = norm2(w)
            if (H(k + 1, k) > 0.0d0) then
                Q(:, k + 1) = w / H(k + 1, k)
            else
                Q(:, k + 1) = 0.0d0
            end if

        end do

    end subroutine Arnoldi_MGS

end module NAFPack_Krylov_method
