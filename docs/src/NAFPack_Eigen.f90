!> Module for eigenvalue and eigenvector computations in NAFPack
module NAFPack_Eigen

    use NAFPack_kinds, only: dp
    use NAFPack_constant, only: TOL_CONVERGENCE_dp, MAX_ITERATION
    use NAFPack_matrix_decomposition, only: QR_decomposition
    use NAFPack_matricielle, only: Identity_n, normalise

    implicit none(type, external)

    private

    public :: Eigen

contains

    !================== Eigen ===============================================================

    !> Computes the eigenvalues and eigenvectors of a matrix A
    !> \[ A * \vec{v} = \lambda * \vec{v} \]
    !> with **A** a square matrix, **λ** the eigenvalue, and **v** the eigenvector.
    !> This subroutine allows you to choose the method for computing eigenvalues and eigenvectors:
    !>
    !> - Power iteration
    !> - QR algorithm (with or without shift)
    !> The default method is Power iteration.
    subroutine Eigen(A, lambda, vp, method, k)
        real(dp), dimension(:, :), intent(in) :: A
        character(LEN=*), optional, intent(in) :: method
        integer, optional, intent(in) :: k
        real(dp), dimension(:, :), optional, intent(out) :: vp
        real(dp), dimension(:), intent(out) :: lambda
        real(dp), dimension(size(A, 1), size(A, 1)) :: A_tmp
        real(dp), dimension(size(A, 1), size(A, 1)) :: vp_tmp
        character(LEN=50) :: base_method
        integer :: N, i, k_max, pos

        if (present(k)) then
            if (k <= 0) stop "ERROR :: k must be a positive integer"
            k_max = k
        else
            k_max = MAX_ITERATION
        end if

        N = size(A, 1)
        if (size(A, 2) /= N) stop "ERROR :: Matrix A not square"

        if (size(lambda, 1) /= N) stop "ERROR :: dimension lambda"
        if (present(vp) .and. (size(vp, 1) /= N .or. size(vp, 2) /= N)) stop "ERROR :: dimension vp"

        if (method == "Power_iteration") then

            A_tmp = A
            do i = 1, N
                call Power_iteration(A_tmp, lambda(i), vp_tmp(i, :), k_max)
                A_tmp = deflation(A_tmp, lambda(i), vp_tmp(i, :), k_max)
            end do

            if (present(vp)) vp = vp_tmp

        else if (index(method, "QR") == 1) then

            if (present(vp)) vp = 0
            if (present(vp)) print*,"WARNING :: No solution for eigenvectors with the QR method"

            pos = index(trim(method), "_Shifted")

            if (pos > 0 .and. pos + 7 == len_trim(method)) then
                base_method = method(:pos - 1)
                call Eigen_QR_Shifted(A, lambda, base_method, N, k_max)
            else
                call Eigen_QR(A, lambda, method, N, k_max)
            end if

        else
            stop "ERROR :: Wrong method for Eigen"
        end if

    end subroutine Eigen

    !> QR algorithm for computing eigenvalues
    !>
    !> This subroutine implements the QR algorithm for computing the eigenvalues of a matrix.
    subroutine Eigen_QR(A, lambda, method, N, k)
        real(dp), dimension(:, :), intent(in) :: A
        character(LEN=*), intent(in) :: method
        integer, intent(in) :: N, k
        real(dp), dimension(:), intent(out) :: lambda
        real(dp), dimension(size(A, 1)) :: lambda_old
        real(dp), dimension(size(A, 1), size(A, 1)) :: A_tmp, Q, R
        real(dp) :: diff
        integer :: i, j

        A_tmp = A

        do i = 1, k

            lambda_old = lambda

            call QR_decomposition(A_tmp, method, Q, R)

            A_tmp = matmul(R, Q)

            diff = abs(A_tmp(2, 1))
            do j = 3, N
                if (maxval(abs(A_tmp(j, 1:j - 1))) > diff) then
                    diff = maxval(abs(A_tmp(j, 1:j - 1)))
                end if
            end do

            if (i == k) then
                print*," WARNING :: non-convergence of the QR Algorithm for eigenvalues "//method
                print*,"convergence = ", diff
                exit
            end if

            if (diff <= TOL_CONVERGENCE_dp) exit
        end do

        ! Extract eigenvalues
        lambda = [(A_tmp(i, i), i=1, N)]

    end subroutine Eigen_QR

    !> Shifted QR algorithm for computing eigenvalues
    !>
    !> This subroutine implements the shifted QR algorithm for computing the eigenvalues of a matrix.
    !> The shift is chosen as the last diagonal element of the matrix.
    subroutine Eigen_QR_Shifted(A, lambda, method, N, k)
        integer, intent(in) :: N, k
        character(LEN=*), intent(in) :: method
        real(dp), dimension(N, N), intent(in) :: A
        real(dp), dimension(N), intent(out) :: lambda
        integer :: i, j
        real(dp), dimension(size(A, 1), size(A, 1)) :: A_tmp, Q, R, Id
        real(dp) :: shift, diff

        A_tmp = A
        Id = Identity_n(N)

        do i = 1, k
            !choice of shift: last diagonal element
            shift = A_tmp(N, N)

            ! Gap : A - µI
            A_tmp = A_tmp - shift * Id

            ! QR Decomposition : A - µI = Q * R
            call QR_decomposition(A_tmp, method, Q, R)

            ! A = RQ + µI
            A_tmp = matmul(R, Q) + shift * Id

            diff = abs(A_tmp(2, 1))
            do j = 3, N
                if (maxval(abs(A_tmp(j, 1:j - 1))) > diff) then
                    diff = maxval(abs(A_tmp(j, 1:j - 1)))
                end if
            end do

            if (i == k) then
                print*,"WARNING :: non-convergence of the Shifted QR Algorithm for eigenvalues ", &
                    trim(method)
                print*,"convergence = ", diff
                exit
            end if

            if (diff <= TOL_CONVERGENCE_dp) exit

        end do

        ! Extract eigenvalues
        lambda = [(A_tmp(i, i), i=1, N)]

    end subroutine Eigen_QR_Shifted

    !> Power iteration method for computing the dominant eigenvalue and eigenvector
    !>
    !> This subroutine implements the power iteration method for finding the dominant eigenvalue and eigenvector of a matrix.
    !> It iteratively computes the eigenvector and eigenvalue until convergence
    subroutine Power_iteration(A, lambda, vp, k)
        real(dp), dimension(:, :), intent(in) :: A
        integer, intent(in) :: k
        real(dp), dimension(:), intent(out) :: vp
        real(dp), intent(out) :: lambda
        real(dp), dimension(size(A, 1)) :: u, vp_tmp, r
        integer :: i, N

        N = size(A, 1)
        call random_number(u)

        u = normalise(u)
        vp_tmp = matmul(A, u)
        lambda = dot_product(vp_tmp, u)
        r = vp_tmp - lambda * u

        do i = 1, k
            u = normalise(vp_tmp)
            vp_tmp = matmul(A, u)
            lambda = dot_product(vp_tmp, u)
            if (norm2(r) <= TOL_CONVERGENCE_dp) exit
            r = vp_tmp - lambda * u
            if (i == k) then
                print*,"WARNING :: non-convergence of the power iteration method"
            end if
        end do

        vp = u

    end subroutine Power_iteration

    !> Deflation method for removing the influence of an eigenvalue and eigenvector
    !>
    !> This function performs deflation on a matrix A by removing the influence of an eigenvalue and its corresponding eigenvector.
    function deflation(A, lambda, vp, k) result(result)

        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(in) :: vp
        real(dp), intent(in) :: lambda
        integer, intent(in) :: k
        real(dp), dimension(size(A, 1), size(A, 1)) :: result
        real(dp), dimension(size(A, 1)) :: wp
        integer :: i, j, N
        real(dp) :: lambda1

        N = size(A, 1)
        result = A

        call Power_iteration(transpose(A), lambda1, wp, k)
        do i = 1, N
            do j = 1, N
                result(i, j) = result(i, j) - (lambda * vp(i) * wp(j)) / dot_product(vp, wp)
            end do
        end do

    end function deflation

end module NAFPack_Eigen
