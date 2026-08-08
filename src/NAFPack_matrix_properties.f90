module NAFPack_matrix_properties

    use NAFPack_kinds, only: dp
    use NAFPack_constant, only: TOL_TEST_dp
    use NAFPack_matricielle, only: Identity_n, Diag
    use NAFPack_Eigen, only: Eigen

    implicit none(type, external)

    private
    public :: is_square_matrix, is_symmetric, is_orthogonal, is_SPD, is_tridiagonal, &
              is_diagonally_dominant, is_non_zero_diagonal

contains

    function is_square_matrix(A) result(is_square)
        real(dp), dimension(:, :), intent(in) :: A
        logical :: is_square

        is_square = (size(A, 1) == size(A, 2))

    end function is_square_matrix

    function is_symmetric(A) result(is_sym)
        real(dp), dimension(:, :), intent(in) :: A
        logical :: is_sym

        is_sym = all(A == transpose(A))

    end function is_symmetric

    function is_orthogonal(A) result(is_orth)
        real(dp), dimension(:, :), intent(in) :: A
        logical :: is_orth

        is_orth = all(abs(matmul(A, transpose(A)) - Identity_n(size(A, 1))) < TOL_TEST_dp)

    end function is_orthogonal

    function is_SPD(A, is_sym) result(is_spd_matrix)
        real(dp), dimension(:, :), intent(in) :: A
        logical, optional, intent(in) :: is_sym
        real(dp), dimension(size(A, 1)) :: lambda
        logical :: is_spd_matrix

        if (present(is_sym)) then
            if (.not. is_sym) then
                is_spd_matrix = .false.
                return
            end if
        else if (.not. is_symmetric(A)) then
            is_spd_matrix = .false.
            return
        end if

        call Eigen(A, lambda, method="Power_iteration")
        if (minval(lambda) < 0) then
            is_spd_matrix = .false.
        else
            is_spd_matrix = .true.
        end if

    end function is_SPD

    function is_tridiagonal(A) result(is_tridiag)
        real(dp), dimension(:, :), intent(in) :: A
        logical :: is_tridiag
        integer :: i, j, N

        N = size(A, 1)
        is_tridiag = .true.

        do i = 1, N
            do j = 1, N
                if (abs(i - j) > 1) then
                    if (abs(A(i, j)) > TOL_TEST_dp) then
                        is_tridiag = .false.
                        return
                    end if
                end if
            end do
        end do

    end function is_tridiagonal

    function is_diagonally_dominant(A) result(is_diag_dom)
        real(dp), dimension(:, :), intent(in) :: A
        logical :: is_diag_dom
        integer :: i, N
        real(dp) :: row_sum

        N = size(A, 1)
        is_diag_dom = .true.

        do i = 1, N
            row_sum = sum(abs(A(i, :))) - abs(A(i, i))
            if (abs(A(i, i)) < row_sum) then
                is_diag_dom = .false.
                return
            end if
        end do

    end function is_diagonally_dominant

    function is_non_zero_diagonal(A) result(is_non_zero_diag)
        real(dp), dimension(:, :), intent(in) :: A
        logical :: is_non_zero_diag

        is_non_zero_diag = .true.

        if (any(abs(Diag(A)) < TOL_TEST_dp)) is_non_zero_diag = .false.

    end function is_non_zero_diagonal

end module NAFPack_matrix_properties
