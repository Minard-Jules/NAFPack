module NAFPack_matrix_tools

    use NAFPack_constant, only: dp
    use NAFPack_matricielle, only: Identity_n, Trace

    implicit none(type, external)

    private
    public :: Faddeev_Leverrier

contains

    subroutine Faddeev_Leverrier(A, c, Ainv, success, check)
        integer, parameter :: dp = kind(1.0d0)
        real(dp), dimension(:, :), intent(in) :: A
        logical, optional, intent(in) :: check
        real(dp), dimension(:), intent(out) :: c
        real(dp), dimension(size(A, 1), size(A, 1)), optional, intent(out) :: Ainv
        logical, optional, intent(out) :: success
        real(dp), dimension(size(A, 1), size(A, 1)) :: Bk, I, B_Nm1, AB
        logical :: do_check
        integer :: N, k

        N = size(A, 1)
        do_check = .true.

        if (present(check)) do_check = check

        if (do_check) then
            print*,"Checking if the matrix A is square and size of c is correct"
            if (size(A, 2) /= N .or. size(c) < N + 1) then
                print*,"Error : Matrix A must be square and size of c must be at least N+1"
                stop
            end if
        end if

        ! Initialization
        I = Identity_n(N)
        c = 0.0_dp
        c(1) = 1.0_dp
        c(2) = -Trace(A)
        Bk = A + c(2) * I

        do k = 2, N
            AB = matmul(A, Bk)
            c(k + 1) = -Trace(AB) / real(k, dp)
            Bk = AB + c(k + 1) * I
            if (k == N - 1 .and. present(Ainv)) B_Nm1 = -Bk
        end do

        if (present(Ainv) .and. present(success)) then
            if (abs(c(N + 1)) < 1.0e-12_dp) then
                success = .false.
                Ainv = 0.0_dp
            else
                success = .true.
                Ainv = B_Nm1 / c(N + 1)
            end if
        else if (present(Ainv)) then
            if (abs(c(N + 1)) < 1.0e-12_dp) then
                Ainv = 0.0_dp
            else
                Ainv = B_Nm1 / c(N + 1)
            end if
        end if

    end subroutine Faddeev_Leverrier

end module NAFPack_matrix_tools
