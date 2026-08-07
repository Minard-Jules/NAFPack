module NAFPack_math_utils

    use NAFPack_kinds, only: sp, isp

    implicit none(type, external)

    private
    public :: sieve_of_eratosthenes
    public :: is_power_of_two, is_power_of_p, power_of_p_exponent

contains

    pure function sieve_of_eratosthenes(N) result(primes)
        integer(isp), intent(in) :: N
        integer(isp), dimension(:), allocatable :: primes
        logical, dimension(:), allocatable :: is_prime
        integer(isp) :: i, j, count_primes, limit, idx

        allocate (is_prime(0:N))
        is_prime = .true.
        is_prime(0:1) = .false.

        limit = int(sqrt(real(N, kind=sp)))

        do i = 2, limit
            if (is_prime(i)) then
                do j = i * i, N, i
                    is_prime(j) = .false.
                end do
            end if
        end do

        count_primes = count(is_prime)
        allocate (primes(count_primes))
        idx = 1
        do i = 2, N
            if (is_prime(i)) then
                primes(idx) = i
                idx = idx + 1
            end if
        end do
        deallocate (is_prime)
    end function sieve_of_eratosthenes

    pure function is_power_of_two(N) result(value)
        integer(isp), intent(in) :: N
        logical :: value

        if (N < 1) then
            value = .false.
        else
            value = (iand(N, N - 1) == 0)
        end if
    end function is_power_of_two

    pure function is_power_of_p(N, p) result(value)
        integer(isp), intent(in) :: N, p
        logical :: value
        integer(isp) :: tmp

        if (N < 1 .or. p < 2) then
            value = .false.
        else
            tmp = N
            do while (mod(tmp, p) == 0)
                tmp = tmp / p
            end do
            value = (tmp == 1)
        end if
    end function is_power_of_p

    pure function power_of_p_exponent(N, p) result(exponent)
        integer(isp), intent(in) :: N, p
        integer(isp) :: exponent, tmp

        if (N < 1 .or. p < 2) then
            exponent = 0
        else
            exponent = 0
            tmp = N
            do while (mod(tmp, p) == 0)
                tmp = tmp / p
                exponent = exponent + 1
            end do
        end if
    end function power_of_p_exponent

end module NAFPack_math_utils
