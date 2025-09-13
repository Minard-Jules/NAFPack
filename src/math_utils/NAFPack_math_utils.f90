module NAFPack_math_utils

    use NAFPack_kinds, only: sp, isp

    implicit none(type, external)

    private
    public :: sieve_of_eratosthenes
    public :: is_power_of_two, power_of_p_exponent

contains

    function sieve_of_eratosthenes(N) result(primes)
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

    function is_power_of_two(N) result(value)
        integer(isp), intent(in) :: N
        logical :: value

        if (N < 1) then
            value = .false.
        else
            value = (iand(N, N - 1) == 0)
        end if
    end function is_power_of_two

    function power_of_p_exponent(N, p) result(exponent)
        integer(isp), intent(in) :: N, p
        integer(isp) :: exponent
        integer(isp) :: tmp

        exponent = 0
        tmp = N
        do while (tmp > 1)
            tmp = tmp / p
            exponent = exponent + 1
        end do
    end function power_of_p_exponent

end module NAFPack_math_utils
