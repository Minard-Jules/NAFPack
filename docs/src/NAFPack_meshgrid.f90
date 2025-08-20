!> Module for creating a meshgrid from two vectors
!>
!> This module provides a subroutine to create a meshgrid.
module NAFPack_meshgrid

    use NAFPack_constant

    implicit none(type, external)

    private
    public :: meshgrid

contains

    !> Make N-dimensional meshgrid from two vectors **x_vector** and **y_vector**
    subroutine meshgrid(x_vector, y_vector, X, Y)
        real(dp), dimension(:), intent(IN) :: x_vector, y_vector
        real(dp), dimension(size(y_vector), size(x_vector)), intent(OUT) :: X, Y
        integer :: sX, sY, i

        sX = size(x_vector)
        sY = size(y_vector)

        do i = 1, sY
            X(i, :) = x_vector
        end do

        do i = 1, sX
            Y(:, i) = y_vector
        end do

    end subroutine meshgrid

end module NAFPack_meshgrid
