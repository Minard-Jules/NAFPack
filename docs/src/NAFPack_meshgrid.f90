!> Module for creating a meshgrid from two vectors
!>
!> This module provides a subroutine to create a meshgrid.
MODULE NAFPack_meshgrid

    USE NAFPack_constant

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: meshgrid
  
    CONTAINS

    !> Make N-dimensional meshgrid from two vectors **x_vector** and **y_vector**
    SUBROUTINE meshgrid(x_vector, y_vector, X, Y)
        REAL(dp), DIMENSION(:), INTENT(IN)   :: x_vector, y_vector
        REAL(dp), DIMENSION(SIZE(y_vector), SIZE(x_vector)), INTENT(OUT)  :: X, Y
        INTEGER :: sX, sY, i

        sX = size(x_vector)
        sY = size(y_vector)

        DO i = 1, sY
            X(i, :) = x_vector
        END DO

        DO i = 1, sX
            Y(:, i) = y_vector
        END DO

    END SUBROUTINE meshgrid
  
END MODULE NAFPack_meshgrid