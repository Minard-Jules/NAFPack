!> Module for creating a meshgrid from two vectors
!>
!> This module provides a subroutine to create a meshgrid.
MODULE NAFPack_meshgrid

    USE NAFPack_constant

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: meshgrid
  
    CONTAINS

    !> returns a meshgrid from two vectors
    SUBROUTINE meshgrid(varx, vary, X, Y)
        !> vector representing the coordinates of a grid.
        REAL(dp), DIMENSION(:), INTENT(IN)   :: varx, vary
        !> output grid coordinates,
        !> `X` is a matrix where each row is `varx`
        REAL(dp), DIMENSION(SIZE(vary), SIZE(varx)), INTENT(OUT)  :: X
        !> output grid coordinates,
        !> `Y` is a matrix where each column is `vary`
        REAL(dp), DIMENSION(SIZE(vary), SIZE(varx)), INTENT(OUT)  :: Y
        INTEGER :: sX, sY, i
    
        sX = size(varx) 
        sY = size(vary)
    
        DO i = 1, sY
            X(i, :) = varx
        END DO

        DO i = 1, sX
            Y(:, i) = vary
        END DO

    END SUBROUTINE
  
END MODULE NAFPack_meshgrid