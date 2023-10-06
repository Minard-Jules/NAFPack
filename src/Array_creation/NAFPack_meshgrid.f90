MODULE NAFPack_meshgrid

    USE NAFPack_constant

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: meshgrid
  
    CONTAINS

    !create a meshgrid
    SUBROUTINE meshgrid(varx, vary, X, Y)

        REAL(dp), DIMENSION(:), INTENT(IN)   :: varx, vary
        REAL(dp), DIMENSION(SIZE(vary), SIZE(varx)), INTENT(OUT)  :: X, Y
        INTEGER :: sX, sY, i
    
        sX = size(varx) ; sY = size(vary)
    
        DO i = 1, sY
            X(i, :) = varx
        END DO

        DO i = 1, sX
            Y(:, i) = vary
        END DO

    END SUBROUTINE
  
END MODULE NAFPack_meshgrid