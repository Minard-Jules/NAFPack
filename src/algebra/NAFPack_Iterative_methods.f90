MODULE NAFPack_Iterative_methods

    USE NAFPack_constant
    USE NAFPack_matrix_decomposition
    USE NAFPack_matricielle

    IMPLICIT NONE
    
    PRIVATE

    PUBLIC :: Jacobi, Gauss_Seidel, SOR, JOR, SIP_ILU, SIP_ICF, SSOR
    
   CONTAINS

    !> Jacobi iterative method
    !>
    !> This subroutine implements the Jacobi method for solving linear systems.
    SUBROUTINE Jacobi(A, b, x0, x)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        REAL(dp), DIMENSION(SIZE(A, 1)), INTENT(OUT) :: x
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x0(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x(i) = x(i) / A(i, i)
        END DO

    END SUBROUTINE Jacobi

    !> Gauss-Seidel iterative method
    !>
    !> This subroutine implements the Gauss-Seidel method for solving linear systems.
    SUBROUTINE Gauss_Seidel(A, b, x0, x)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        REAL(dp), DIMENSION(SIZE(A, 1)), INTENT(OUT) :: x
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x(i) = x(i) / A(i, i)
        END DO

    END SUBROUTINE Gauss_Seidel

    !> Successive Over-Relaxation (SOR) iterative method
    !>
    !> This subroutine implements the SOR method for solving linear systems.
    SUBROUTINE SOR(A, b, x0, x, omega)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        REAL(dp), INTENT(IN) :: omega
        REAL(dp), DIMENSION(SIZE(A, 1)), INTENT(OUT) :: x
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x(i) = omega*(x(i) / A(i, i) - x0(i)) + x0(i)
        END DO

    END SUBROUTINE SOR

    !> Jacobi over-relaxation (JOR) iterative method
    !>
    !> This subroutine implements the Jacobi over-relaxation method for solving linear systems.
    SUBROUTINE JOR(A, b, x0, x, omega)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        REAL(dp), INTENT(IN) :: omega
        REAL(dp), DIMENSION(SIZE(A, 1)), INTENT(OUT) :: x
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x0(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x(i) = omega*(x(i) / A(i, i) - x0(i)) + x0(i)
        END DO

    END SUBROUTINE JOR

    !> strongly implicit procedure (SIP) method (or stone's method)
    !>
    !> This subroutine implements the SIP method for solving linear systems.
    !> It uses the incomplete LU decomposition of the matrix A.
    SUBROUTINE SIP_ILU(L, U, x0, x, r, omega)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: L, U
        REAL(dp), DIMENSION(:), INTENT(IN) :: r, x0
        REAL(dp), DIMENSION(:), INTENT(OUT) :: x
        REAL(dp), INTENT(IN) :: omega
        REAL(dp), DIMENSION(SIZE(L, 1)) :: y, z

        y = forward(L, r)
        
        z = backward(U, y)
        
        x = x0 + omega * z

    END SUBROUTINE SIP_ILU

    !> strongly implicit procedure (SIP) method (or stone's method)
    !>
    !> This subroutine implements the SIP method for solving linear systems.
    !> It uses the incomplete Cholesky decomposition of the matrix A.
    SUBROUTINE SIP_ICF(L, x0, x, r, omega)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: L
        REAL(dp), DIMENSION(:), INTENT(IN) :: r, x0
        REAL(dp), DIMENSION(:), INTENT(OUT) :: x
        REAL(dp), INTENT(IN) :: omega
        REAL(dp), DIMENSION(SIZE(L, 1)) :: y, z

        y = forward(L, r)
    
        z = backward(TRANSPOSE(L), y)
        
        x = x0 + omega * z

    END SUBROUTINE SIP_ICF

    !> Symmetric successive Over-Relaxation (SSOR) iterative method
    !>
    !> This subroutine implements the SSOR method for solving linear systems.
    SUBROUTINE SSOR(A, b, x0, x, omega)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x0
        REAL(dp), INTENT(IN) :: omega
        REAL(dp), DIMENSION(SIZE(A, 1)), INTENT(OUT) :: x
        REAL(dp), DIMENSION(SIZE(A, 1)) :: x_tmp
        INTEGER :: i, N

        N = SIZE(A, 1)

        ! forward
        DO i = 1, N
            x_tmp(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x_tmp(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x0(i+1:N))
            x_tmp(i) = omega*(x_tmp(i) / A(i, i) - x0(i)) + x0(i)
        END DO

        ! backward
        DO i = N, 1, -1
            x(i) = b(i) - DOT_PRODUCT(A(i, 1:i-1), x_tmp(1:i-1)) - DOT_PRODUCT(A(i, i+1:N), x(i+1:N))
            x(i) = omega*(x(i) / A(i, i) - x_tmp(i)) + x_tmp(i)
        END DO

    END SUBROUTINE SSOR


END MODULE NAFPack_Iterative_methods