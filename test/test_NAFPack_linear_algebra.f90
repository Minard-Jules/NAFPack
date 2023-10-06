MODULE test_NAFPack_linear_algebra

    USE NAFPack_linear_algebra
    USE NAFPack_constant

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: test_linear_algebra
    PUBLIC :: test_linear_system

    CONTAINS

!================ Linear System =========================================================== 

    SUBROUTINE test_methode_direct(A, b, x, method, stat)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x
        CHARACTER(LEN = *), INTENT(IN) :: method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(SIZE(x)) :: diff_x, x_tmp
    
        x_tmp = direct_methode(A, b, method = method)
    
        diff_x = x - x_tmp
        IF (MAXVAL(ABS(diff_x)) < epsi_test) THEN
            WRITE(*,'(A,T40,A,A)') green_color//method, " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//method, " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF
    
    END SUBROUTINE test_methode_direct

    SUBROUTINE test_iterative_methods(A, b, x, method, stat)
        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        REAL(dp), DIMENSION(:), INTENT(IN) :: b, x
        CHARACTER(LEN = *), INTENT(IN) :: method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(SIZE(x)) :: diff_x, x_tmp
    
        x_tmp = Iterative_methods(A, b, method = method)
        
        diff_x = x - x_tmp
        IF (MAXVAL(ABS(diff_x)) < epsi_test) THEN
            WRITE(*,'(A,T40,A,A)') green_color//method, " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//method, " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF
    
    END SUBROUTINE test_iterative_methods
    

    SUBROUTINE test_linear_system(stat)

        LOGICAL, INTENT(INOUT) :: stat
        
        INTEGER, PARAMETER :: N = 3

        REAL(dp), DIMENSION(N, N) :: A
        REAL(dp), DIMENSION(N) :: b
        REAL(dp), DIMENSION(N) :: x

        A = RESHAPE([1, -1, 2, -1, 5, -4, 2, -4, 6], [N, N])
        b = [1, 1, 2]
        x = [0, 1, 1]

        WRITE(*,'(A, A)') cyan_color//"Direct methode"// reset_color
        PRINT*, " "

        !==================================================================
        !Gauss method
        CALL test_methode_direct(A, b, x, "Gauss", stat)

        !==================================================================
        !Gauss pivot method
        CALL test_methode_direct(A, b, x, "Gauss_pivot", stat)

        !==================================================================
        !LU decomposition method
        CALL test_methode_direct(A, b, x, "A_LU", stat)

        !==================================================================
        !LDU decomposition method
        CALL test_methode_direct(A, b, x, "A_LDU", stat)

        !==================================================================
        !Cholesky decomposition method
        CALL test_methode_direct(A, b, x, "Cholesky", stat)

        !==================================================================
        !QR_Householder decomposition method
        CALL test_methode_direct(A, b, x, "QR_Householder", stat)

        !==================================================================
        !QR_Givens decomposition method
        CALL test_methode_direct(A, b, x, "QR_Givens", stat)

        !==================================================================
        !QR_Gram_Schmidt_Classical decomposition method
        CALL test_methode_direct(A, b, x, "QR_Gram_Schmidt_Classical", stat)

        !==================================================================
        !QR_Gram_Schmidt_Modified decomposition method
        CALL test_methode_direct(A, b, x, "QR_Gram_Schmidt_Modified", stat)


        !================Iterative_methods=============================
        
        PRINT*, " "
        PRINT*, " "
        WRITE(*,'(A, A)') cyan_color//"Iterative_methods"// reset_color
        PRINT*, " "

        A = RESHAPE([5, 3, 1, -1, 8, 1, 2, -2, 4], [N, N])
        b = [12, -25, 6]
        x = [1, -3, 2]

        !==================================================================
        !Jacobi method
        CALL test_iterative_methods(A, b, x, "Jacobi", stat)

        !==================================================================
        !Gauss_Seidel method
        CALL test_iterative_methods(A, b, x, "Gauss_Seidel", stat)

        !==================================================================
        !SOR method
        CALL test_iterative_methods(A, b, x, "SOR", stat)

    END SUBROUTINE test_linear_system

!================ Eigen ===================================================================    

    SUBROUTINE test_Eigen_vecteur_propre(A, method, stat)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), INTENT(IN) :: method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION (SIZE(A, 1), SIZE(A, 2)) :: vp
        REAL(dp), DIMENSION (SIZE(A, 1)) :: lambda, diff
        INTEGER :: i
        LOGICAL :: verif_Eigen = .TRUE.

        WRITE(*,'(A, A)') cyan_color//"Eigen value"// reset_color
        PRINT*, " "

        CALL Eigen(A, lambda, vp = vp, method = method)

        DO i = 1, SIZE(A, 1)
            diff = MATMUL(A, vp(i, :)) - lambda(i) * vp(i, :)
            IF(MAXVAL(ABS(diff)) > epsi_test) verif_Eigen = .FALSE.
        END DO

        IF(verif_Eigen)THEN
            WRITE(*,'(A,T40,A,A)') green_color//method, " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//method, " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF

    END SUBROUTINE test_Eigen_vecteur_propre

    SUBROUTINE test_Eigen_without_vecteur_propre(A, method, stat)

        REAL(dp), DIMENSION(:, :), INTENT(IN) :: A
        CHARACTER(LEN = *), INTENT(IN) :: method
        LOGICAL, INTENT(INOUT) :: stat
        REAL(dp), DIMENSION(SIZE(A, 1)) :: lambda_test, lambda, diff
        CHARACTER(LEN = 25) :: method_test
        LOGICAL :: verif_Eigen = .TRUE.

        method_test = "Power_iteration"


        IF(stat)THEN
            verif_Eigen = .FALSE.
        ELSE
            lambda_test = [-4, 2, 1]
            CALL Eigen(A, lambda, method = method)

            diff = lambda_test - lambda
            IF(MAXVAL(ABS(diff)) > epsi_test) verif_Eigen = .FALSE.
        END IF

        IF(verif_Eigen)THEN
            WRITE(*,'(A,T40,A,A)') green_color//method, " :: OK"// reset_color
        ELSE
            WRITE(*,'(A,T40,A)') red_color//method, " :: ECHEC"// reset_color
            stat = .TRUE.
        END IF


    END SUBROUTINE test_Eigen_without_vecteur_propre

    SUBROUTINE test_linear_algebra(stat)

        LOGICAL, INTENT(INOUT) :: stat
        LOGICAL :: stat_tmp = .FALSE.

        INTEGER, PARAMETER :: N = 3
        REAL(dp), DIMENSION (N,N) :: A

        A=RESHAPE([0, 3, -2, 2, -2, 2, -1, 0, 1], [N, N])
        CALL test_Eigen_vecteur_propre(A, "Power_iteration", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Householder", stat_tmp)
        CALL test_Eigen_without_vecteur_propre(A, "QR_Givens", stat_tmp)
        !CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Classical", stat_tmp)
        !CALL test_Eigen_without_vecteur_propre(A, "QR_Gram_Schmidt_Modified", stat_tmp)
        IF(stat_tmp) stat = stat_tmp

    END SUBROUTINE test_linear_algebra

END MODULE test_NAFPack_linear_algebra