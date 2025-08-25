module NAFPack_matrix_market

    use NAFPack_constant, only: dp

    implicit none(type, external)

    private
    public :: readMatrixMarket

contains

    subroutine readMatrixMarket(filename, A)
        character(LEN=*), intent(in) :: filename
        real(dp), dimension(:, :), allocatable, intent(out) :: A

        integer :: i, j, nrows, ncols, nnz, ios, k, row_idx, col_idx
        character(LEN=256) :: line, header
        character(LEN=50) :: object, format_type, field, symmetry
        integer :: unit
        real(dp) :: val

        unit = 10

        ! Open the file
        open (NEWUNIT=unit, FILE=filename, STATUS='old', ACTION='read', IOSTAT=ios)
        if (ios /= 0) then
            print*,'Error: opening file: ', filename
            stop
        end if

        ! read the header
        read (unit, '(A)', IOSTAT=ios) header
        if (ios /= 0) then
            print*,'Error: reading header Matrix Market'
            close (unit)
            stop
        end if

        ! Check Matrix Market format
        if (index(header, '%%MatrixMarket') == 0) then
            print*,'Error: file is not in Matrix Market format'
            close (unit)
            return
        end if

        ! Parse the header to extract information
        read (header, *, iostat=ios) object, object, format_type, field, symmetry
        if (ios /= 0) then
            print*,'Error reading Matrix Market header'
            close (unit)
            stop
        end if

        ! Skip comment lines and find the size line
        do
            read (unit, '(A)', IOSTAT=ios) line
            if (ios /= 0) then
                print*,'Error: malformed file'
                close (unit)
                stop
            end if
            if (line(1:1) /= '%') then
                exit
            end if
        end do

        ! Read matrix dimensions
        read (line, *, iostat=ios) nrows, ncols, nnz
        if (ios /= 0) then
            print*,'Error: unable to read dimensions'
            close (unit)
            return
        end if

        allocate (A(nrows, ncols))
        A = 0.0_dp

        ! Read matrix entries
        if (trim(format_type) == 'coordinate') then
            ! Coordinate format (sparse)
            do k = 1, nnz
                read (unit, *, iostat=ios) row_idx, col_idx, val
                if (ios /= 0) then
                    print*,'Error: reading entry', k
                    exit
                end if

                A(row_idx, col_idx) = val

                ! If the matrix is symmetric, also fill A(j,i)
                if (trim(symmetry) == 'symmetric' .and. row_idx /= col_idx) then
                    A(col_idx, row_idx) = val
                end if
            end do

        else if (trim(format_type) == 'array') then
            ! Dense format (array)
            do j = 1, ncols
                do i = 1, nrows
                    read (unit, *, iostat=ios) val
                    if (ios /= 0) then
                        print*,'Error: reading element (', i, ',', j, ')'
                        close (unit)
                        return
                    end if
                    A(i, j) = val
                end do
            end do

        else
            print*,'Error: unsupported format: ', trim(format_type)
            close (unit)
            return
        end if

        close (unit)
    end subroutine readMatrixMarket

end module NAFPack_matrix_market
