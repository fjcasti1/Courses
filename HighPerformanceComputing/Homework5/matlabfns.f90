module matlabfns
implicit none
contains
    subroutine meshgrid(xgv, ygv, zgv, X, Y, Z)
    ! Found from :
    ! https://stackoverflow.com/questions/37090616/fortran-equivalent-of-the-matlab-meshgrid-function
        implicit none
        real,intent(in)   :: xgv(:), ygv(:), zgv(:)
        real,intent(out)  :: X(:,:,:), Y(:,:,:), Z(:,:,:)
        integer           :: sX, sY, sZ, i

        sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)

        do i=1,sZ
        X(:,:,i) = spread( xgv, 1, sY )
        Y(:,:,i) = spread( ygv, 2, sX )
        enddo ! i
        do i=1,sX
        Z(i,:,:) = spread( zgv, 1, sY)
        enddo ! i
    end subroutine

    SUBROUTINE linspace(z, l, k, n)
    ! Found from : 
    ! http://folk.uio.no/gunnarw/CSE-FL/Fortran/PART5/cse-fl.pdf
        IMPLICIT NONE
        !// Argument declarations
        COMPLEX, DIMENSION(n)   :: z
        INTEGER                 :: l
        INTEGER                 :: k
        INTEGER                 :: n
        !// Local variables
        INTEGER                 :: i
        REAL                    :: d, x, y

        x = FLOAT(k)
        y = FLOAT(l)
        d = (x - y) / n
        PRINT *, d, k, l, n
        z(1) = FLOAT(l)
        DO i = 2, n-1
        z(i) = z(i-1) + d
        END DO
        z(1) = y
        z(n) = x
        RETURN
    END SUBROUTINE linspace
end module matlabfns
