module matlabfns
  use precision
  implicit none
  contains
    SUBROUTINE linspace(x, xMin, xMax, N)
    ! Found from : 
    ! http://folk.uio.no/gunnarw/CSE-FL/Fortran/PART5/cse-fl.pdf
      IMPLICIT NONE
      !// Argument declarations
      integer  :: N, j
      real(DP) :: xMin, xMax, h
      real(DP) :: x(1:N)

      !// Local variables
      h = (xMax - xMin) / (N-1)
      DO j = 1, N
        x(j) = h*(j-1)+xMin
      END DO
    END SUBROUTINE linspace

    subroutine meshgrid(x,y,xgv,ygv, ierr)
        !..............................................................................
        !meshgrid generate mesh grid over a rectangular domain of [xmin xmax, ymin, ymax]
        ! Inputs:
        !     xgv, ygv are grid vectors in form of full grid data
        ! Outputs:
        !     X and Y are matrix each of size [ny by nx] contains the grid data.
        !     The coordinates of point (i,j) is [X(i,j), Y(i,j)]
        !     ierr: The error flag
        !     """
        !     # Example
        !     # call meshgrid(X, Y, [0.,1.,2.,3.],[5.,6.,7.,8.])
        !     # X
        !     # [0.0, 1.0, 2.0, 3.0,
        !     #  0.0, 1.0, 2.0, 3.0,
        !     #  0.0, 1.0, 2.0, 3.0,
        !     #  0.0, 1.0, 2.0, 3.0]
        !     #
        !     #Y
        !     #[ 5.0, 5.0, 5.0, 5.0,
        !     #  6.0, 6.0, 6.0, 6.0,
        !     #  7.0, 7.0, 7.0, 7.0,
        !     #  8.0, 8.0, 8.0, 8.0]
        !..............................................................................
        ! Rev 0.2, Feb 2018
        ! New feature added: xgv and ygv as full grid vector are accepted now

        ! Arguments
        real(DP), intent(out), allocatable  :: x(:,:)
        real(DP), intent(out), allocatable  :: y(:,:)
        real(DP), intent(in)                :: xgv(:) ! x grid vector [start, stop, step] or [start, stop]
        real(DP), intent(in),  optional     :: ygv(:) ! y grid vector [start, stop, step] or [start, stop]
        integer,  intent(out), optional     :: ierr   ! the error value

        ! Local variables
        integer:: sv
        integer:: nx
        integer:: ny
        logical:: only_xgv_available

        ! Initial setting
        only_xgv_available  = .false.
        sv=0 !Assume no error

        nx=size(xgv, dim=1)

        if (present(ygv)) then
            ny = size(ygv, dim=1)
        else
            only_xgv_available=.true.
            ny=nx
        end if

        allocate(x(ny,nx),y(ny,nx),stat=sv)
        if (sv /=0) then
            print*, "allocataion erro in meshgrid"
            stop
        end if

        x(1,:)    = xgv
        x(2:ny,:) = spread(xgv, dim=1, ncopies=ny-1)

        if (only_xgv_available) then
            y=transpose(x)
        else
            y(:,1)    = ygv
            y(:,2:nx) = spread(ygv,dim=2,ncopies=nx-1)
        end if

        if (present(ierr)) then
            ierr=sv
        end if

    end subroutine meshgrid
end module matlabfns
