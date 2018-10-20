module henondim
!  HENONDIM - Implement Algorithm D to determine the dimension of the basin of
!  infinity of the Henon map.
!
!  SYNOPSIS
!    use henondim
!
!  DESCRIPTION
!    Write a paragraph that summarizes your computational approach to implement
!    Algorithm D.
!
!  PUBLIC ROUTINES DEFINED
!    Brief list of procedures that can be called from outside this module.
!
!  REQUIRED DEPENDENCIES
!    precision - to define KINDs for single and double floating-point.
!    And any other modules that you define.
!
!  REVISION HISTORY
!   10/15/18 - First implementation.
!
!  PROGRAMMER
!   Francisco Castillo, fjcasti1@asu.edu
!
  use matlabfns
  use myStructModule
  use precision
  implicit none
  real(DP), parameter:: LOCKOUT=100
  real(DP), parameter:: BOXMIN=-3.0, BOXXMAX=3.0  ! box limits
  integer, parameter:: NGRID=4096  ! number of points on a side
  contains
!-------------------------------------------------------------------------------
  subroutine henon_map(N,xMin,xMax,yMin,yMax,params)
!  Iterate the Henon map for one or more initial conditions up to
!  MAXITER times or until the orbit of one of the points is further than
!  LOCKOUT units away from the origin.  (It suffices to check whether the
!  |X_n| > LOCKOUT on the nth iteration.)
!
!  Argument declarations go here.  Explicitly declare all arguments and give
!  each of them an INTENT.
!
!
!  Declare local variables here.  You may declare the map parameters A and B
!  as Fortran PARAMETERS or simple local variables at your option.
!
    implicit none
    integer  :: N, j
    real(DP) :: xMin, xMax, yMin, yMax
    real(DP) :: x(N), y(N)
    real(DP), dimension (:,:), allocatable :: xGrid, yGrid
    real(DP), dimension (:,:), allocatable :: xGrid2, yGrid2
    real(DP), dimension (:,:), allocatable :: X0, aux 
    integer, dimension (:,:), allocatable  :: basin
    logical, dimension (:)  , allocatable  :: vbasin
    type(paramStruct) :: params
    logical :: ix(N*N)

    call linspace(x,xMin,xMax,N)
    call linspace(y,yMin,yMax,N)
!    print*, x
  !!!  print*, y
    call meshgrid(xGrid,yGrid,x,y)
!!!    print*, " "
!!!    do j=1,N
!!!    print*, xGrid(j,:)
!!!    end do
!!!    print*, " "
!!!    do j=1,N
!!!    print*, yGrid(j,:)
!!!    end do
    xGrid2 = reshape(xGrid,(/ N*N,1 /) )
    yGrid2 = reshape(yGrid,(/ N*N,1 /) )
  !!!  print*, " "
  !!!  print*, xGrid2
  !!!  print*, " "
  !!!  print*, yGrid2
    allocate(X0(N*N,2))
    allocate(aux(N*N,1))
    allocate(basin(N,1))
    basin=1
    X0(:,1)=xGrid2(:,1)
    X0(:,2)=yGrid2(:,1)
    print *, " "
    print *, "INITIAL:"
    do j=1,N*N
      print*, X0(j,:) 
    end do
    ix=all(abs(X0)<params%L,2)
    print*, all(abs(X0)<params%L,2)
!!!    do j=1,N*N
!!!      print*, X0(j,:)
!!!    end do
    print *, "K = ",params%K
    print *, "L = ",params%L
    do j=1,params%K
      where(ix.eq. .True.)
        aux(:,1) = X0(:,1)
        X0(:,1) = params%a-X0(:,1)**2d0+params%b*X0(:,2)
        X0(:,2) = aux(:,1)
!!!        aux(:,1) = X0(:,1)
!!!        X0(:,1) = X0(:,2)**2d0
!!!        X0(:,2) = aux(:,1)
        ix=all(abs(X0)<params%L,2)
      end where
      print *, ix
    end do
    do j=1,N*N
      print*, X0(j,:)
    end do
    print *, " "
    print *, params%L
    print *, "Final: ", ix

    basin = abs(reshape(ix,(/ N,N /) ))

    print *, " "
    do j=1,N
      print*, basin(j,:) 
    end do
  !  return
  end subroutine henon_map
!-------------------------------------------------------------------------------
!!!  subroutine basin_alg(n,xMin,xMax,yMin,yMax,params)
  
! BASIN_ALG - implements Algorithm D.
!
!  Argument declarations go here.  Explicitly declare all arguments and give
!  each of them an INTENT.
!
!
!  Local variables go here.
!

!  The rest of your code goes here.
!
!!! return
!!! end subroutine basin_alg
!-------------------------------------------------------------------------------
end module henondim
