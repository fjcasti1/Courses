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
  use lsq
  use precision
  implicit none
  real(DP), parameter:: LOCKOUT=100
  real(DP), parameter:: BOXMIN=-3.0, BOXMAX=3.0  ! box limits
!  integer, parameter:: NGRID=4096  ! number of points on a side
  integer, parameter:: NGRID=100  ! number of points on a side
  integer, parameter:: K=100  ! maximum number of iterations 
  contains
!-------------------------------------------------------------------------------
  subroutine henon_map(basin,eps)
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
    integer  :: j
    real(DP) :: x(NGRID), y(NGRID)
    real(DP), dimension (:,:), allocatable :: xGrid, yGrid
    real(DP), dimension (:,:), allocatable :: X0, aux 
    real(DP) :: a = 2.12, b = -0.3, eps
    logical, intent(inout) :: basin(NGRID,NGRID)
    logical :: ix(NGRID*NGRID)

    call linspace(x,BOXMIN+eps,BOXMAX+eps,NGRID)
    call linspace(y,BOXMIN,BOXMAX,NGRID)
!    print*, x
  !!!  print*, y
    call meshgrid(xGrid,yGrid,x,y)
!!!    print*, " "
!!!    do j=1,NGRID
!!!    print*, xGrid(j,:)
!!!    end do
!!!    print*, " "
!!!    do j=1,NGRID
!!!    print*, yGrid(j,:)
!!!    end do
!!rm    xGrid2 = reshape(xGrid,(/ NGRID*NGRID,1 /) )
!!rm    yGrid2 = reshape(yGrid,(/ NGRID*NGRID,1 /) )
  !!!  print*, " "
  !!!  print*, xGrid2
  !!!  print*, " "
  !!!  print*, yGrid2
    allocate(X0(NGRID*NGRID,2))
    allocate(aux(NGRID*NGRID,1))
    X0(:,1) = reshape(xGrid,(/ NGRID*NGRID /) )
    X0(:,2) = reshape(yGrid,(/ NGRID*NGRID /) )
!!rm    X0(:,1)=xGrid2(:,1)
!!rm    X0(:,2)=yGrid2(:,1)
!!    print *, " "
!!    print *, "INITIAL:"
!    do j=1,NGRID*NGRID
!      print*, X0(j,:) 
!    end do
    ix=all(abs(X0)<LOCKOUT,2)
!    print*, all(abs(X0)<LOCKOUT,2)
!!!    do j=1,NGRID*NGRID
!!!      print*, X0(j,:)
!!!    end do
!    print *, "K = ",K
!    print *, "L = ",LOCKOUT
    do j=1,K
      where(ix.eq. .True.)
        aux(:,1) = X0(:,1)
        X0(:,1) = a-X0(:,1)**2d0+b*X0(:,2)
        X0(:,2) = aux(:,1)
!!!        aux(:,1) = X0(:,1)
!!!        X0(:,1) = X0(:,2)**2d0
!!!        X0(:,2) = aux(:,1)
        ix=all(abs(X0)<LOCKOUT,2)
      end where
!      print *, ix
    end do
!    do j=1,NGRID*NGRID
!      print*, X0(j,:)
!    end do
!    print *, " "
!    print *, LOCKOUT
!    print *, "Final: ", ix
    ix = .not. ix
    basin = reshape(ix,(/ NGRID,NGRID /) )

!    print *, " "
!    do j=1,NGRID
!      print*, basin(j,:) 
!    end do
  !  return
  end subroutine henon_map
!-------------------------------------------------------------------------------
  subroutine basin_alg(Neps,N,eps)
     
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
    integer  :: j, Neps, ierr
    real(DP) :: eps(Neps), param(2), N(Neps)
    logical  :: basin(NGRID,NGRID)
    logical  :: basinP(NGRID,NGRID), basinM(NGRID,NGRID)

    do j=1,Neps
      call henon_map(basin,0d0)  
      call henon_map(basinP,eps(j))  
      call henon_map(basinM,-eps(j))  
      N(j) = count(basin.ne.basinP.or.basin.ne.basinM)
    end do

!!!    print *, " "
!!!    do j=1,NGRID
!!!      print*, basin(j,:) 
!!!    end do
!!!    print *, " "
!!!    do j=1,NGRID
!!!      print*, basinP(j,:) 
!!!    end do
!!!    print *, " "
!!!    do j=1,NGRID
!!!      print*, basinM(j,:) 
!!!    end do
    print *, " "
    print *, "N = "
    do j=1,Neps
      print*, N(j) 
    end do
    print *, " "
!    print *, count(basin.ne.basinP)
!    print *, count(basin.ne.basinM)
!    print *, count(basin.ne.basinP.or.basin.ne.basinM)
    call linfit('Power',Neps,eps,N,param,ierr)    
    if(ierr==0) print *, param
    return
  end subroutine basin_alg
!-------------------------------------------------------------------------------
end module henondim
