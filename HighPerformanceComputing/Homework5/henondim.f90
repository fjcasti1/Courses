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
  integer, parameter:: NGRID=4096  ! number of points on a side
!  integer, parameter:: NGRID=1000  ! number of points on a side
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
    real(DP), dimension (:)  , allocatable :: x, y
    real(DP), dimension (:,:), allocatable :: xGrid, yGrid
    real(DP), dimension (:,:), allocatable :: X0, aux 
    real(DP) :: a = 2.12, b = -0.3, eps
    logical, intent(inout) :: basin(NGRID,NGRID)
    logical, dimension(:), allocatable :: ix

    allocate(x(NGRID))
    allocate(y(NGRID))
    allocate(xGrid(NGRID,NGRID))
    allocate(yGrid(NGRID,NGRID))
    allocate(X0(NGRID*NGRID,2))
    allocate(aux(NGRID*NGRID,1))
    allocate(ix(NGRID*NGRID))

    call linspace(x,BOXMIN+eps,BOXMAX+eps,NGRID)
    call linspace(y,BOXMIN,BOXMAX,NGRID)
    call meshgrid(xGrid,yGrid,x,y)
    X0(:,1) = reshape(xGrid,(/ NGRID*NGRID /) )
    X0(:,2) = reshape(yGrid,(/ NGRID*NGRID /) )
    ix=all(abs(X0)<LOCKOUT,2)
    do j=1,K
      where(ix.eq. .True.)
        aux(:,1) = X0(:,1)
        X0(:,1) = a-X0(:,1)**2d0+b*X0(:,2)
        X0(:,2) = aux(:,1)
        ix=all(abs(X0)<LOCKOUT,2)
      end where
    end do
    ix = .not. ix
    basin = reshape(ix,(/ NGRID,NGRID /) )
  end subroutine henon_map
!-------------------------------------------------------------------------------
  subroutine henon_compare(N,eps,basin)
    real(DP), intent(in)    :: eps
    real(DP), intent(inout) :: N
    logical, intent(in)  :: basin(NGRID,NGRID)
    logical, dimension(:,:), allocatable  :: basinP, basinM
    
    allocate(basinP(NGRID,NGRID))
    allocate(basinM(NGRID,NGRID))

    call henon_map(basinP,eps)  
    call henon_map(basinM,-eps)  
    N = count(basin.ne.basinP.or.basin.ne.basinM)
  end subroutine henon_compare

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
    logical, dimension(:,:), allocatable  :: basin
    
    allocate(basin(NGRID,NGRID))

    call henon_map(basin,0d0)  
    !$omp parallel do 
    do j=1,Neps
      call henon_compare(N(j),eps(j),basin)
    enddo
    !$omp end parallel do

   ! print *, " "
   ! print *, "N = "
   ! do j=1,Neps
   !   print*, N(j) 
   ! end do
   ! print *, " "

    call linfit('Power',Neps,eps,N,param,ierr)    
    if(ierr==0) print *, "d = ",1+param(2)
    return
  end subroutine basin_alg
!-------------------------------------------------------------------------------
end module henondim
