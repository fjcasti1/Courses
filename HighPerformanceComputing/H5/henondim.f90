module henondim
!  HENONDIM - Implement Algorithm D to determine the dimension of the basin of
!  infinity of the Henon map.
!
!  SYNOPSIS
!    use henondim
!
!  DESCRIPTION
!    The program testing this algorithm (test_henon) will pass a 1D array of
!    epsilons, the dimension of this array, and the dimension d. The basin is
!    calculated using the henon_map subroutine. Then the code, when running in
!    parallel, will thread assigning a different epsilon to each thread. Next,
!    for each epsilon, the subroutine henon_compare will compute basinP and
!    basinM (the basins obtained of the grid shifted by + or - epsilon) and
!    count the number of epsilon-uncertain grid points. To finish, it will use
!    the subroutine linfit of the module lsq, to obtain the dimension (if the
!    no erro flah is returned), which will be then passed to the testing
!    program test_henon.
!
!  PUBLIC ROUTINES DEFINED
!   - henon_map
!   - basin_alg
!   - basin_compare
!
!  REQUIRED DEPENDENCIES
!    precision - to define KINDs for single and double floating-point.
!          lsq - to do the power fitting needed at the end of th Algorithm D. 
!
!  REVISION HISTORY
!   10/19/18 - Code works using constructed Matlab functions linspace and
!              meshgrid. Also the subroutine reshape since 2D arrays basins
!              are used.
!   10/24/18 - Removed all Matlab functions and reshaping. Adapted code to work
!              with 1D array basins of dimension NGRID^2. Included subroutine
!              basin_compare in order to use OpenMP. The subroutine will make
!              sure each thread has individual copies of the basins to work with
!              them.
!
!  PROGRAMMER
!   Francisco Castillo-Carrasco, fjcasti1@asu.edu
!
  use lsq
  use precision
  implicit none
  real(DP), parameter:: LOCKOUT=100
  real(DP), parameter:: BOXMIN=-3.0d0, BOXMAX=3.0d0  ! box limits
  integer, parameter:: NGRID=4096  ! number of points on a side
  integer, parameter:: K=100  ! maximum number of iterations 
  contains
!-------------------------------------------------------------------------------
  subroutine henon_map(basin,eps)
!  Iterate the Henon map for one or more initial conditions up to
!  MAXITER times or until the orbit of one of the points is further than
!  LOCKOUT units away from the origin.  (It suffices to check whether the
!  |X_n| > LOCKOUT on the nth iteration.)
!
! ARGUMENTS - basin : logical, 1D array of dimension (NGRID*NGRID), True will
!                     mean that the grid point tends to infinity.
!           - eps   : real, eps represents the magnitude of the shifting
!                     of the initial grid.
!
! PROCEDURE : We define the initial grid using array constructors and do loops
!             and store it in the 2D array X0. This array will evolve according
!             to the equations giveni in the Henon map. The grid points will
!             be filtered out using the commands where and all.

    implicit none
    logical , intent(inout) :: basin(NGRID*NGRID)
    real(DP), intent(in)    :: eps

    integer  :: j
    real(DP) :: a = 2.12d0, b = -0.3d0, h
    real(DP), dimension (:)  , allocatable :: x, y
    real(DP), dimension (:)  , allocatable :: aux 
    real(DP), dimension (:,:), allocatable :: X0 

    allocate(x(NGRID))
    allocate(y(NGRID))
    allocate(X0(NGRID*NGRID,2))
    allocate(aux(NGRID*NGRID))

    h = (BOXMAX-BOXMIN)/(NGRID-1)
    x = (/ (h*(j-1)+BOXMIN+eps, j=1,NGRID)  /)
    y = (/ (h*(j-1)+BOXMIN, j=1,NGRID)  /)

    do j=1,NGRID
      X0(1+(j-1)*NGRID:j*NGRID,1)=x(j)
      X0(1+(j-1)*NGRID:j*NGRID,2)=y(:)
    end do

    basin=all(abs(X0)>LOCKOUT,2)

    do j=1,K
      where(basin.eq. .False.)
        aux(:) = X0(:,1)
        X0(:,1) = a-X0(:,1)**2d0+b*X0(:,2)
        X0(:,2) = aux(:)
        basin=all(abs(X0)>LOCKOUT,2)
      end where
    end do

  end subroutine henon_map
!-------------------------------------------------------------------------------
  subroutine basin_compare(N,eps,basin)
! BASIN_COMPARE - counts the number of epsilon uncertain grid points.
!
! ARGUMENTS - N     : real, N represents the number of epsilon-uncertain points
!                     for the given value of epsilon, eps.
!           - eps   : real, eps represents the magnitude of the shifting of the
!                     initial grid.
!           - basin : logical, array of dimension NGRID^2, it is the original
!                     unshifted basin, used to compare with the shifted ones and
!                     obtain the number of epsilon-uncertain grid points.
!
! PROCEDURE : We simply use the henon_map subroutine to calculate the shifted
!             basins and the function count to obtain the total number of grid
!             points whose basins don't match. This is, the total number of
!             epsilon-uncertain grid points.

    implicit none
    real(DP), intent(in)    :: eps
    real(DP), intent(inout) :: N
    logical, intent(in)     :: basin(NGRID*NGRID)
    logical, dimension(:), allocatable  :: basinP, basinM
    integer :: j

    allocate(basinP(NGRID*NGRID))
    allocate(basinM(NGRID*NGRID))

    call henon_map(basinP,eps)  
    call henon_map(basinM,-eps)  
    N = count(basin.ne.basinP.or.basin.ne.basinM)
!    N=0
!    do j=1,NGRID*NGRID
!      if (basin(j).ne.basinP(j).or.basin(j).ne.basinM(j)) then
!        N=N+1
!      end if
!    end do
  end subroutine basin_compare
!-------------------------------------------------------------------------------
  subroutine basin_alg(Neps,eps,d)
     
! BASIN_ALG - implements Algorithm D.
!
! ARGUMENTS - Neps : integer, dimension of 1D arrays N and eps.
!           - eps  : real, 1D array eps which represents the magnitude of the
!                    shifting of the initial grid.
!           - d    : real, d represents the dimension of the basin of attraction. 
!
! PROCEDURE : We obtain the basin of attraction of the Henon map for the
!             unshifted grid, basin. Then, for multiple epsilons, we compare
!             the basins obtained out of shifted grids by a magnitued epsilon.
!             This is done in the subroutine basin_compare in such a way that
!             can be parallelizable using OpenMP. Such subroutine returns then
!             two 1D arrays, the number of epsilon uncertain grid points N,
!             and the different epsilon used. Those two arrays are passed to
!             the linfit subroutine to obtain the dimension of the basin of
!             attraction.

    integer,  intent(in)  :: Neps
    real(DP), intent(inout) :: eps(Neps), d

    integer  :: j, ierr
    real(DP) :: param(2), N(Neps)
    logical, dimension(:), allocatable  :: basin
    
    allocate(basin(NGRID*NGRID))

    call henon_map(basin,0.d0)  
    !$omp parallel do 
    do j=1,Neps
      call basin_compare(N(j),eps(j),basin)
    enddo
    !$omp end parallel do

    call linfit('Power',Neps,eps,N,param,ierr)    
    if(ierr==0) d=1+param(2)
    return
  end subroutine basin_alg
!-------------------------------------------------------------------------------
end module henondim
