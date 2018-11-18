module henondim
!  HENONDIM - Implement Algorithm D to determine the dimension of the basin of
!  infinity of the Henon map. Uses MPI.
!
!  SYNOPSIS
!    use henondim
!
!  DESCRIPTION
!    This module contains the necessary subroutines to calculate the number of
!    epsilon-uncertain grid points of the basin calculated for the Henon map.
!    It is parallelized using MPI by assigning columns to each processing element.
!    The code calling this module must get the results and add them together
!    using mpi_reduce.
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
!   11/13/18 - Adapted code to use MPI parallelism. However, each processing
!              element creates the whole domain and selects only a part to work
!              with. It would be better to only create the corresponding section
!              of the domain.
!   11/17/18 - Only the corresponding section of the domain is created and
!              worked with by each processing element, depending on its rank.
!
!  PROGRAMMER
!   Francisco Castillo-Carrasco, fjcasti1@asu.edu
!
  use precision
  implicit none
  contains
!-------------------------------------------------------------------------------
  subroutine henon_map(basin,basindim,xextent,me,npes,params,K,LOCKOUT,eps)
!  Iterate the Henon map for one or more initial conditions up to
!  MAXITER times or until the orbit of one of the points is further than
!  LOCKOUT units away from the origin.  (It suffices to check whether the
!  |X_n| > LOCKOUT on the nth iteration.)
!
! ARGUMENTS - basin   : logical array which dimension depends of the rank of the
!                       processor of dimension (the number of colums times NGRIDy), 
!                       it is the original unshifted basin, used to compare with
!                       the shifted ones and obtain the number of epsilon-uncertain
!                       grid points. True will mean that the grid point tends
!                       to infinity.
!           - basindim: integer, dimension of the basin being worked with by the
!                       processing unit
!           - xextent : real, array with the following parameters:
!                                         Xmin, Xmax, Ymin, Ymax, NGRIDx, NGRIDy
!           - me      : integer, the MPI rank
!           - npes    : integer, the Number of Processing Elements (> 0)
!           - params  : real, parameters for the henon map. 
!           - K       : integer, maxmimum number of iterations for the henon map 
!           - LOCKOUT : real, limit where we separate infinity points from the rest
!           - eps     : real, the magnitude of the shifting of the initial grid.
!
! PROCEDURE : We define the initial grid using array constructors and do loops
!             and store it in the 2D array X0. This array will evolve according
!             to the equations giveni in the Henon map. The grid points will
!             be filtered out using the commands where and all.

    use, intrinsic::iso_fortran_env
    implicit none
    integer,  intent(in)  :: me, npes, K, basindim
    real(DP), intent(in)  :: xextent(6), params(2), eps, LOCKOUT
    logical,  intent(out) :: basin(basindim)

    integer  :: j, N, M, NGRIDx, NGRIDy
    real(DP) :: xMin, xMax, yMin, yMax, hx, hy
    real(DP), dimension(:)  , allocatable :: x, y, auxme
    real(DP), dimension(:,:), allocatable :: X0me

    xMin   = xextent(1)
    xMax   = xextent(2)
    yMin   = xextent(3)
    yMax   = xextent(4)
    NGRIDx = xextent(5) 
    NGRIDy = xextent(6) 

    hx = (xMax-xMin)/(NGRIDx-1)
    hy = (yMax-yMin)/(NGRIDy-1)
    N  = NGRIDx/npes
    M  = MOD(int(NGRIDx),npes)

    allocate(y(NGRIDy))
    y  = [(hy*(j-1)+yMin, j=1,NGRIDy)]

    if (me.lt.M) then
      allocate(x(N+1))
      allocate(X0me((N+1)*NGRIDy,2))
      allocate(auxme((N+1)*NGRIDy))
      x  = [(hx*(j-1)+xMin+eps, j=1+(N+1)*me,(N+1)*(me+1))]
      do j=1,N+1
        X0me(1+(j-1)*NGRIDy:j*NGRIDy,1) = x(j)
        X0me(1+(j-1)*NGRIDy:j*NGRIDy,2) = y(:)
      end do
    else
      allocate(x(N))
      allocate(X0me(N*NGRIDy,2))
      allocate(auxme(N*NGRIDy))
      x  = [(hx*(j-1)+xMin+eps, j=1+N*me+M,N*(me+1)+M)]
      do j=1,N
        X0me(1+(j-1)*NGRIDy:j*NGRIDy,1) = x(j)
        X0me(1+(j-1)*NGRIDy:j*NGRIDy,2) = y(:)
      end do
    endif

    basin=all(abs(X0me)>LOCKOUT,2)

    do j=1,K
      where(basin.eq..False.)
        auxme     = X0me(:,1)
        X0me(:,1) = params(1)-X0me(:,1)**2d0+params(2)*X0me(:,2)
        X0me(:,2) = auxme(:) 
        basin = all(abs(X0me)>LOCKOUT,2)
      end where
    end do

  end subroutine henon_map
!-------------------------------------------------------------------------------
  subroutine basin_compare(basin0,basindim,xextent,me,npes,params,K,LOCKOUT,eps,Nuncert)
! BASIN_COMPARE - counts the number of epsilon uncertain grid points.
!
! ARGUMENTS - basin0  : logical array which dimension depends of the rank of the
!                       processor of dimension (the number of colums times NGRIDy), 
!                       it is the original unshifted basin, used to compare with
!                       the shifted ones and obtain the number of epsilon-uncertain
!                       grid points.
!           - basindim: integer, dimension of the basin being worked with by the
!                       processing unit
!           - xextent : real, array with the following parameters:
!                                         Xmin, Xmax, Ymin, Ymax, NGRIDx, NGRIDy
!           - me      : integer, the MPI rank
!           - npes    : integer, the Number of Processing Elements (> 0)
!           - params  : real, parameters for the henon map. 
!           - K       : integer, maxmimum number of iterations for the henon map 
!           - LOCKOUT : real, limit where we separate infinity points from the rest
!           - eps     : real, the magnitude of the shifting of the initial grid.
!           - Nuncert : integer, the number of uncertain points
! ARGUMENTS - N     : real, N represents the number of epsilon-uncertain points
!                     for the given value of epsilon, eps.
!           - eps   : real, eps represents the magnitude of the shifting of the
!                     initial grid.
!
! PROCEDURE : We simply use the henon_map subroutine to calculate the shifted
!             basins and the function count to obtain the total number of grid
!             points whose basins don't match. This is, the number of
!             epsilon-uncertain grid points for the section that the given
!             processor is working with.
    implicit none
    integer,  intent(in)  :: me, npes, K, basindim
    integer,  intent(out) :: Nuncert
    real(DP), intent(in)  :: xextent(6), params(2), eps, LOCKOUT
    logical,  intent(in)  :: basin0(basindim)

    logical, dimension(:), allocatable :: basinP, basinM

    allocate(basinP(basindim))
    allocate(basinM(basindim))
! Obtain the basin for the shifted grid an amount eps to the right
    call henon_map(basinP,basindim,xextent,me,npes,params,K,LOCKOUT,eps)
! Obtain the basin for the shifted grid an amount eps to the left
    call henon_map(basinM,basindim,xextent,me,npes,params,K,LOCKOUT,-eps)
! Compare the basins and obtain the number of uncertain grid points.
    Nuncert = count(basin0.ne.basinP.or.basin0.ne.basinM)
  end subroutine basin_compare
!-------------------------------------------------------------------------------
  subroutine basin_alg(xextent,me,npes,params,K,LOCKOUT,eps,epsdim,uncert)
     
! BASIN_ALG - implements Algorithm D, using MPI.
!
! ARGUMENTS - xextent : real, array with the following parameters:
!                                         Xmin, Xmax, Ymin, Ymax, NGRIDx, NGRIDy
!           - me      : integer, the MPI rank
!           - npes    : integer, the Number of Processing Elements (> 0)
!           - params  : real, parameters for the henon map. 
!           - K       : integer, maxmimum number of iterations for the henon map 
!           - LOCKOUT : real, limit where we separate infinity points from the rest
!           - eps     : real, 1D array eps which represents the magnitude of the
!                        shifting of the initial grid.
!           - epsdim  : integer, dimension of the array of epsilons. 
!           - uncert  : integer, 1D array which represents the number of uncertain points

! PROCEDURE : We obtain the basin of attraction for the Henon map for the
!             unshifted grid, basin. Then, for multiple epsilons, we compare
!             the basins obtained out of shifted grids by a magnitued epsilon.
!             This is done in the subroutine basin_compare. Such subroutine
!             returns then the number of epsilon uncertain grid points uncert.
!             MPI is implemented by calculating the number of columns each
!             processor has to deal with. 
    implicit none
    integer, intent(in)  :: me, npes, K, epsdim
    integer, intent(inout) :: uncert(epsdim)
    real(DP), intent(in) :: xextent(6), params(2), eps(epsdim), LOCKOUT

    integer  :: j, ierr, basindim, N, M, NGRIDx, NGRIDy
    real(DP) :: d
    logical,  dimension(:)  , allocatable :: basin

    NGRIDx = xextent(5) 
    NGRIDy = xextent(6) 
! The following section assigns the number of columns that each processor has to
! work with
    N  = NGRIDx/npes
    M  = MOD(int(NGRIDx),npes)
    if (me.lt.M) then
      allocate(basin((N+1)*NGRIDy))
    else
      allocate(basin(N*NGRIDy))
    endif
! Calculate the unshifted basin
    call henon_map(basin,size(basin),xextent,me,npes,params,K,LOCKOUT,0d0)

    do j=1,epsdim
      call basin_compare(basin,size(basin),xextent,me,npes,params,K,LOCKOUT,eps(j),uncert(j))
    enddo

    return
  end subroutine basin_alg
!-------------------------------------------------------------------------------
end module henondim
