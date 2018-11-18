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
!  use lsq
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
! ARGUMENTS - basin : logical, 1D array of dimension (NGRID*NGRID), True will
!                     mean that the grid point tends to infinity.
!           - eps   : real, eps represents the magnitude of the shifting
!                     of the initial grid.
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

!    write(OUTPUT_UNIT, 100) me, npes
!    write(OUTPUT_UNIT, 101) me, params
!    write(OUTPUT_UNIT, 102) me, xextent
!    write(OUTPUT_UNIT, 103) me, N
!    write(OUTPUT_UNIT, 104) me, M
!    do j=1,size(x)
!      write(OUTPUT_UNIT, 105) me, j, x(j)
!    end do
!    do j=1,NGRIDy
!      write(OUTPUT_UNIT, 106) me, j, y(j)
!    end do
!    do j=1,size(X0me,1)
!      write(OUTPUT_UNIT, 107) me, j, X0me(j,1), X0me(j,2)
!    end do
!    write(OUTPUT_UNIT, 108) me, hx, hy

    basin=all(abs(X0me)>LOCKOUT,2)

    do j=1,K
      where(basin.eq..False.)
        auxme     = X0me(:,1)
        X0me(:,1) = params(1)-X0me(:,1)**2d0+params(2)*X0me(:,2)
        X0me(:,2) = auxme(:) 
        basin = all(abs(X0me)>LOCKOUT,2)
      end where
    end do

!100 format('process ', i0, ': process count: ', i0)
!101 format('process ', i0, ': parameters: ',2f8.2)
!102 format('process ', i0, ': xextent: ', 6f8.2)
!103 format('process ', i0, ': N: ', i0)
!104 format('process ', i0, ': M: ', i0)
!105 format('process ', i0, ': x(', i0,'): ', f6.2)
!106 format('process ', i0, ': y(', i0,'): ', f6.2)
!107 format('process ', i0, ': X0me(', i0,',:): ', f6.2,' , ', f6.2)
!108 format('process ', i0, ': hx=', f6.2,', hy=',f6.2)



 
!    allocate(x(NGRIDx))
!    allocate(y(NGRIDy))
!    allocate(X0(NGRIDx*NGRIDy,2))
!
!    hx = (xMax-xMin)/(NGRIDx-1)
!    hy = (yMax-yMin)/(NGRIDy-1)
!    x  = [(hx*(j-1)+xMin+eps, j=1,NGRIDx)]
!    y  = [(hy*(j-1)+yMin, j=1,NGRIDy)]
!    N  = NGRIDx*NGRIDy/npes
!    M  = MOD(int(NGRIDx*NGRIDy),npes)
!
!    do j=1,NGRIDx
!      X0(1+(j-1)*NGRIDy:j*NGRIDy,1)=x(j)
!      X0(1+(j-1)*NGRIDy:j*NGRIDy,2)=y(:)
!    end do
!
!    if (me.lt.M) then
!      allocate(X0me(N+1,2))
!      allocate(auxme(N+1))
!      X0me=X0((N+1)*me+1:(N+1)*(me+1),:)
!    else
!      allocate(X0me(N,2))
!      allocate(auxme(N))
!      X0me=X0(N*me+M+1:N*(me+1)+M,:)
!    endif
!
!    basin=all(abs(X0me)>LOCKOUT,2)
!
!    do j=1,K
!      where(basin.eq..False.)
!        auxme     = X0me(:,1)
!        X0me(:,1) = params(1)-X0me(:,1)**2d0+params(2)*X0me(:,2)
!        X0me(:,2) = auxme(:) 
!        basin = all(abs(X0me)>LOCKOUT,2)
!      end where
!    end do
  end subroutine henon_map
!-------------------------------------------------------------------------------
  subroutine basin_compare(basin0,basindim,xextent,me,npes,params,K,LOCKOUT,eps,Nuncert)
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
!    use, intrinsic::iso_fortran_env
    implicit none
    integer,  intent(in)  :: me, npes, K, basindim
    integer,  intent(out) :: Nuncert
    real(DP), intent(in)  :: xextent(6), params(2), eps, LOCKOUT
    logical,  intent(in)  :: basin0(basindim)

    logical, dimension(:), allocatable :: basinP, basinM

    allocate(basinP(basindim))
    allocate(basinM(basindim))

    call henon_map(basinP,basindim,xextent,me,npes,params,K,LOCKOUT,eps)
    call henon_map(basinM,basindim,xextent,me,npes,params,K,LOCKOUT,-eps)
    Nuncert = count(basin0.ne.basinP.or.basin0.ne.basinM)
  end subroutine basin_compare
!-------------------------------------------------------------------------------
  subroutine basin_alg(xextent,me,npes,params,K,LOCKOUT,eps,epsdim,uncert)
     
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
!    use, intrinsic::iso_fortran_env
    implicit none
    integer, intent(in)  :: me, npes, K, epsdim
    integer, intent(out) :: uncert(epsdim)
    real(DP), intent(in) :: xextent(6), params(2), eps(epsdim), LOCKOUT

    integer  :: j, ierr, basindim, N, M, NGRIDx, NGRIDy
    real(DP) :: d
    logical,  dimension(:)  , allocatable :: basin
    

    NGRIDx = xextent(5) 
    NGRIDy = xextent(6) 
 
    N  = NGRIDx*NGRIDy/npes
    M  = MOD(int(NGRIDx*NGRIDy),npes)
    if (me.lt.M) then
      allocate(basin(N+1))
    else
      allocate(basin(N))
    endif

    call henon_map(basin,size(basin),xextent,me,npes,params,K,LOCKOUT,0d0)

    do j=1,epsdim
      call basin_compare(basin,size(basin),xextent,me,npes,params,K,LOCKOUT,eps(j),uncert(j))
    enddo

!    call linfit('Power',Neps,eps,N,param,ierr)    
!    if(ierr==0) d=1+param(2)
    return
  end subroutine basin_alg
    


!    print *, size(eps)
!    print *, K 
!    print *, LOCKOUT 
!    do j=1,size(eps)
!      print*, eps(j)
!    end do

!    write(OUTPUT_UNIT, 100) me, npes
!    write(OUTPUT_UNIT, 101) me, params
!    write(OUTPUT_UNIT, 102) me, xextent
!    write(OUTPUT_UNIT, 103) me, N
!    write(OUTPUT_UNIT, 104) me, M
!    do i=1,NGRIDx
!      write(OUTPUT_UNIT, 105) me, i, x(i)
!      !print *, x(i)
!    end do
!    do j=1,NGRIDy
!      write(OUTPUT_UNIT, 106) me, j, y(j)
!      !print *, y(j)
!    end do
!    do i=1,NGRIDx*NGRIDy
!      write(OUTPUT_UNIT, 107) me, i, X0(i,1)
!    end do
!    do i=1,NGRIDx*NGRIDy
!      write(OUTPUT_UNIT, 108) me, i, X0(i,2)
!    end do

!    do i=1,size(X0me,1)
!      write(OUTPUT_UNIT, 111) me, i, X0me(i,1), i, X0me(i,2)
!    end do
!    write(OUTPUT_UNIT, 112) me, hx, hy

!100 format('process ', i0, ': process count: ', i0)
!101 format('process ', i0, ': parameters: ',2f8.2)
!102 format('process ', i0, ': xextent: ', 6f8.2)
!103 format('process ', i0, ': N: ', i0)
!104 format('process ', i0, ': M: ', i0)
!105 format('process ', i0, ': x(', i0,'): ', f6.2)
!106 format('process ', i0, ': y(', i0,'): ', f6.2)
!107 format('process ', i0, ': X0(', i0,',1): ', f6.2)
!108 format('process ', i0, ': X0(', i0,',2): ', f6.2)
!109 format('process ', i0, ': less than: ', i0)
!110 format('process ', i0, ': greater than or equal to: ', i0)
!111 format('process ', i0, ': X0me(', i0,',1): ', f6.2,', X0me(', i0,',2): ', f6.2)
!112 format('process ', i0, ': hx=', f6.2,', hy=',f6.2)

!  Place dummy values in the first 3 elements of UNCERT

!   uncert=[me, 2*me, 3*me]  ! depends on rank and number of processors
!
!   return
!   end subroutine basin_dim
!-------------------------------------------------------------------------------
end module henondim
