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
!   Your name, email@asu.edu
!

    use precision
    implicit none
    real(DOUBLE), parameter:: LOCKOUT=100
    real(DOUBLE), parameter:: BOXMIN=-3.0, BOXXMAX=3.0  ! box limits
    integer, parameter:: NGRID=4096  ! number of points on a side
    contains
!-------------------------------------------------------------------------------
    subroutine henon_map(argument list here)
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
    return
    end subroutine henon_map
!-------------------------------------------------------------------------------
    subroutine basin_alg(argument list here)
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
    return
    end subroutine basin_alg
!-------------------------------------------------------------------------------
    end module basindim

