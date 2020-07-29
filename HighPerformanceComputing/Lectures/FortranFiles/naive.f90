!  NAIVE - naive method for doing forward substitution.
!  Replace this module with one called SUBSTITUTION that implements
!  a vectorized version of both the forward and backward substitution
!  algorithms.  The vectorized versions should access memory more
!  efficiently than this version.
 
    module naive
    use precision
    implicit none
    contains
!------------------------------------------------------------------------------
    subroutine forward_naive(a, n, b, ierr)
!  FORWARD_NAIVE - naive implementation of forward substitution.
!  We seek to solve a linear system of the form Ax=b, where A is a
!  square lower-triangular matrix.
!
!  This routine is VERY INEFFICIENT for large matrices, because the innermost
!  loop traverses the rows of the matrix A instead of the columns.
!  Arguments:
!  A  : an NxN matrix, assumed lower triangular (the elements of the upper
!       triangle are not accessed).
!  N : the size of the problem
!  B :=: On input, the right-hand side to be solved for; on return, holds
!      the solution, X.
!  IERR := set to 0 unless an error occurs (i.e., the system is exactly
!     singular).

    integer, intent(in):: n
    real(DP), intent(in):: a(n,n)
    real(DP), intent(inout):: b(n)
    integer, intent(out):: ierr

!  Local variables

    integer:: j,k

!  Start with the first equation and forward substitute into the remainder

    ierr = 0
    do k = 1, n
       if(a(k,k).eq.0) goto 911  ! exactly singular problem
       do j = 1, k-1
          b(k) = b(k) - a(k,j)*b(j)
       enddo
       b(k) = b(k)/a(k,k)
    enddo  ! or END DO
    return

!  Oops

911 continue
    ierr = -1
    return
    end subroutine forward_naive
!------------------------------------------------------------------------------
    end module naive
