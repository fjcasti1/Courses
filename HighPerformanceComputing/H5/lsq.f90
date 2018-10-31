module lsq
!  LSQ - Interface module for LAPACK QR routines, and a special-purpose
!  subroutine for simple Y versus X linear regressions.
!
!  SYNOPSIS
!    use lsq
!
!  DESCRIPTION
!   The subroutine linfit will calculate the parameters a and b for: 
!     -- Linear fitting: y=a·x+b.
!     -- Power fitting: y=a·x^b, using logarithms to linearize.
!   The code chooses which one to apply based on the given entry "model".
!
!  PUBLIC ROUTINES DEFINED
!    linfit - Y versus X linear regression.
!
!  DEPENDENCIES
!    precision - defines KINDs for single- and double-precision floating point.
!
!  REVISION HISTORY
!    10/20/18 - Final version, including both fittings.
!
!  PROGRAMMER
!    Francisco Castillo-Carrasco, fjcasti1@asu.edu
!
  use precision
  implicit none
!
!  Generic interface (for single- or double-precision) for LAPACK
!  QR factorizations SGELS and DGELS.
!
  interface xgels
    subroutine sgels(tr,n,k,nrhs,x,ldx,y,ldy,work,lwork,info)
        import
        character, intent(in):: tr
        integer, intent(in):: n, k, ldx, ldy, nrhs
        integer, intent(inout):: lwork
        real(SP), intent(inout):: x(ldx,k), y(ldy,nrhs), work(lwork)
        integer, intent(out):: info
     end subroutine sgels
!
     subroutine dgels(tr, n, k, nrhs, x, ldx, y, ldy, work, lwork, info)
        import
        character, intent(in):: tr
        integer, intent(in):: n, k, ldx, ldy, nrhs
        integer, intent(inout):: lwork
        real(DP), intent(inout):: x(ldx,k), y(ldy,nrhs), work(lwork)
        integer, intent(out):: info
     end subroutine dgels
  end interface
!
  contains
!----------
  subroutine linfit(model, n, x, y, param, ierr)
      !  ARGUMENTS - model : character variable, expecting "Linear" or "Power", not
      !                      case sensitive.
      !            - n     : dimension of the arrays x and y.
      !            - x     : real 1D array of dimension n, independent variable.
      !            - y     : real 1D array of dimension n, y=f(x).
      !            - param : real 1D array of dimention 2, contains the parameters
      !                      of the fitting.
      !            - ierr  : error flag, 0 if no errors found.
      !
      ! Comments: When using the LAPACK libraries, the code uses it once with
      ! lwork=-1 to find the optimal lwork and a second time to find the with
      ! the optimal lwork to solve the system.
!
  character(*), intent(in) :: model
  integer , intent(in)  :: n
  real(DP), intent(in)  :: x(n)
  real(DP), intent(inout):: y(n)
  real(DP), intent(out) :: param(2)
  integer , intent(out) :: ierr
  integer  :: lwork
  integer  :: info
  real(DP), dimension(:), allocatable :: work
  real(DP) :: A(n,2), B(n)
  lwork=-1

  select case(model(1:1))
  case('L','l')
    A(:,1)=1
    A(:,2)=x
    allocate(work(1))
    call dgels('N', n, 2, 1, A, n, y, n, work, lwork, info)
    lwork=work(1)
    deallocate(work)

    allocate(work(lwork))
    call dgels('N', n, 2, 1, A, n, y, n, work, lwork, info)
    param(1:2)=y(1:2)
    ierr=info
  case('P','p')
    A(:,1)=1
    A(:,2)=log(x)
    B=log(y)
    allocate(work(1))
    call dgels('N', n, 2, 1, A, n, B, n, work, lwork, info)
    lwork=work(1)
    deallocate(work)

    allocate(work(lwork))
    call dgels('N', n, 2, 1, A, n, B, n, work, lwork, info)
    param(1)=exp(B(1))
    param(2)=B(2)
    ierr=info
  case default
    ierr=-1
    print*, "ERROR, MODEL UNKNOWN"
  end select
!
  end subroutine linfit
!-----------------------
end module lsq
