   module lsq
!  LSQ - Interface module for LAPACK QR routines, and a special-purpose
!  subroutine for simple Y versus X linear regressions.
!
!  SYNOPSIS
!    use lsq
!
!  DESCRIPTION
!    Write a one-sentence description of LINFIT.
!
!  PUBLIC ROUTINES DEFINED
!    linfit - Y versus X linear regression.
!
!  DEPENDENCIES
!    precision - defines KINDs for single- and double-precision floating point.
!
!  REVISION HISTORY
!    10/15/18 - first implementation.
!
!  PROGRAMMER
!    Your name, email@asu.edu
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
!  LINFIT - briefly describe the input and output arguments here.
!
   character(*), intent(in):: model
   integer, intent(in):: n
   real(DP), intent(in):: x(n), y(n)
   real(DP), intent(out):: param(2)
   integer, intent(out):: ierr
!
!  your code here
!
   end subroutine linfit
!-----------------------
   end module lsq
