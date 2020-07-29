module substitution
  use precision
  use, intrinsic :: iso_fortran_env
  implicit none
  contains
    subroutine forward_subs(L, n, b, ierr)
        integer,  intent(in)    :: n
        real(DP), intent(in)    :: L(n,n)
        real(DP), intent(inout) :: b(n)
        integer,  intent(out)   :: ierr
        integer :: j 
!        write(OUTPUT_UNIT,*) L
!        write(OUTPUT_UNIT,*) L(1,1), L(2,2), L(3,3)
        do j=1,n-1
          if (L(j,j).eq.0d0) goto 911
          b(j)=b(j)/L(j,j)
          b(j+1:n)=b(j+1:n)-b(j)*L(j+1:n,j)
        end do
        b(n)=b(n)/L(n,n)
        return
911 continue
        ierr = -1
        return
    end subroutine forward_subs

    subroutine backward_subs(U, n, b, ierr)
        integer,  intent(in)    :: n
        real(DP), intent(in)    :: U(n,n)
        real(DP), intent(inout) :: b(n)
        integer,  intent(out)   :: ierr
        integer :: j 

        do j=n,2,-1
          if (U(j,j).eq.0) goto 911
          b(j)=b(j)/U(j,j)
          b(1:j-1)=b(1:j-1)-b(j)*U(1:j-1,j)
        end do
        b(1)=b(1)/U(1,1)
        return
911 continue
        ierr = -1
        return
    end subroutine backward_subs 
end module substitution
