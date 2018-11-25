module jacobimodule
    use precision
    use, intrinsic::iso_fortran_env
    contains
!------------------------------------------------------------------------------
      function rhs(x,y,M,N)
        implicit none
        integer :: M, N, i, j
        real(DP), dimension(0:M), intent(in) :: x
        real(DP), dimension(0:N), intent(in) :: y
        real(DP), dimension(0:M,0:N) :: rhs
  
        do i=0,M
          do j=0,N
            rhs(i,j) = exp(x(i)+y(j))*( (x(i)**2d0+3*x(i))*(y(j)**2d0-y(j)) + &
                                            (y(j)**2d0+3*y(j))*(x(i)**2d0-x(i)) )
          enddo
        enddo
      end function rhs
  !------------------------------------------------------------------------------
      function uexact(x,y,M,N)
        implicit none
        integer :: M, N, i, j
        real(DP), dimension(0:M), intent(in) :: x
        real(DP), dimension(0:N), intent(in) :: y
        real(DP), dimension(0:M,0:N) :: uexact
  
        do i=0,M
          do j=0,N
            uexact(i,j)=exp(x(i)+y(j))*(x(i)**2d0-x(i))*(y(j)**2d0-y(j))
          enddo
        enddo
      end function uexact
  !------------------------------------------------------------------------------

end module jacobimodule
