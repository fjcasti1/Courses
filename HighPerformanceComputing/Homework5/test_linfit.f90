program test_linfit
  use henondim
  use myStructModule
  use lsq
  use precision
  implicit none
  real(DP):: x(20), y(20), noise(20), param(2)
  integer:: j,ierr
  type(paramStruct) :: params
  
  params%K=100
  params%L=100
  params%a=2.12
  params%b=-0.3
  call random_number(noise) ! uniformly distributed in [0, 1]
  do j=1,20
  x(j) = j
  y(j) = 1.0 + 2*x(j) + 0.05*(noise(j) - 0.5)
  enddo
  
  call linfit('Linear',20,x,y,param,ierr)
!  if(ierr == 0) write(6,*) param
  

  call random_number(noise) ! uniformly distributed in [0, 1]
  do j=1,20
  x(j) = j
  y(j) = 2*x(j)**3 + 0.05*(noise(j) - 0.5)
  enddo

  call linfit('Power',20,x,y,param,ierr)
!  if(ierr == 0) write(6,*) param

  call henon_map(5,-3d0,3d0,-3d0,3d0,params)
end program test_linfit
