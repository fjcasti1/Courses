program test_linfit
  use henondim
  use lsq
  use precision
  implicit none
  real :: t0, tf
  real(DP):: x(20), y(20), noise(20), param(2)
  integer :: j, ierr, Neps
!  integer , dimension(:), allocatable  :: N
  real(DP), dimension(:), allocatable  :: eps, N
  logical, dimension(:,:), allocatable :: basin
  
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
  Neps=10
  allocate(basin(NGRID,NGRID))
  allocate(N(Neps),eps(Neps))
  eps = (/(2d0**j, j=-21,-12) /)
!  print *, eps
 !>> Construct arrays N and eps  
  call cpu_time(t0)
  call basin_alg(Neps,N,eps)
  call cpu_time(tf)
  print *, " "
  print *, "Time = ",tf-t0," seconds."

  
end program test_linfit
