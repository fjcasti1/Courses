program test_henon
  ! This program just tries the code written for the assigment.
  use henondim
  use lsq
  use precision
  implicit none
  real(DP) :: d
  integer :: j,Neps
  real(DP), dimension(:), allocatable  :: eps, N
  
  Neps=10
  allocate(N(Neps),eps(Neps))
  eps = (/(2d0**j, j=-21,-12) /)

  call basin_alg(Neps,eps,d)
  print *, " "
  print *, "d =",d
  
end program test_henon
