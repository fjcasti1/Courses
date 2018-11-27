program jacobimpi
    use, intrinsic::iso_fortran_env
    use mpi
    use precision
 !   use jacobimodule
    implicit none
    integer :: k
    integer, parameter :: ROOT=0, INTERNAL_ERROR=1
    integer  :: me, npes, beta, k
    integer  :: M, N, i, j, iter, ierr, nargs
    real(DP) :: xMin, xMax, yMin, yMax, hx, hy
    real(DP) :: delta, bcleft, bcright, bcup, bcdown
    real(DP), dimension(:),   allocatable :: x, y
    real(DP), dimension(:,:), allocatable :: f, sol, u, uold, ucomp
    real(DP), dimension(:,:), allocatable :: metest

    real(DP) :: xyextent(4), BCsData(4), deltaConv
    integer  :: gridData(2), iterMax, q, Mme, Nme
    integer  :: status(MPI_STATUS_SIZE)
    character(20)  :: filename='jacobi.txt'
    character(200) :: errmsg='none'
    namelist/jacobiparams/ xyextent, gridData, iterMax, q, deltaConv, BCsData

  
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, npes, ierr)

    if(me.eq.ROOT) then
      nargs = command_argument_count()
      if(nargs.gt.0) call get_command_argument(1,filename)
      open(unit=4, file=filename, status='old', iostat=ierr, iomsg=errmsg)
      if(ierr.ne.0) goto 911
      read(4,jacobiparams,iostat=ierr,iomsg=errmsg)
      if(ierr.ne.0) goto 911
      close(4,iostat=ierr, iomsg=errmsg)
      if(ierr.ne.0) goto 911
    endif

    call mpi_bcast(xyextent,4, MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)
    call mpi_bcast(gridData,2, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
    call mpi_bcast(iterMax,1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
    call mpi_bcast(q,1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
    call mpi_bcast(deltaConv,1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)
    call mpi_bcast(BCsData,4, MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)
  
    xMin = xyextent(1)
    xMax = xyextent(2)
    yMin = xyextent(3)
    yMax = xyextent(4)
    M = gridData(1)
    N = gridData(2)
    bcleft  = BCsData(1)
    bcright = BCsData(2)
    bcup    = BCsData(3)
    bcdown  = BCsData(4)

!!!    print*, " "
!!!    print*, "--------------------------------------"
!!!    print*, "|------------ PARAMETERS ------------|"
!!!    print*, "--------------------------------------"
!!!    print*, "| Number of processes : ", npes
!!!    print*, "| Process number : ", me
!!!    print*, "--------------------------------------"
!!!    print*, "| xMin = ", xMin
!!!    print*, "| xMax = ", xMax
!!!    print*, "| yMin = ", yMin
!!!    print*, "| yMax = ", yMax
!!!    print*, "| M = ", M
!!!    print*, "| N = ", N
!!!    print*, "| iterMax = ", iterMax
!!!    print*, "| q = ", q
!!!    print*, "| deltaConv = ", deltaConv
!!!    print*, "| BCleft  = ", bcleft
!!!    print*, "| BCright = ", bcright
!!!    print*, "| BCup    = ", bcup
!!!    print*, "| BCdown  = ", bcdown
!!!    print*, "--------------------------------------"

    call mpi_barrier(MPI_COMM_WORLD,ierr)

!'''''''''''''''''''''''''''''''''''''''''''''''''''''!
!'''''''' IMPLEMENT CHECK FOR PERFECT SQUARE '''''''''!
!'''''''''''''''''''''''''''''''''''''''''''''''''''''!
    beta = sqrt(DBLE(npes))
    Mme = (M+1)/beta+1
    Nme = (N+1)/beta+1
    allocate(x(0:Mme))
    allocate(y(0:Nme))
    allocate(metest(0:Mme,0:Nme))

    metest = me

    print*, " "
    print*, "--------------------------------------"
    print*, "|------------ GRID CHECK ------------|"
    print*, "--------------------------------------"
    print*, "| Number of processes : ", npes
    print*, "| Process number : ", me
    print*, "--------------------------------------"
    print*, "| Mme = ", Mme
    print*, "| Nme = ", Nme
    print*, "| "
    print*, "| grid test : "
    do j=0,Nme
      print*, metest(:,j)
    enddo
    print*, "--------------------------------------"
    
    if (MOD(me,2)==1) then
      call MPI_Sendrecv(metest(1,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,90,&
          metest(0,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,91,MPI_COMM_WORLD,status,ierr)

      if(MOD((me+1),beta).ne.0) call MPI_Sendrecv(metest(Mme-1,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,92,&
          metest(Mme,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,93,MPI_COMM_WORLD,status,ierr)
    else
      call MPI_Sendrecv(metest(Mme-1,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,91,&
          metest(Mme,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,90,MPI_COMM_WORLD,status,ierr)

      if(MOD(me,beta).ne.0) call MPI_Sendrecv(metest(1,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,93,&
          metest(0,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,92,MPI_COMM_WORLD,status,ierr)
    endif

    k = FLOOR(DBLE(me/beta))

    if (MOD(k,2)==1) then
      call MPI_Sendrecv(metest(1:Mme-1,1),Mme-1,MPI_DOUBLE,me-beta,94,&
              metest(1:Mme-1,0),Mme-1,MPI_DOUBLE,me-beta,95,MPI_COMM_WORLD,status,ierr)

      if(k+1.lt.beta) call MPI_Sendrecv(metest(1:Mme-1,Nme-1),Mme-1,MPI_DOUBLE,me+beta,96,&
            metest(1:Mme-1,Nme),Mme-1,MPI_DOUBLE,me+beta,97,MPI_COMM_WORLD,status,ierr)
    else
      call MPI_Sendrecv(metest(1:Mme-1,Nme-1),Mme-1,MPI_DOUBLE,me+beta,95,&
            metest(1:Mme-1,Nme),Mme-1,MPI_DOUBLE,me+beta,94,MPI_COMM_WORLD,status,ierr)

      if(k.ne.0) call MPI_Sendrecv(metest(1:Mme-1,1),Mme-1,MPI_DOUBLE,me-beta,97,&
              metest(1:Mme-1,0),Mme-1,MPI_DOUBLE,me-beta,96,MPI_COMM_WORLD,status,ierr)
    endif

    call mpi_barrier(MPI_COMM_WORLD,ierr)

    print*, " "
    print*, "--------------------------------------"
    print*, "|----------- GRID CHECK 2 ------------|"
    print*, "--------------------------------------"
    print*, "| Number of processes : ", npes
    print*, "| Process number : ", me
    print*, "| Floor(me/beta) : ", k
    print*, "--------------------------------------"
    print*, "| Mme = ", Mme
    print*, "| Nme = ", Nme
    print*, "| "
    print*, "| grid test : "
    do j=0,Nme
      print*, metest(:,j)
    enddo
    print*, "--------------------------------------"

!!    print *, " I am MPI process ", me, "of ", npes, ". On my left is ", left, &
!!                                " and on my right is ", right
!!    print *, " "

!!!!    hx = (xMax-xMin)/M
!!!!    hy = (yMax-yMin)/N
!!!!    
!!!!    allocate(x(0:M))
!!!!    allocate(y(0:N))
!!!!    allocate(f(0:M,0:N))
!!!!    allocate(sol(0:M,0:N))
!!!!    allocate(u(0:M,0:N))
!!!!    allocate(uold(0:M,0:N))
!!!!  
!!!!    x = [(hx*(j-1)+xMin, j=1, M+1)]
!!!!    y = [(hy*(j-1)+yMin, j=1, N+1)]
!!!!  
!!!!    f   = rhs(x,y)
!!!!    sol = uexact(x,y)
!!!!
!!!!    uold(1:M-1,1:N-1) = 0
!!!!    uold(0,:) = g
!!!!    uold(M,:) = g
!!!!    uold(:,0) = g
!!!!    uold(:,N) = g
!!!!    ucomp = uold
!!!!  
!!!!!!!    print*, " "
!!!!!!!    print*, "Vector x:"
!!!!!!!    print*, " "
!!!!!!!    do i=0,M
!!!!!!!      print*, x(i)
!!!!!!!    enddo
!!!!!!!    print*, " "
!!!!!!!    print*, "Vector y:"
!!!!!!!    print*, " "
!!!!!!!    do i=0,N
!!!!!!!      print*, y(i)
!!!!!!!    enddo
!!!!!!!    print*, " "
!!!!!!!    print*, "RHS:"
!!!!!!!    print*, " "
!!!!!!!    do j=0,N
!!!!!!!      print*, f(:,j)  
!!!!!!!    enddo
!!!!!!!    print*, " "
!!!!!!!    print*, "EXACT SOL:"
!!!!!!!    print*, " "
!!!!!!!    do j=0,N
!!!!!!!      print*, sol(:,j)  
!!!!!!!    enddo
!!!!!!!    print*, " "
!!!!

!!!!    
!!!!    call jacobiIter(u,f,sol,M,N,iterMax,delta,deltaConv,q,hx,hy,iter)
!!!!
!!!!!!!    print*, " "
!!!!!!!    print*, "NUMERICAL SOLUTION"
!!!!!!!    print*, " "
!!!!!!!    print*, "k: ", k
!!!!!!!    print*, "ucomp:"
!!!!!!!    do j=0,N
!!!!!!!      print*, ucomp(:,j)  
!!!!!!!    enddo
!!!!!!!    print*, " "
!!!!!!!    print*, "u:"
!!!!!!!    do j=0,N
!!!!!!!      print*, u(:,j)  
!!!!!!!    enddo
!!!!!!!    print*, " "
!!!!!!!    print*, "Iterations: ", iter
!!!!!!!    print*, " "
!!!!!!!    print*, "delta : ", InfNorm(u-ucomp)
!!!!!!!    print*, " "
!!!!!!!    print*, "c : ", InfNorm(u-sol)
!!!!
!!!!    print*, " "
!!!!    print*, "---------"
!!!!    print*, "SOLUTION"
!!!!    print*, "---------"
!!!!    print*, "Iterations: ", iter
!!!!    print*, "delta : ", InfNorm(u-ucomp)
!!!!    print*, "c : ", InfNorm(u-sol)

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call mpi_finalize(ierr)
    stop


911 continue
    if(me.eq.ROOT) write(ERROR_UNIT,999) me, errmsg
999 format('MPI process ', i0, ': ', a)    
    call mpi_abort(MPI_COMM_WORLD, INTERNAL_ERROR, ierr)
    
!------------------------------------------------------------------------------
!---------------------------------- CONTAINS ----------------------------------
!------------------------------------------------------------------------------
    contains
!------------------------------------------------------------------------------
      subroutine jacobiIter(u,f,sol,M,N,iterMax,delta,deltaConv,q,hx,hy,iter)
        implicit none
        integer,  intent(in)    :: M, N, iterMax, q
        integer,  intent(out)   :: iter
        real(DP), intent(in)    :: f(0:M,0:N), sol(0:M,0:N)
        real(DP), intent(in)    :: deltaConv, hx, hy
        real(DP), intent(inout) :: u(0:M,0:N), delta
  
        integer  :: i, j, k
        real(DP) :: uold(0:M,0:N), ucomp(0:M,0:N)
  
        uold  = u
        ucomp = u
  
        do iter=1,iterMax
          do i=1,M-1
            do j=1,N-1
              u(i,j) = 0.25d0*(uold(i+1,j)+uold(i-1,j)+uold(i,j+1)+uold(i,j-1)&
                                                                  -hx*hy*f(i,j))
            enddo
          enddo
          uold = u
          if (MOD(iter,q)==0) then 
            k=iter
            delta = InfNorm(u-ucomp)
    !!!        print*, " "
    !!!        print*, "TESTING:"
    !!!        print*, "k = ", k
    !!!        print*, "delta = ", delta 
            if(delta.lt.deltaConv) exit 
            if(iter.lt.iterMax) ucomp = u
          endif
        enddo
      end subroutine jacobiIter
!------------------------------------------------------------------------------
      function InfNorm(V)
          implicit none
          real(DP), dimension(0:,0:), intent(in) :: V
          real(DP) :: InfNorm
          InfNorm = maxval(maxval(abs(V),1))
      end function InfNorm
      function rhs(x,y)
        implicit none
        real(DP), dimension(0:M), intent(in) :: x
        real(DP), dimension(0:N), intent(in) :: y
        real(DP), dimension(0:M,0:N) :: rhs
  
        integer :: i, j
  
        do i=0,M
          do j=0,N
            rhs(i,j) = exp(x(i)+y(j))*( (x(i)**2d0+3*x(i))*(y(j)**2d0-y(j)) + &
                                            (y(j)**2d0+3*y(j))*(x(i)**2d0-x(i)) )
          enddo
        enddo
      end function rhs
  !------------------------------------------------------------------------------
      function uexact(x,y)
        implicit none
        real(DP), dimension(0:M), intent(in) :: x
        real(DP), dimension(0:N), intent(in) :: y
        real(DP), dimension(0:M,0:N) :: uexact
  
        integer :: i, j
  
        do i=0,M
          do j=0,N
            uexact(i,j)=exp(x(i)+y(j))*(x(i)**2d0-x(i))*(y(j)**2d0-y(j))
          enddo
        enddo
      end function uexact
  !------------------------------------------------------------------------------

end program jacobimpi
!!!    if (MOD(me,2)==1) then
!!!      call MPI_Sendrecv(me,1,MPI_INTEGER,me-1,90,&
!!!          left,1,MPI_INTEGER,me-1,91,MPI_COMM_WORLD,status,ierr)
!!!
!!!      if (me<npes-1) call MPI_Sendrecv(me,1,MPI_INTEGER,me+1,92,&
!!!          right,1,MPI_INTEGER,me+1,93,MPI_COMM_WORLD,status,ierr)
!!!    else
!!!      call MPI_Sendrecv(me,1,MPI_INTEGER,me+1,91,&
!!!          right,1,MPI_INTEGER,me+1,90,MPI_COMM_WORLD,status,ierr)
!!!
!!!      if (me>0) call MPI_Sendrecv(me,1,MPI_INTEGER,me-1,93,&
!!!          left,1,MPI_INTEGER,me-1,92,MPI_COMM_WORLD,status,ierr)
!!!    endif
