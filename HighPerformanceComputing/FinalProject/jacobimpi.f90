program jacobimpi
    use, intrinsic::iso_fortran_env
    use mpi
    use precision
 !   use jacobimodule
    implicit none
    integer, parameter :: ROOT=0, INTERNAL_ERROR=1
    integer  :: me, npes, beta, k, k1, k2
    integer  :: M, N, i, j, iter, ierr, nargs, auxM, auxN
    real(DP) :: xMin, xMax, yMin, yMax, hx, hy
    real(DP) :: delta, deltaALL, c, bcleft, bcright, bcup, bcdown
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

    if(me.eq.ROOT) then
      print*, " "
      print*, "--------------------------------------"
      print*, "|------------ PARAMETERS ------------|"
      print*, "--------------------------------------"
      print*, "| Number of processes : ", npes
      print*, "| Process number : ", me
      print*, "--------------------------------------"
      print*, "| xMin = ", xMin
      print*, "| xMax = ", xMax
      print*, "| yMin = ", yMin
      print*, "| yMax = ", yMax
      print*, "| M = ", M
      print*, "| N = ", N
      print*, "| iterMax = ", iterMax
      print*, "| q = ", q
      print*, "| deltaConv = ", deltaConv
      print*, "| BCleft  = ", bcleft
      print*, "| BCright = ", bcright
      print*, "| BCup    = ", bcup
      print*, "| BCdown  = ", bcdown
      print*, "--------------------------------------"
    endif

    hx = (xMax-xMin)/M
    hy = (yMax-yMin)/N

!'''''''''''''''''''''''''''''''''''''''''''''''''''''!
!'''''''' IMPLEMENT CHECK FOR PERFECT SQUARE '''''''''!
!'''''''''''''''''''''''''''''''''''''''''''''''''''''!
    beta = sqrt(DBLE(npes))
    auxM = (M+1)/beta
    auxN = (N+1)/beta
    k1 = MOD(me,beta)
    k2 = FLOOR(DBLE(me/beta))
  
    if (k1.eq.0) then
      allocate(x(0:auxM))
      x = [(hx*(j-1)+xMin, j=1, auxM+1)]
    elseif (k1.eq.auxM) then
      allocate(x(0:auxM))
      x = [(hx*(j-1)+xMin, j=auxM*k1, auxM*(k1+1))]
    else
      allocate(x(0:auxM+1))
      x = [(hx*(j-1)+xMin, j=auxM*k1, auxM*(k1+1)+1)]
    endif

    if (k2.eq.0) then
      allocate(y(0:auxN))
      y = [(hy*(j-1)+yMin, j=1, auxN+1)]
    elseif (k2.eq.auxN) then
      allocate(y(0:auxN))
      y = [(hy*(j-1)+yMin, j=auxN*k2, auxN*(k2+1))]
    else
      allocate(y(0:auxN+1))
      y = [(hy*(j-1)+yMin, j=auxN*k2, auxN*(k2+1)+1)]
    endif
    
    Mme = size(x)-1
    Nme = size(y)-1

    allocate(f(0:Mme,0:Nme))
    allocate(sol(0:Mme,0:Nme))
    allocate(u(0:Mme,0:Nme))
    allocate(uold(0:Mme,0:Nme))
    allocate(ucomp(0:Mme,0:Nme))
 
    f   = rhs(x,y,Mme,Nme)
    sol = uexact(x,y,Mme,Nme)

  ! Initialization of u and Boundary Conditions
    u = 0d0  ! Initialization of u
    if (k1.eq.0) u(0,:)=bcleft
    if (k1.eq.(beta-1)) u(Mme,:)=bcright
    if (k2.eq.0) u(:,0)=bcdown
    if (k2.eq.(beta-1)) u(:,Nme)=bcup

!    call jacobiIter(u,f,size(x)-1,size(y)-1,iterMax,delta,deltaConv,q,hx,hy,iter,k2,beta)
    uold  = u
    ucomp = u

    do iter=1,iterMax
!    do iter=1,200
      do i=1,Mme-1
        do j=1,Nme-1
          u(i,j) = (hy**2d0*(uold(i-1,j)+uold(i+1,j))+&
                    hx**2d0*(uold(i,j-1)+uold(i,j+1))-&
                    hx**2d0*hy**2d0*f(i,j))/(2*(hx**2d0+hy**2d0))
        enddo
      enddo

!      print*, "process : ", me, " iter = ", iter
!      call updateGhostCells(me,u,Mme,Nme,k2,beta) 
      if (MOD(me,2)==1) then
        call MPI_Sendrecv(u(1,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,90,&
            u(0,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,91,MPI_COMM_WORLD,status,ierr)
        if(MOD((me+1),beta).ne.0) call MPI_Sendrecv(u(Mme-1,1:Nme-1),Nme-1,&
                MPI_DOUBLE,me+1,92,u(Mme,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,93,&
                                                      MPI_COMM_WORLD,status,ierr)
      else
        call MPI_Sendrecv(u(Mme-1,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,91,&
            u(Mme,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,90,MPI_COMM_WORLD,status,ierr)
        if(MOD(me,beta).ne.0) call MPI_Sendrecv(u(1,1:Nme-1),Nme-1,MPI_DOUBLE,&
                me-1,93,u(0,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,92,MPI_COMM_WORLD,&
                                                                      status,ierr)
      endif

      if (MOD(k2,2)==1) then
        call MPI_Sendrecv(u(1:Mme-1,1),Mme-1,MPI_DOUBLE,me-beta,94,&
            u(1:Mme-1,0),Mme-1,MPI_DOUBLE,me-beta,95,MPI_COMM_WORLD,status,ierr)
        if(k2+1.lt.beta) call MPI_Sendrecv(u(1:Mme-1,Nme-1),Mme-1,MPI_DOUBLE,&
                      me+beta,96,u(1:Mme-1,Nme),Mme-1,MPI_DOUBLE,me+beta,97,&
                                                      MPI_COMM_WORLD,status,ierr)
      else
        call MPI_Sendrecv(u(1:Mme-1,Nme-1),Mme-1,MPI_DOUBLE,me+beta,95,&
                  u(1:Mme-1,Nme),Mme-1,MPI_DOUBLE,me+beta,94,MPI_COMM_WORLD,&
                                                                      status,ierr)
        if(k2.ne.0) call MPI_Sendrecv(u(1:Mme-1,1),Mme-1,MPI_DOUBLE,me-beta,97,&
                  u(1:Mme-1,0),Mme-1,MPI_DOUBLE,me-beta,96,MPI_COMM_WORLD,&
                                                                    status,ierr)
      endif

      uold = u
      if (MOD(iter,q)==0) then 
        k=iter
        delta = InfNorm(u-ucomp)
        call MPI_ALLreduce(delta,deltaALL,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
        if(deltaALL.lt.deltaConv) exit 
        if(iter.lt.iterMax) ucomp = u
      endif
    enddo
    
    call MPI_reduce(InfNorm(u-sol),c,1,MPI_DOUBLE,MPI_MAX,ROOT,MPI_COMM_WORLD,ierr)

    if (me.eq.ROOT) then
      print*, " "
      print*, "---------"
      print*, "SOLUTION"
      print*, "---------"
      print*, "Iterations: ", iter
      print*, "delta : ", deltaALL 
      print*, "c : ", c 
    endif

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
      subroutine jacobiIter(u,f,Mme,Nme,iterMax,delta,deltaConv,q,hx,hy,iter,k2,beta)
        implicit none
        integer,  intent(in)    :: Mme, Nme, iterMax, q, k2, beta
        integer,  intent(out)   :: iter
        real(DP), intent(in)    :: f(0:Mme,0:Nme)
        real(DP), intent(in)    :: deltaConv, hx, hy
        real(DP), intent(inout) :: u(0:Mme,0:Nme), delta
  
        integer  :: i, j, k
        real(DP) :: uold(0:Mme,0:Nme), ucomp(0:Mme,0:Nme)
  
        uold  = u
        ucomp = u
  
        do iter=1,iterMax
          do i=1,Mme-1
            do j=1,Nme-1
              u(i,j) = (hy**2d0*(uold(i-1,j)+uold(i+1,j))+&
                        hx**2d0*(uold(i,j-1)+uold(i,j+1))-&
                        hx**2d0*hy**2d0*f(i,j))/(2*(hx**2d0+hy**2d0))
            enddo
          enddo

!          call updateGhostCells(me,u,Mme,Nme,k2,beta) 

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
      subroutine updateGhostCells(me,u,Mme,Nme,k,beta)
        implicit none
        integer,  intent(in)    :: me, Mme, Nme, k, beta
        real(DP), intent(inout) :: u(0:Mme,0:Nme)

        integer :: ierr

        if (MOD(me,2)==1) then
          call MPI_Sendrecv(u(1,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,90,&
              u(0,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,91,MPI_COMM_WORLD,status,ierr)

          if(MOD((me+1),beta).ne.0) call MPI_Sendrecv(u(Mme-1,1:Nme-1),Nme-1,&
                  MPI_DOUBLE,me+1,92,u(Mme,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,93,&
                                                        MPI_COMM_WORLD,status,ierr)
        else
          call MPI_Sendrecv(u(Mme-1,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,91,&
              u(Mme,1:Nme-1),Nme-1,MPI_DOUBLE,me+1,90,MPI_COMM_WORLD,status,ierr)

          if(MOD(me,beta).ne.0) call MPI_Sendrecv(u(1,1:Nme-1),Nme-1,MPI_DOUBLE,&
                  me-1,93,u(0,1:Nme-1),Nme-1,MPI_DOUBLE,me-1,92,MPI_COMM_WORLD,&
                                                                        status,ierr)
        endif

        if (MOD(k,2)==1) then
          call MPI_Sendrecv(u(1:Mme-1,1),Mme-1,MPI_DOUBLE,me-beta,94,&
              u(1:Mme-1,0),Mme-1,MPI_DOUBLE,me-beta,95,MPI_COMM_WORLD,status,ierr)

          if(k+1.lt.beta) call MPI_Sendrecv(u(1:Mme-1,Nme-1),Mme-1,MPI_DOUBLE,&
                        me+beta,96,u(1:Mme-1,Nme),Mme-1,MPI_DOUBLE,me+beta,97,&
                                                        MPI_COMM_WORLD,status,ierr)
        else
          call MPI_Sendrecv(u(1:Mme-1,Nme-1),Mme-1,MPI_DOUBLE,me+beta,95,&
                    u(1:Mme-1,Nme),Mme-1,MPI_DOUBLE,me+beta,94,MPI_COMM_WORLD,&
                                                                        status,ierr)

          if(k.ne.0) call MPI_Sendrecv(u(1:Mme-1,1),Mme-1,MPI_DOUBLE,me-beta,97,&
                    u(1:Mme-1,0),Mme-1,MPI_DOUBLE,me-beta,96,MPI_COMM_WORLD,&
                                                                      status,ierr)
        endif

      end subroutine updateGhostCells
!------------------------------------------------------------------------------
      function InfNorm(V)
          implicit none
          real(DP), dimension(0:,0:), intent(in) :: V
          real(DP) :: InfNorm

          InfNorm = maxval(maxval(abs(V),1))
      end function InfNorm
!------------------------------------------------------------------------------
      function rhs(x,y,M,N)
        implicit none
        integer , intent(in) :: M, N
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
      function uexact(x,y,M,N)
        implicit none
        integer , intent(in) :: M, N
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





!!!    print*, " "
!!!    print*, "--------------------------------------"
!!!    print*, "|------------ GRID CHECK ------------|"
!!!    print*, "--------------------------------------"
!!!    print*, "| Number of processes : ", npes
!!!    print*, "| Process number : ", me
!!!    print*, "--------------------------------------"
!!!    print*, "| Mme = ", Mme
!!!    print*, "| Nme = ", Nme
!!!    print*, "| "
!!!    print*, "| grid test : "
!!!    do j=0,Nme
!!!      print*, metest(:,j)
!!!    enddo
!!!    print*, "--------------------------------------"
    

!!!    print*, " "
!!!    print*, "--------------------------------------"
!!!    print*, "|----------- GRID CHECK 2 ------------|"
!!!    print*, "--------------------------------------"
!!!    print*, "| Number of processes : ", npes
!!!    print*, "| Process number : ", me
!!!    print*, "| Floor(me/beta) : ", k
!!!    print*, "--------------------------------------"
!!!    print*, "| Mme = ", Mme
!!!    print*, "| Nme = ", Nme
!!!    print*, "| "
!!!    print*, "| grid test : "
!!!    do j=0,Nme
!!!      print*, metest(:,j)
!!!    enddo
!!!    print*, "--------------------------------------"

!!!    print *, " I am MPI process ", me, "of ", npes, ". On my left is ", left, &
!!!                                " and on my right is ", right
!!!    print *, " "

!    print*, " "
!    print*, "--------------------------------------"
!    print*, "|------------ GRID CHECK ------------|"
!    print*, "--------------------------------------"
!    print*, "| Number of processes : ", npes
!    print*, "| Process number : ", me
!    print*, "--------------------------------------"
!    print*, "| Mme = ", Mme
!    print*, "| Nme = ", Nme
!    print*, "| k1 = ", k1
!    print*, "| k2 = ", k2
!    print*, "| size(x) = ", size(x)
!    print*, "| size(y) = ", size(y)
!    print*, "| size(f) = ", size(f,1), ",  ", size(f,2)
!    print*, "| size(sol) = ", size(sol,1), ",  ", size(sol,2)
!    print*, "| "
!    print*, "| vector x : "
!    print*, x(:)
!    print*, "| vector y : "
!    print*, y(:)
!    print*, "--------------------------------------"
