program jacobimpi
    use, intrinsic::iso_fortran_env
    use precision

    implicit none
    integer, parameter :: ROOT=1, INTERNAL_ERROR=1
    integer  :: me, npes, beta, k, k1, k2
    integer  :: M, N, i, j, iter, ierr, nargs, auxM, auxN
    real(DP) :: xMin, xMax, yMin, yMax, hx, hy
    real(DP) :: delta[*], deltaALL[*], c[*], c_ALL, bcleft, bcright, bcup, bcdown
    real(DP), dimension(:),   allocatable :: x, y
    real(DP), dimension(:,:), allocatable :: f, sol, u[:], uold, ucomp
    real(DP), dimension(:,:), allocatable :: metest

    real(DP) :: xyextent(4)[*], BCsData(4)[*], deltaConv[*]
    integer  :: gridData(2)[*], iterMax[*], q[*], Mme, Nme
    character(20)  :: filename='jacobi.txt'
    character(200) :: errmsg='none'
    namelist/jacobiparams/ xyextent, gridData, iterMax, q, deltaConv, BCsData

    deltaALL = 0
    c_ALL = 0
    me = this_image() 
    npes = num_images()
    sync all

    if (me==ROOT) then
      nargs = command_argument_count()
      if(nargs.gt.0) call get_command_argument(1,filename)
      open(unit=4, file=filename, status='old', iostat=ierr, iomsg=errmsg)
      if(ierr.ne.0) goto 911
      read(4,jacobiparams,iostat=ierr,iomsg=errmsg)
      if(ierr.ne.0) goto 911
      close(4,iostat=ierr, iomsg=errmsg)
      if(ierr.ne.0) goto 911

      do i=2,npes
        xyextent(:)[i] = xyextent(:)
        gridData(:)[i] = gridData(:)
        iterMax[i]  = iterMax
        q[i] = q
        deltaConv[i]  = deltaConv
        BCsData(:)[i] = BCsData(:)
      enddo
    endif
    sync all

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
      print*, "| Number of images : ", npes
      print*, "| Image number : ", me
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
    k1 = MOD(me-1,beta)
    k2 = FLOOR(DBLE((me-1)/beta))
  
    if (k1.eq.0) then
      allocate(x(0:auxM))
      x = [(hx*(j-1)+xMin, j=1, auxM+1)]
    elseif (k1.eq.beta-1) then
      allocate(x(0:auxM))
      x = [(hx*(j-1)+xMin, j=auxM*k1, auxM*(k1+1))]
    else
      allocate(x(0:auxM+1))
      x = [(hx*(j-1)+xMin, j=auxM*k1, auxM*(k1+1)+1)]
    endif

    if (k2.eq.0) then
      allocate(y(0:auxN))
      y = [(hy*(j-1)+yMin, j=1, auxN+1)]
    elseif (k2.eq.beta-1) then
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
    allocate(u(0:Mme,0:Nme)[*])
    allocate(uold(0:Mme,0:Nme))
    allocate(ucomp(0:Mme,0:Nme))
 
    f   = rhs(x,y,Mme,Nme)
    sol = uexact(x,y,Mme,Nme)

  ! Initialization of u and Boundary Conditions
    u = 0d0  ! Initialization of u
    if (k1.eq.0) u(0,:)=bcleft
    if (k1.eq.beta-1) u(Mme,:)=bcright
    if (k2.eq.0) u(:,0)=bcdown
    if (k2.eq.beta-1) u(:,Nme)=bcup

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

      sync all
      ! Update left ghost cells
      if (k1.ne.0) then ! No left subdomain
        if (k1.eq.1) then ! The left subdomain has lower x-dimension
          u(0,1:Nme-1) = u(Mme-2,1:Nme-1)[me-1]
        elseif (k1.eq.beta-1) then ! The left subdomain has higher x-dimension
          u(0,1:Nme-1) = u(Mme,1:Nme-1)[me-1]
        else ! The left subdomain has same x-dimension
          u(0,1:Nme-1) = u(Mme-1,1:Nme-1)[me-1]
        endif
      endif
      ! Update right ghost cells
      if (k1.ne.beta-1) then ! No right subdomain
        u(Mme,1:Nme-1) = u(1,1:Nme-1)[me+1] ! No need to check dimensions of
                     ! the different subdomains since they all start as 0,1,...
      endif
      ! Update bottom ghost cells
      if (k2.ne.0) then ! No bottom subdomain
        if (k2.eq.1) then ! The bottom subdomain has lower y-dimension
          u(1:Mme-1,0) = u(1:Mme-1,Nme-2)[me-beta]
        elseif (k1.eq.beta-1) then ! The bottom subdomain has higher y-dimension
          u(1:Mme-1,0) = u(1:Mme-1,Nme)[me-beta]
        else ! The bottom subdomain has same y-dimension
          u(1:Mme-1,0) = u(1:Mme-1,Nme-1)[me-beta]
        endif
      endif
      ! Update top ghost cells
      if (k2.ne.beta-1) then ! No top subdomain
        u(1:Mme-1,Nme) = u(1:Mme-1,1)[me+beta] ! No need to check dimensions of
                     ! the different subdomains since they all start as 0,1,...
      endif

      uold = u
      if (MOD(iter,q)==0) then 
        delta = InfNorm(u-ucomp)
        sync all
        if (me.eq.ROOT) then
          do i=1,npes
            deltaALL = MAX(delta[i],deltaALL)
          enddo
          do i=2,npes
            deltaALL[i] = deltaALL
          enddo
        endif
        if(deltaALL.lt.deltaConv) exit 
        if(iter.lt.iterMax) ucomp = u
      endif
    enddo

    c = InfNorm(u-sol)
    sync all

    if (me.eq.ROOT) then
      do i=1,npes
        c_ALL = MAX(c[i],c_ALL)
      enddo

      print*, " "
      print*, "---------"
      print*, "SOLUTION"
      print*, "---------"
      print*, "Iterations: ", iter
      print*, "delta : ", deltaALL
      print*, "c : ", c_ALL 
    endif

    stop

911 continue
    if(me.eq.ROOT) write(ERROR_UNIT,999) me, errmsg
999 format('MPI process ', i0, ': ', a)    
   
    
!------------------------------------------------------------------------------
!---------------------------------- CONTAINS ----------------------------------
!------------------------------------------------------------------------------
    contains
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
