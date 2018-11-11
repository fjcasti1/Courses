!  MPIHENON - Example main program for the MPI version of the Henon basin
!  uncertainty calculation.  All values are taken from a file named
!  on the command line (or 'henon.txt' otherwise).

    use, intrinsic::iso_fortran_env  ! Fortran 2003 and later
    use mpi
    use henondim
    use lsq
    use precision
    ! use henondim or whatever module contains your actual Henon basin code
    implicit none
    integer, parameter:: ROOT=0
    integer, parameter:: NEPS=10
    integer, parameter:: INTERNAL_ERROR=1  ! shell error code; must be nonzero
!    integer, parameter:: DP=kind(1.0d0)  ! double precision

!  Local variables

    integer:: j, ierr, nargs, ierrfit
    real(DP),parameter:: HALF=0.5

!  These are the epsilon values requested in previous assignments, and you
!  will use the same values for this program.

    real(DP) :: epsilon(NEPS)=[(HALF**j, j=12,21)]

!  Map parameters, as a 2-vector

    real(DP) :: henonparams(2)  !  parameters a ,b
    integer      :: K        !  max number of iterations
    real(DP) :: LOCKOUT  !  threshold

!  Other variables

!  The first four values determine the X and Y limits of the grid,
!  and the second two values are the number of grid points.  XNUM and YNUM
!  are intended to be integer valued, but they are declared as double
!  precision so that they can be sent as a single message.  (As integer
!  values, XNUM and YNUM would have to be sent either as a separate
!  message, or all these values would have to be wrapped as a user-defined
!  type and registered with the MPI library.)  

    real(DP):: xextent(6)  ! xmin, xmax, ymin, ymax, xnum, ynum

!  Other variables

    character(200):: errmsg='none'  ! sufficiently long for typical messages
    character(80):: filename='henon.txt'  ! default input file name
    integer:: uncert(NEPS)  ! number of uncertain points for each epsilon
    integer:: uncertsum(NEPS)  ! total number of uncertain points
    integer:: me  ! my MPI rank (0 = root process)
    integer:: npes  ! Number of Processing ElementS (> 0), i.e., the
    real(DP) :: d, C(2)
    namelist/basinparams/ henonparams, xextent, K, LOCKOUT

!  Initialize.

    call mpi_init(ierr)  ! should be first executable statement
    call mpi_comm_rank(MPI_COMM_WORLD, me, ierr)  ! my rank
    call mpi_comm_size(MPI_COMM_WORLD, npes, ierr) ! number of processes

!  Start the Henon basin computation.

    uncert=0

!  Read data on the root process and broadcast.  We MUST NOT allow the Fortran
!  runtime to kill the root process if the input data file cannot be opened or
!  if there is an i/o error.  Otherwise, the root process exits without
!  calling MPI_FINALIZE, the other processes hang while waiting from data 
!  from the now-dead root process.
!  General rule: !  ALL i/o statements in an MPI Fortran program
!  should include the IOSTAT= and IOMSG= specifiers, and an appropriate
!  branch taken to display the error message and kill the entire MPI program.
!  Similarly, ALL exceptions in a C++ program must be caught in a global
!  catch{} block that takes analogous action.

    if(me.eq.ROOT) then
       nargs=command_argument_count()
       if(nargs.gt.0) call get_command_argument(1,filename)
       open(unit=4, file=filename, status='old', iostat=ierr, iomsg=errmsg)
       if(ierr.ne.0) goto 911
       read(4,basinparams,iostat=ierr,iomsg=errmsg)
       if(ierr.ne.0) goto 911
       close(4,iostat=ierr,iomsg=errmsg)  ! close the data file
       if(ierr.ne.0) goto 911
    endif

!  Broadcast the map parameters to everyone.
!  All processes must call MPI_BCAST to receive the broadcast messages.

    call mpi_bcast(henonparams, 2, MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)

!  Likewise, broadcast the grid extent and number of grid points in 
!  each direction (as floating-point numbers).

    call mpi_bcast(xextent, size(xextent), MPI_DOUBLE, ROOT, MPI_COMM_WORLD, &
      ierr)

    call mpi_bcast(K, 1, MPI_INT, ROOT, MPI_COMM_WORLD, ierr)

    call mpi_bcast(LOCKOUT, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)

!  Divide up the grid according to the number of processes and calculate
!  the number of uncertain points for each epsilon on our portion of the grid.
!  Include whatever other parameters that you need in this call.
!  You should include appropriate code to handle unexpected errors by
!  writing a message to ERROR_UNIT and calling MPI_ABORT.
    
!    call basin_dim(xextent, me, npes, henonparams, K, LOCKOUT, &
!                                          epsilon, size(epsilon), uncert)

    call basin_alg(xextent, me, npes, henonparams, K, LOCKOUT, &
                                          epsilon, size(epsilon), uncert)

!  Sum up the total number of uncertain points from each process for each
!  epsilon and place the results on the root.  MPI computes a global sum
!  by adding corresponding elements of UNCERT from each process.
!  Other global operations include MPI_PROD (product), MPI_MAX (maximum of each
!  element), and MPI_MIN (minimum of each element).  It's also possible to
!  define your own, but we won't cover that topic here.

    call mpi_reduce(uncert, uncertsum, NEPS, MPI_INTEGER, MPI_SUM, ROOT, &
      MPI_COMM_WORLD, ierr)

!  Let the root process report the results on the standard output.
!  UNCERTSUM is undefined on the other processes.

    if(me.eq.ROOT) then
       do j=1, NEPS
         write(OUTPUT_UNIT,110,iostat=ierr,iomsg=errmsg) epsilon(j),uncertsum(j)
       enddo
       call linfit('Power',size(epsilon),epsilon,DBLE(uncertsum),C,ierrfit)
       if(ierrfit==0) then
         d=1+C(2)
         write(OUTPUT_UNIT,111,iostat=ierr,iomsg=errmsg) d
       endif
    endif
110 format(es13.6, 2x, i7)
111 format('d = ', f10.8)

!  Clean up and exit.

    call mpi_barrier(MPI_COMM_WORLD, ierr)  ! wait for everyone to finish
    call mpi_finalize(ierr)  ! should be last executable statement before STOP
    stop
!---------   End of normal computation

!  Exception handling:
!  Abort on i/o errors, because if we can't read the data on the root process, 
!  then we can't compute a useful result.
!  Advice: when composing error messages, it is helpful to include the rank
!  of the MPI process, so you will know who is complaining.
!  Always write error messages to the standard error unit, not to the standard
!  output, as MPI launchers usually create separate error log files.
!  In such cases, your error messages won't get mixed in with other
!  computational output.  Error log files should be empty upon successful
!  completion or (if you desire) contain only debugging messages.

911 continue
    if(me.eq.ROOT) write(ERROR_UNIT,999) me, errmsg
999 format('MPI process ', i0, ': ', a)
    call mpi_abort(MPI_COMM_WORLD, INTERNAL_ERROR, ierr)
    end

!---------CUT HERE---------------------------------------------------------
!!!   subroutine basin_dim(xextent, me, npes, params, K, LOCKOUT, eps, epsdim, uncert)
!!!!    call basin_alg(xextent, me, npes, henonparams, K, LOCKOUT, &
!!!!                                          epsilon, size(epsilon), uncert)
!!!!  BASIN_DIM - dummy subroutine that should be deleted here and replaced with
!!!!  a suitably modified routine in module henondim.
!!!!  This version simply prints each process's copy of the inputs to
!!!!  the standard output.  The outputs are likely to be jumbled when
!!!!  run on multiple processors.
!!!
!!!    use, intrinsic::iso_fortran_env
!!!    implicit none
!!!    integer, parameter::DP=kind(1.0d0)
!!!    real(DP), intent(in)::xextent(6)
!!!    real(DP),intent(in):: params(2)
!!!    integer, intent(in)::me, npes
!!!    integer, intent(out)::uncert(3)
!!!    real(DP) :: hx, hy, LOCKOUT
!!!    real(DP), dimension(:)  , allocatable :: x, y, auxm, auxme
!!!    real(DP), dimension(:,:), allocatable :: x0, X0me
!!!    logical, dimension(:), allocatable :: basin
!!!    integer :: i, j, N, M, NGRIDx, NGRIDy, K, epsdim
!!!    real(DP) :: eps(epsdim) 
!!!    
!!!    NGRIDx = xextent(5) 
!!!    NGRIDy = xextent(6) 
!!! 
!!!    allocate(x(NGRIDx))
!!!    allocate(y(NGRIDy))
!!!    allocate(X0(NGRIDx*NGRIDy,2))
!!!
!!!    hx = (xextent(2)-xextent(1))/(NGRIDx-1)
!!!    hy = (xextent(4)-xextent(3))/(NGRIDy-1)
!!!    x  = [(hx*(j-1)+xextent(1), j=1,NGRIDx)]
!!!    y  = [(hy*(j-1)+xextent(3), j=1,NGRIDy)]
!!!    N  = NGRIDx*NGRIDy/npes
!!!    M  = MOD(int(NGRIDx*NGRIDy),npes)
!!!
!!!    do j=1,NGRIDx
!!!      X0(1+(j-1)*NGRIDy:j*NGRIDy,1)=x(j)
!!!      X0(1+(j-1)*NGRIDy:j*NGRIDy,2)=y(:)
!!!    end do
!!!
!!!    if (me.lt.M) then
!!!      allocate(X0me(N+1,2))
!!!      allocate(auxme(N+1))
!!!      allocate(basin(N+1))
!!!      X0me=X0((N+1)*me+1:(N+1)*(me+1),:)
!!!    else
!!!      allocate(X0me(N,2))
!!!      allocate(auxme(N))
!!!      allocate(basin(N))
!!!      X0me=X0(N*me+M+1:N*(me+1)+M,:)
!!!    endif
!!!
!!!!    print *, size(eps)
!!!!    print *, K 
!!!!    print *, LOCKOUT 
!!!!    do j=1,size(eps)
!!!!      print*, eps(j)
!!!!    end do
!!!    basin=all(abs(X0me)>LOCKOUT,2)
!!!
!!!!    do j=1,K
!!!!      where(basin.eq..False.)
!!!!        auxme = X0me(:,1)
!!!!        X0me(:,1) = params(1)-X0me(:,1)**2d0+params(2)*X0me(:,2)
!!!!        X0me(:,2) = auxme(:) 
!!!!        basin = all(abs(X0me)>LOCKOUT,2)
!!!!      end where
!!!!    end do
!!!
!!!!    write(OUTPUT_UNIT, 100) me, npes
!!!!    write(OUTPUT_UNIT, 101) me, params
!!!!    write(OUTPUT_UNIT, 102) me, xextent
!!!!    write(OUTPUT_UNIT, 103) me, N
!!!!    write(OUTPUT_UNIT, 104) me, M
!!!!    do i=1,NGRIDx
!!!!      write(OUTPUT_UNIT, 105) me, i, x(i)
!!!!      !print *, x(i)
!!!!    end do
!!!!    do j=1,NGRIDy
!!!!      write(OUTPUT_UNIT, 106) me, j, y(j)
!!!!      !print *, y(j)
!!!!    end do
!!!!    do i=1,NGRIDx*NGRIDy
!!!!      write(OUTPUT_UNIT, 107) me, i, X0(i,1)
!!!!    end do
!!!!    do i=1,NGRIDx*NGRIDy
!!!!      write(OUTPUT_UNIT, 108) me, i, X0(i,2)
!!!!    end do
!!!
!!!!    do i=1,size(X0me,1)
!!!!      write(OUTPUT_UNIT, 111) me, i, X0me(i,1), i, X0me(i,2)
!!!!    end do
!!!!    write(OUTPUT_UNIT, 112) me, hx, hy
!!!
!!!100 format('process ', i0, ': process count: ', i0)
!!!101 format('process ', i0, ': parameters: ',2f8.2)
!!!102 format('process ', i0, ': xextent: ', 6f8.2)
!!!103 format('process ', i0, ': N: ', i0)
!!!104 format('process ', i0, ': M: ', i0)
!!!105 format('process ', i0, ': x(', i0,'): ', f6.2)
!!!106 format('process ', i0, ': y(', i0,'): ', f6.2)
!!!107 format('process ', i0, ': X0(', i0,',1): ', f6.2)
!!!108 format('process ', i0, ': X0(', i0,',2): ', f6.2)
!!!109 format('process ', i0, ': less than: ', i0)
!!!110 format('process ', i0, ': greater than or equal to: ', i0)
!!!111 format('process ', i0, ': X0me(', i0,',1): ', f6.2,', X0me(', i0,',2): ', f6.2)
!!!112 format('process ', i0, ': hx=', f6.2,', hy=',f6.2)
!!!
!!!!  Place dummy values in the first 3 elements of UNCERT
!!!
!!!   uncert=[me, 2*me, 3*me]  ! depends on rank and number of processors
!!!
!!!   return
!!!   end subroutine basin_dim
