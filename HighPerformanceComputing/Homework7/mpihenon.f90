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

!  Local variables

    integer:: j, ierr, nargs, ierrfit
    real(DP),parameter:: HALF=0.5

!  These are the epsilon values requested in previous assignments, and you
!  will use the same values for this program.

    real(DP) :: epsilon(NEPS)=[(HALF**j, j=12,21)]

!  Map parameters, as a 2-vector

    real(DP) :: henonparams(2)  !  parameters a ,b
    integer  :: K        !  max number of iterations
    real(DP) :: LOCKOUT  !  threshold

!  Other variables

!  The first four values determine the X and Y limits of the grid,
!  and the second two values are the number of grid points.  XNUM and YNUM
!  are intended to be integer valued, but they are declared as double
!  precision so that they can be sent as a single message.  (As integer
!  values, XNUM and YNUM would have to be sent either as a separate
!  message, or all these values would have to be wrapped as a user-defined
!  type and registered with the MPI library.)  

    real(DP):: xextent(6)  ! xmin, xmax, ymin, ymax, Nx, Ny

!  Other variables

    character(200):: errmsg='none'  ! sufficiently long for typical messages
    character(80):: filename='henon.txt'  ! default input file name
    integer:: uncert(NEPS)  ! number of uncertain points for each epsilon
    integer:: uncertsum(NEPS)  ! total number of uncertain points
    integer:: me  ! my MPI rank (0 = root process)
    integer:: npes  ! Number of Processing ElementS (> 0)
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

!  Once more, broadcast the max number of iterations K and the lockout parameter LOCKOUT.

    call mpi_bcast(K, 1, MPI_INT, ROOT, MPI_COMM_WORLD, ierr)
    call mpi_bcast(LOCKOUT, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)

!  The basi_alg subroutine will divide up the grid according to the number of 
!  processes and calculate the number of uncertain points for each epsilon
!  on the correspoinding protion of the grid for that processor. It will return
!  an array of uncertain points (same dimension as the array of epsilons).

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

