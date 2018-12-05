# HOMEWORK_5 Makefile

OBJs  =  precision.o jacobimpi.o 
mpif90comp = mpiifort
MKLLIB = -mkl # =parallel -liomp5 -lpthread -axavx 
switch = -O3 -qopenmp -heap-arrays# -g -traceback -check arg_temp_created -check all#-qopenmp-report=2# -fpe0 -check all -check bounds -g#-diag-file-append=intel_build.log -fpe0

# Makefile
.phony: clean

jacobimpi: $(OBJs)
	module load intel-mpi/2018x; \
	$(mpif90comp) $(switch) -o jacobimpi.exe $(OBJs) $(MKLLIB) 
%.o: %.f90
	module load intel-mpi/2018x; \
	$(mpif90comp) -c $(switch) $<
%.mod: %.o %.f90
	$(mpif90comp) -c $(switch) $<
clean:
	rm -f *.out
	rm -f *.err
	rm -f *.log
	rm -f *.mod
	rm -f *.o
	rm -f *.exe

# End of Makefile
