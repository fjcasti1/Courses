# HOMEWORK_5 Makefile

OBJs =  precision.o jacobicoarray.o 
f90comp = ifort
MKLLIB = -mkl # =parallel -liomp5 -lpthread -axavx 
switch = -O3 -heap-arrays -coarray=shared# -g -traceback -check arg_temp_created -check all#-qopenmp-report=2# -fpe0 -check all -check bounds -g#-diag-file-append=intel_build.log -fpe0

# Makefile
.phony: clean

jacobicoarray: $(OBJs)
	module load intel-mpi/2018x; \
	$(f90comp) $(switch) -o jacobicoarray.exe $(OBJs) $(MKLLIB) 
%.o: %.f90
	module load intel-mpi/2018x; \
	$(f90comp) -c $(switch) $<
%.mod: %.o %.f90
	$(f90comp) -c $(switch) $<
clean:
	rm -f *.out
	rm -f *.err
	rm -f *.log
	rm -f *.mod
	rm -f *.o
	rm -f *.exe

# End of Makefile
