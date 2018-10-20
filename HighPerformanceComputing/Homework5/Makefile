# HOMEWORK_5 Makefile

OBJs =  precision.o lsq.o myStructModule.o matlabfns.o henondim.o test_linfit.o
f90comp = ifort
MKLLIB = -mkl # =parallel -liomp5 -lpthread -axavx 
switch = -O3  # -fpe0 -check all -check bounds -g#-diag-file-append=intel_build.log -fpe0

# Makefile
.phony: clean

test_linfit: $(OBJs)
	module load intel/2018x; \
	$(f90comp) $(switch) -o test_linfit.exe $(OBJs) $(MKLLIB) 
%.o: %.f90
	module load intel/2018x; \
	$(f90comp) -c $(switch) $<
%.mod: %.o %.f90
	$(f90comp) -c $(switch) $<
clean:
	rm -f *.out
	rm -f *.log
	rm -f *.mod
	rm -f *.o
	rm -f *.exe

# End of Makefile
