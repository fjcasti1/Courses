# HOMEWORK_3

#OBJs = createDebug.o libraryGKE.o knife_edge.o
OBJsBackward =  precision.o substitution.o testBackward.o
OBJsForward =  precision.o substitution.o testForward.o
f90comp = ifort
switch = -O3 #-fpe0 -check all -check bounds -g#-diag-file-append=intel_build.log -fpe0


# Makefile
.phony: clean

testForward: $(OBJsForward)
	module load intel/2018x; \
	$(f90comp) $(switch) -o testForward.exe $(OBJsForward)
testBackward: $(OBJsBackward)
	module load intel/2018x; \
	$(f90comp) $(switch) -o testBackward.exe $(OBJsBackward) 
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
