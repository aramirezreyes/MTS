#FC=ifort  # Compilation
#LD=ifort  # Linker
FC=gfortran  # Compilation
LD=gfortran  # Linker
LIBS=  -L./libraries -lsparsekit
INCLUDE= -I./libraries 
EXE=pw1mts.exe

OBJS = DFLib.o pw1mts.o

all:   $(OBJS) 
	@echo "Links edition"
	$(LD) $(INCLUDE) $(LDFLAGS) $(OBJS) $(LIBS)  -o $(EXE) 

pw1mts.o : pw1mts.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c pw1mts.f90

DFLib.o : DFLib.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c DFLib.f90

clean: 
	rm -f *.o *.mod *.dat *~ $(EXE) 
 
%.o : %.mod 

