default: all

main.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o main.o -c main.cpp

model.o: Model.cpp Model.h
	mpicxx -std=c++11 -Wall -o model.o -c Model.cpp

burgers.o: Burgers.cpp Burgers.h Model.h
	mpicxx -std=c++11 -Wall -o burgers.o -c Burgers.cpp

compile: main.o  model.o burgers.o
	mpicxx -o HPC main.o  model.o burgers.o -O3 -ffast-math -funroll-loops -march=native -ftree-vectorize

.PHONY: clean # Specify that ’clean’ is not a real file
	target

# Entries:            ax,  ay, b,  c,  Nt,  Nx,  Ny  Px Py <user Inputs: h or v> 

diff: compile
	mpiexec -np 1 HPC 0.0 0.0 0.0 1.0 4000 2001 2001 1 1
 
advx: compile
	mpiexec -np 1 HPC 1.0 0.0 0.0 0.0 4000 2001 2001 1 1

advy: compile
	mpiexec -np 1 HPC 0.0 1.0 0.0 0.0 4000 2001 2001 1 1 

burgers: compile
	mpiexec -np 1 HPC 1.0 0.5 1.0 0.02 4000 2001 2001 1 1


# Entries:            ax,  ay, b,  c,  Nt,  Nx,  Ny  Px Py <user Inputs: h or v> 

diffp: compile
	mpiexec -np 2 HPC 0.0 0.0 0.0 1.0 4000 2001 2001 2 1
 
advxp: compile
	mpiexec -np 2 HPC 1.0 0.0 0.0 0.0 4000 2001 2001 2 1

advyp: compile
	mpiexec -np 2 HPC 0.0 1.0 0.0 0.0 4000 2001 2001 2 1 

burgersp: compile
	mpiexec -np 2 HPC 1.0 0.5 1.0 0.02 4000 2001 2001 2 1


# Entries:            ax,  ay, b,  c,  Nt,  Nx,  Ny  Px Py <user Inputs: h or v> 



diffn: compile
	mpiexec -np 36 HPC 0.0 0.0 0.0 1.0 4000 2001 2001 6 6
 
advxn: compile
	mpiexec -np 36 HPC 1.0 0.0 0.0 0.0 4000 2001 2001 6 6

advyn: compile
	mpiexec -np 36 HPC 0.0 1.0 0.0 0.0 4000 2001 2001 6 6 

burgersn: compile
	mpiexec -np 36 HPC 1.0 0.5 1.0 0.02 4000 2001 2001 6 6  


# Entries:            ax,  ay, b,  c,  Nt,  Nx,  Ny  Px Py <user Inputs: h or v> 

clean:
	rm -f *.o HPC   # Clean up (and ignore any errors)

all: diff advx advy burgers diffp advxp advyp burgersp diffn advxn advyn burgersn clean