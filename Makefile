default: all

main.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o main.o -c main.cpp

model.o: Model.cpp Model.h
	mpicxx -std=c++11 -Wall -o model.o -c Model.cpp

burgers.o: Burgers.cpp Burgers.h Model.h
	mpicxx -std=c++11 -Wall -o burgers.o -c Burgers.cpp

compile: main.o  model.o burgers.o
	mpicxx -o my_prog main.o  model.o burgers.o

.PHONY: clean # Specify that ’clean’ is not a real file
	target

diff: compile
	mpiexec -np 2 my_prog 1 0.5 1 0.02 2 1

clean:
	rm -f *.o my_prog   # Clean up (and ignore any errors)

all: diff clean