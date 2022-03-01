CC=icc

CFLAGS=-g -Wall -O0 -qmkl=sequential -DMKL_ILP64 -I"${MKLROOT}/include" 

READ=read_mesh

all: cli.o

#pre-build:
#	$(source /opt/intel/oneapi/setvars.sh)
	
cli.o: full.o tucker.o $(READ).o solver.o cli.cpp
	$(CC) $(CFLAGS) full.o tucker.o $(READ).o solver.o cli.cpp -o cli.o
	
full.o: full.cpp
	$(CC) $(CFLAGS) -c full.cpp -o full.o

tucker.o: tucker.cpp
	$(CC) $(CFLAGS) -c tucker.cpp -o tucker.o
	
$(READ).o: $(READ).cpp
	$(CC) $(CFLAGS) -c $(READ).cpp -o $(READ).o
	
solver.o: solver.cpp
	$(CC) $(CFLAGS) -c solver.cpp -o solver.o
	
tensor-test.o: tensor-test.cpp full.o tucker.o $(READ).o solver.o
	$(CC) $(CFLAGS) full.o tucker.o $(READ).o solver.o tensor-test.cpp -o tensor-test.o

mesh-test.o: mesh-test.cpp $(READ).o
	$(CC) $(CFLAGS) $(READ).o mesh-test.cpp -o mesh-test.o
	
clean:
	rm -rf *.o
