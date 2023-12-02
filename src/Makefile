CC=icpx

CFLAGS=-std=c++14 -g -Wall -O3 -qopenmp -qmkl=sequential -DMKL_ILP64 -I"${MKLROOT}/include"

all: cli.o
	
#pre-build:
#	$(source /opt/intel/oneapi/setvars.sh)

cli.o: full.o tucker.o mesh.o solver.o cli.cpp
	$(CC) $(CFLAGS) full.o tucker.o mesh.o solver.o cli.cpp -o cli.o

full.o: full.cpp
	$(CC) $(CFLAGS) -c full.cpp -o full.o

tucker.o: tucker.cpp
	$(CC) $(CFLAGS) -c tucker.cpp -o tucker.o

mesh.o: mesh.cpp
	$(CC) $(CFLAGS) -c mesh.cpp -o mesh.o

solver.o: solver.cpp
	$(CC) $(CFLAGS) -c solver.cpp -o solver.o

tensor-test.o: tensor-test.cpp full.o tucker.o mesh.o solver.o
	$(CC) $(CFLAGS) full.o tucker.o mesh.o solver.o tensor-test.cpp -o tensor-test.o

mesh-test.o: mesh-test.cpp mesh.o
	$(CC) $(CFLAGS) mesh.o mesh-test.cpp -o mesh-test.o

clean:
	rm -rf *.o

