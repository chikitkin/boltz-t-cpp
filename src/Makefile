CC=icc

CFLAGS=-g -Wall -O0 -mkl=sequential -DMKL_ILP64 -I"${MKLROOT}/include" -I"/home/egor/eigen-3.4.0" -I"/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/MEDCOUPLING/include" -I"/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/medfile/include" -I"/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/hdf5/include" -L"/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/MEDCOUPLING/lib" -L"/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/medfile/lib" -L"/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/hdf5/lib"

all: pre-build cli.o

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
	
mesh2.o: mesh2.cpp
	$(CC) $(CFLAGS) mesh2.cpp -o mesh2.o
	
solver.o: solver.cpp
	$(CC) $(CFLAGS) -c solver.cpp -o solver.o
	
tensor-test.o: tensor-test.cpp full.o tucker.o solver.o
	$(CC) $(CFLAGS) full.o tucker.o solver.o tensor-test.cpp -o tensor-test.o

mesh-test.o: mesh-test.cpp mesh.o
	$(CC) $(CFLAGS) mesh.o mesh-test.cpp -o mesh-test.o
	
clean:
	rm -rf *.o
