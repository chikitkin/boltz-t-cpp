CC=icc

CFLAGS=-g -Wall -O0 -mkl=sequential -DMKL_ILP64 -I"${MKLROOT}/include" -I"/home/egor/eigen-3.4.0"

all: cli.o

#pre-build:
#	$(source /opt/intel/oneapi/setvars.sh)
	
cli.o: full.o tucker.o read_starcd.o solver.o cli.cpp
	$(CC) $(CFLAGS) full.o tucker.o read_starcd.o solver.o cli.cpp -o cli.o
	
full.o: full.cpp
	$(CC) $(CFLAGS) -c full.cpp -o full.o

tucker.o: tucker.cpp
	$(CC) $(CFLAGS) -c tucker.cpp -o tucker.o
	
read_starcd.o: read_starcd.cpp
	$(CC) $(CFLAGS) -c read_starcd.cpp -o read_starcd.o
	
mesh2.o: mesh2.cpp
	$(CC) $(CFLAGS) mesh2.cpp -o mesh2.o
	
solver.o: solver.cpp
	$(CC) $(CFLAGS) -c solver.cpp -o solver.o
	
tensor-test.o: tensor-test.cpp full.o tucker.o solver.o
	$(CC) $(CFLAGS) full.o tucker.o solver.o tensor-test.cpp -o tensor-test.o

mesh-test.o: mesh-test.cpp read_starcd.o
	$(CC) $(CFLAGS) read_starcd.o mesh-test.cpp -o mesh-test.o
	
clean:
	rm -rf *.o
