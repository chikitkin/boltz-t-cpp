CC=icc

CFLAGS=-g -Wall -O0 -mkl=sequential -DMKL_ILP64 -I"${MKLROOT}/include"

all: run.o

run.o: full.o tucker.o mesh.o solver.o run.cpp
	$(CC) $(CFLAGS) full.o tucker.o mesh.o solver.o run.cpp -o run.o
	
full.o: full.cpp
	$(CC) $(CFLAGS) -c full.cpp -o full.o

tucker.o: tucker.cpp
	$(CC) $(CFLAGS) -c tucker.cpp -o tucker.o
	
mesh.o: mesh.cpp
	$(CC) $(CFLAGS) -c mesh.cpp -o mesh.o
	
solver.o: solver.cpp
	$(CC) $(CFLAGS) -c solver.cpp -o solver.o
	
tensor-test.o: tensor-test.cpp tensor.o
	$(CC) $(CFLAGS) tensor.o tensor-test.cpp -o tensor-test.o
	
mesh-test.o: mesh-test.cpp mesh.o
	$(CC) $(CFLAGS) mesh.o mesh-test.cpp -o mesh-test.o
	
clean:
	rm -rf *.o
