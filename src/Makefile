CC=icc

CFLAGS=-g -Wall -O0 -mkl=sequential -DMKL_ILP64 -I"${MKLROOT}/include"

all: run.o

run.o: tensor.o mesh.o solver.o run.cpp
	$(CC) $(CFLAGS) tensor.o mesh.o solver.o run.cpp -o run.o
	
tensor.o: tensor.cpp
	$(CC) $(CFLAGS) -c tensor.cpp -o tensor.o
	
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