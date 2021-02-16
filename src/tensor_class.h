#ifndef TENSOR_CLASS_H_
#define TENSOR_CLASS_H_

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>
using namespace std;

double *svd_trunc(int m, int n, double *a, double eps, int &r);

double **compress(int n1, int n2, int n3,
		double *a, double eps, int &r1, int &r2, int &r3);

double **qr(int m, int n, double *a);

class Tensor {
public:
	// Constructor
	Tensor(int n1_, int n2_, int n3_);
	// Zero tensor with given ranks
	Tensor(int n1_, int n2_, int n3_, int r1_, int r2_, int r3_);
	// Compress a tensor with given accuracy
	Tensor(int n1_, int n2_, int n3_, double *a, double eps);
	// Create a rank-1 tensor from given factors
	Tensor(int n1_, int n2_, int n3_, double *u1_, double *u2_, double *u3_);
	// Copy constructor
	Tensor(const Tensor& t);
	// Destructor
	~Tensor();

	// Print ranks
	void get_r();
	vector<int> shape() const;
	// Get element
	double At(int i1, int i2, int i3);

	// Orthogonalize factors with QR
	void orthogonalize();
	// Recompress tensor
	void round(double eps);
	double *full();
	// Compute sum of all elements
	double sum() const;

	// Element-wise summation
	friend Tensor add(const Tensor& t1, const Tensor& t2);
	// Element-wise multiplication
	friend Tensor mult(const Tensor& t1, const Tensor& t2);
	Tensor& rmult(const double alpha);
	Tensor& divide(const Tensor& t);
	Tensor& operator =(const Tensor& t);
	friend Tensor reflect(const Tensor& t, char axis);

private:
	int I(int i1, int i2, int i3);
	vector<int> multiI(int I);

	double* g;
	// sizes along each dimension;
	const int n1, n2, n3;
	// u ranks;
	int r1, r2, r3;
	// us;
	double* u1;
	double* u2;
	double* u3;
};

Tensor operator +(const Tensor& t1, const Tensor& t2);
Tensor operator *(const Tensor& t1, const Tensor& t2);

#endif /* TENSOR_CLASS_H_ */
