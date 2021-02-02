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

double *svd_trunc(MKL_INT m, MKL_INT n, double *a, double eps, MKL_INT &r);

double **compress(MKL_INT n1, MKL_INT n2, MKL_INT n3,
		double *a, double eps, MKL_INT &r1, MKL_INT &r2, MKL_INT &r3);

double **qr(MKL_INT m, MKL_INT n, double *a);

class Tensor {
public:
	// Constructor
	Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_);
	// Zero tensor with given ranks
	Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, MKL_INT r1_, MKL_INT r2_, MKL_INT r3_);
	// Compress a tensor with given accuracy
	Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, double *a, double eps);
	// Create a rank-1 tensor from given factors
	Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, double *u1_, double *u2_, double *u3_);
	// Copy constructor
	Tensor(const Tensor& t);
	// Destructor
	~Tensor();

	// Print ranks
	void get_r();
	vector<int> shape() const;
	// Get element
	double At(MKL_INT i1, MKL_INT i2, MKL_INT i3);

	// Orthogonalize factors with QR
	void orthogonalize();
	// Recompress tensor
	void round(double eps);
	double *full();
	// Compute sum of all elements
	double sum();

	// Element-wise summation
	friend Tensor add(const Tensor& t1, const Tensor& t2);
	// Element-wise multiplication
	friend Tensor mult(const Tensor& t1, const Tensor& t2);
	Tensor& rmult(const double alpha);
	Tensor& divide(const Tensor& t);
	Tensor& operator =(const Tensor& t);

private:
	int I(int i1, int i2, int i3);
	vector<int> multiI(int I);

	double* g;
	// sizes along each dimension;
	const MKL_INT n1, n2, n3;
	// u ranks;
	MKL_INT r1, r2, r3;
	// us;
	double* u1;
	double* u2;
	double* u3;
};

Tensor operator +(const Tensor& t1, const Tensor& t2);
Tensor operator *(const Tensor& t1, const Tensor& t2);

#endif /* TENSOR_CLASS_H_ */
