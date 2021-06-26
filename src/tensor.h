#ifndef TENSOR_H_
#define TENSOR_H_

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include <math.h>
#include <algorithm>
using namespace std;

double *svd_trunc(int m, int n, double *a, double eps, int &r);
double *svd_trunc_rmax(int m, int n, double *a, int rmax);

double **compress(int n1, int n2, int n3,
		double *a, double eps, int &r1, int &r2, int &r3, int rmax = 1e+6);

double **qr(int m, int n, const double *a);

class Tensor {
public:
	// Constructors
	// Default constructor
	Tensor();
	// Zero tensor with given ranks
	Tensor(int n1_, int n2_, int n3_, int r1_, int r2_, int r3_);
	// Compress a tensor with given accuracy
	Tensor(int n1_, int n2_, int n3_, double *a, double eps = 1e-14);
	// Create a rank-1 tensor from given factors
	Tensor(int n1_, int n2_, int n3_, double *u1_, double *u2_, double *u3_);
	// Copy constructor
	Tensor(const Tensor& t);
	// Copy assignment operator
	Tensor& operator =(const Tensor& t);
	// Destructor
	~Tensor();

	// Print ranks
	vector<int> n() const;
	vector<int> r() const;
	// Print tensor
	friend ostream& operator << (ostream &out, const Tensor& t);
	// Get element
	double At(int i1, int i2, int i3) const;

	// Orthogonalize factors with QR
	void orthogonalize();
	// Recompress tensor
	void round(double eps = 1e-14, int rmax = 1e+6);
	double *full() const;
	// Compute sum of all elements
	double sum() const;
	double norm();

	friend Tensor operator -(const Tensor& t);
	// Element-wise summation
	friend Tensor operator +(const Tensor& t1, const Tensor& t2);
	friend Tensor operator -(const Tensor& t1, const Tensor& t2);
	// Element-wise multiplication
	friend Tensor operator *(const Tensor& t1, const Tensor& t2);
	friend Tensor operator *(const double alpha, const Tensor& t);
	friend Tensor operator *(const Tensor& t, const double alpha);
	friend Tensor operator /(const Tensor& t1, const Tensor& t2);
	friend Tensor reflect(const Tensor& t, char axis);

	friend Tensor round_t(const Tensor& t, double tol = 1e-14, int rmax = 1e+6);

private:
	int I(int i1, int i2, int i3);
	vector <int> multiI(int I);

	double* g;
	// sizes along each dimension;
	int n1, n2, n3; // TODO was const
	// u ranks;
	int r1, r2, r3;
	// us;
	double* u1;
	double* u2;
	double* u3;
};

#endif /* TENSOR_H_ */
