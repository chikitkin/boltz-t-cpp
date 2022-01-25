#ifndef TUCKER_H_
#define TUCKER_H_

#include "header.h"

double *svd_trunc(int m, int n, double *a, double eps, int &r);
double *svd_trunc_rmax(int m, int n, double *a, int rmax);

double **compress(int n1, int n2, int n3,
		double *a, double eps, int &r1, int &r2, int &r3, int rmax = 1e+6);

double **qr(int m, int n, const double *a);

class Tucker {
public:
	// Constructors
	// Default constructor
	Tucker();
	// Zero tensor with given ranks
	Tucker(int n1_, int n2_, int n3_, int r1_, int r2_, int r3_);
	// Compress a tensor with given accuracy
	Tucker(int n1_, int n2_, int n3_, double *a, double eps = 1e-14);
	// Create a rank-1 tensor from given factors
	Tucker(int n1_, int n2_, int n3_, double *u1_, double *u2_, double *u3_);
	// Copy constructor
	Tucker(const Tucker& t);
	// Copy assignment operator
	Tucker& operator =(const Tucker& t);
	// Destructor
	~Tucker();

	// Print ranks
	std::vector<int> n() const;
	std::vector<int> r() const;
	// Print tensor
	friend std::ostream& operator << (std::ostream &out, const Tucker& t);
	// Get element
	double At(int i1, int i2, int i3) const;

	// Orthogonalize factors with QR
	void orthogonalize();
	// Recompress tensor
	void round(double eps = 1e-14, int rmax = 1000000);
	double *full() const;
	// Compute sum of all elements
	double sum() const;
	double norm();

	friend Tucker operator -(const Tucker& t);
	// Element-wise summation
	friend Tucker operator +(const Tucker& t1, const Tucker& t2);
	friend Tucker operator -(const Tucker& t1, const Tucker& t2);
	// Element-wise multiplication
	friend Tucker operator *(const Tucker& t1, const Tucker& t2);
	friend Tucker operator *(const double alpha, const Tucker& t);
	friend Tucker operator *(const Tucker& t, const double alpha);
	friend Tucker operator /(const Tucker& t1, const Tucker& t2);
	friend Tucker reflect(const Tucker& t, char axis);

	friend Tucker round_t(const Tucker& t, double tol, int rmax);

private:
	int I(int i1, int i2, int i3);
	std::vector <int> multiI(int I);

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

#endif /* TUCKER_H_ */
