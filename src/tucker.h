#ifndef TUCKER_H_
#define TUCKER_H_

#include "header.h"

REAL *svd_trunc(int m, int n, REAL *a, REAL eps, int &r);
REAL *svd_trunc_rmax(int m, int n, REAL *a, int rmax);

REAL **compress(int n1, int n2, int n3,
		REAL *a, REAL eps, int &r1, int &r2, int &r3, int rmax = 1e+6);

REAL **qr(int m, int n, const REAL *a);

class Tucker {
public:
	// Constructors
	// Default constructor
	Tucker();
	// Zero tensor with given ranks
	Tucker(int n1_, int n2_, int n3_, int r1_=1, int r2_=1, int r3_=1);
	// Compress a tensor with given accuracy
	Tucker(int n1_, int n2_, int n3_, REAL *a, REAL eps=1e-14);
	// Create a rank-1 tensor from given factors
	Tucker(int n1_, int n2_, int n3_, REAL *u1_, REAL *u2_, REAL *u3_);
	// Copy constructor
	Tucker(const Tucker& t);
	// Copy assignment operator
	Tucker& operator =(const Tucker& t);
	// Destructor
	~Tucker();

	// Print ranks
	std::vector<int> n() const;
	std::vector<int> r() const;
	REAL compression() const;
	// Print tensor
	friend std::ostream& operator << (std::ostream &out, const Tucker& t);
	// Get element
	REAL At(int i1, int i2, int i3) const;

	// Orthogonalize factors with QR
	void orthogonalize();
	// Recompress tensor
	void round(REAL tol=1e-14, int rmax=1000000);
	REAL *full() const;
	// Compute sum of all elements
	REAL sum() const;
	REAL norm();

	friend Tucker operator -(const Tucker& t);
	// Element-wise summation
	friend Tucker operator +(const Tucker& t1, const Tucker& t2);
	friend Tucker operator -(const Tucker& t1, const Tucker& t2);
	// Element-wise multiplication
	friend Tucker operator *(const Tucker& t1, const Tucker& t2);
	friend Tucker operator *(const REAL alpha, const Tucker& t);
	friend Tucker operator *(const Tucker& t, const REAL alpha);
	friend Tucker operator /(const Tucker& t1, const Tucker& t2);
	friend Tucker reflect(const Tucker& t, char axis);

	friend Tucker round_t(const Tucker& t, REAL tol, int rmax);
/*
	std::string to_string();
	friend Tucker from_string(std::string &tensor_string);
*/
private:
	int I(int i1, int i2, int i3);
	std::vector <int> multiI(int I);

	REAL* g;
	// sizes along each dimension;
	int n1, n2, n3; // TODO was const
	// u ranks;
	int r1, r2, r3;
	// us;
	REAL* u1;
	REAL* u2;
	REAL* u3;
};

#endif /* TUCKER_H_ */
