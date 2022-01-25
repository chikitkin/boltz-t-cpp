#ifndef FULL_H_
#define FULL_H_

#include "header.h"

class Full {
public:
	// Constructors
	// Default constructor
	Full();
	// Zero tensor with given ranks
	Full(int n1_, int n2_, int n3_, int r1_, int r2_, int r3_);
	// Compress a tensor with given accuracy
	Full(int n1_, int n2_, int n3_, double *a, double eps = 1e-14);
	// Create a rank-1 tensor from given factors
	Full(int n1_, int n2_, int n3_, double *u1_, double *u2_, double *u3_);
	// Copy constructor
	Full(const Full& t);
	// Copy assignment operator
	Full& operator =(const Full& t);
	// Destructor
	~Full();

	// Print ranks
	std::vector<int> n() const;
	std::vector<int> r() const;
	// Print tensor
	friend std::ostream& operator << (std::ostream &out, const Full& t);
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

	friend Full operator -(const Full& t);
	// Element-wise summation
	friend Full operator +(const Full& t1, const Full& t2);
	friend Full operator -(const Full& t1, const Full& t2);
	// Element-wise multiplication
	friend Full operator *(const Full& t1, const Full& t2);
	friend Full operator *(const double alpha, const Full& t);
	friend Full operator *(const Full& t, const double alpha);
	friend Full operator /(const Full& t1, const Full& t2);
	friend Full reflect(const Full& t, char axis);

	friend Full round_t(const Full& t, double tol, int rmax);

private:
	int I(int i1, int i2, int i3);
	std::vector <int> multiI(int I);

	double* g;
	// sizes along each dimension;
	int n1, n2, n3; // TODO was const
	// u ranks;
	int r1, r2, r3;
};

#endif /* FULL_H_ */
