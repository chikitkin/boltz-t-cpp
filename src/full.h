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
	Full(int n1_, int n2_, int n3_, REAL *a, REAL eps = 1e-14);
	// Create a rank-1 tensor from given factors
	Full(int n1_, int n2_, int n3_, REAL *u1_, REAL *u2_, REAL *u3_);
	// Copy constructor
	Full(const Full& t);
	// Copy assignment operator
	Full& operator =(const Full& t);
	// Destructor
	~Full();

	// Print ranks
	std::vector<int> n() const;
	std::vector<int> r() const;
	REAL compression() const;
	// Print tensor
	friend std::ostream& operator << (std::ostream &out, const Full& t);
	// Get element
	REAL At(int i1, int i2, int i3) const;
	
	// Orthogonalize factors with QR
	void orthogonalize();
	// Recompress tensor
	void round(REAL tol = 1e-14, int rmax = 1000000);
	REAL *full() const;
	// Compute sum of all elements
	REAL sum() const;
	REAL norm();

	friend Full operator -(const Full& t);
	// Element-wise summation
	friend Full operator +(const Full& t1, const Full& t2);
	friend Full operator -(const Full& t1, const Full& t2);
	// Element-wise multiplication
	friend Full operator *(const Full& t1, const Full& t2);
	friend Full operator *(const REAL alpha, const Full& t);
	friend Full operator *(const Full& t, const REAL alpha);
	friend Full operator /(const Full& t1, const Full& t2);
	friend Full reflect(const Full& t, char axis);
	friend Full minmod(const Full& t1, const Full& t2);

	friend Full round_t(const Full& t, REAL tol, int rmax);
/*
	std::string to_string();
	friend Full from_string(std::string &tensor_string);
*/
private:
	int I(int i1, int i2, int i3);
	std::vector <int> multiI(int I);

	REAL* g;
	// sizes along each dimension;
	int n1, n2, n3; // TODO was const
	// u ranks;
	int r1, r2, r3;
};

#endif /* FULL_H_ */
