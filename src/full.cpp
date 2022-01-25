#include "header.h"
#include "full.h"

// Constructors
// Default constructor
Full::Full()
: n1(1), n2(1), n3(1)
{
//	std::cout << "Default constructor called" << std::endl;
	r1 = 1;
	r2 = 1;
	r3 = 1;
	g = new double[1]();
}
// Zero tensor with given ranks
Full::Full(int n1_, int n2_, int n3_, int r1_, int r2_, int r3_)
: n1(n1_), n2(n2_), n3(n3_)
{
//	std::cout << "Zero constructor called" << std::endl;
	// Create zero tensor
	r1 = r1_;
	r2 = r2_;
	r3 = r3_;
	g = new double[n1 * n2 * n3](); // () IMPORTANT
}
// Compress a tensor with given accuracy
Full::Full(int n1_, int n2_, int n3_, double *a, double eps)
: n1(n1_), n2(n2_), n3(n3_)
{
//	std::cout << "Compress constructor called" << std::endl;
	r1 = 1;
	r2 = 1;
	r3 = 1;
	g = new double[n1 * n2 * n3];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, g, 1);
}
// Create a rank-1 tensor from given factors
Full::Full(int n1_, int n2_, int n3_, double *u1_, double *u2_, double *u3_)
: n1(n1_), n2(n2_), n3(n3_)
{
	r1 = 1;
	r2 = 1;
	r3 = 1;
	g = new double[n1 * n2 * n3];

	for (int i1 = 0; i1 < n1; ++i1) {
		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i3 = 0; i3 < n3; ++i3) {
				g[i1 * n2 * n3 + i2 * n3 + i3] = u1_[i1] * u2_[i2] * u3_[i3];
			}
		}
	}
}
// Copy constructor
Full::Full(const Full& t)
: n1(t.n1), n2(t.n2), n3(t.n3)
{
//	std::cout << "Copy constructor called" << std::endl;
	r1 = t.r1;
	r2 = t.r2;
	r3 = t.r3;

	g = new double[n1 * n2 * n3];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, t.g, 1, g, 1);
}
// Copy assignment operator
Full& Full::operator =(const Full& t)
{
//	std::cout << "Copy assignment operator called" << std::endl;
	if (this == &t)
		return *this;

	double* g_tmp = new double[t.n1 * t.n2 * t.n3];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t.n1 * t.n2 * t.n3, 1, t.g, 1, g_tmp, 1);
	delete [] g;
	g = g_tmp;

	n1 = t.n1;
	n2 = t.n2;
	n3 = t.n3;

	r1 = t.r1;
	r2 = t.r2;
	r3 = t.r3;

	return *this;
}
// Destructor
Full::~Full()
{
//	std::cout << "Destructor called" << std::endl;
	delete [] g;
}

std::vector<int> Full::n() const
{
	std::vector<int> n = {n1, n2, n3};
	return n;
}

std::vector<int> Full::r() const
{
	std::vector<int> ranks = {r1, r2, r3};
	return ranks;
}
// Print tensor
std::ostream& operator << (std::ostream &out, const Full& t)
{
	out << "This is a full 3D tensor with \n";
	out << "r1 = " << t.r1 << ", n1 = " << t.n1 << "\n";
	out << "r2 = " << t.r2 << ", n2 = " << t.n2 << "\n";
	out << "r3 = " << t.r3 << ", n3 = " << t.n3;

	return out;
}
// Get element
double Full::At(int i1, int i2, int i3) const
{
	return g[i1 * n2 * n3 + i2 * n3 + i3];
}
// Orthogonalize factors with QR
void Full::orthogonalize()
{
}
// Recompress tensor
void Full::round(double eps, int rmax)
{
	r1 = 1;
	r2 = 1;
	r3 = 1;
}

double *Full::full() const
{
	double *res = new double[n1*n2*n3];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1 * n2 * n3, 1, g, 1, res, 1);

	return res;
}
// Compute sum of all elements
double Full::sum() const
{
	double res = 0.0;
	for (int i1 = 0; i1 < n1; ++i1) {
		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i3 = 0; i3 < n3; ++i3) {
				res += g[i1 * n2 * n3 + i2 * n3 + i3];
			}
		}
	}
	return res;
}

double Full::norm()
{

	double norm = LAPACKE_dlange (LAPACK_ROW_MAJOR, 'F',
			n1 * n2 * n3, 1, g, 1);

	return norm;
}

Full operator -(const Full& t)
{
	return (-1.0) * t;
}

Full operator +(const Full& t1, const Full& t2) // TODO return reference (or not)
{
	// check that shapes are equal
	if (t1.n() != t2.n()){
		std::cout << "Different shapes in sum!" << std::endl;
		exit(-1);
	}
	Full result(t1.n1, t1.n2, t1.n3, t1.r1+t2.r1, t1.r2+t2.r2, t1.r3+t2.r3);

	int n1 = t1.n1;
	int n2 = t1.n2;
	int n3 = t1.n3;

	for (int i1 = 0; i1 < n1; ++i1) {
		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i3 = 0; i3 < n3; ++i3) {
				result.g[i1 * n2 * n3 + i2 * n3 + i3] = (t1.g[i1 * n2 * n3 + i2 * n3 + i3] + t2.g[i1 * n2 * n3 + i2 * n3 + i3]);
			}
		}
	}

	return result;
}

Full operator -(const Full& t1, const Full& t2) // TODO return reference (or not)
{
	// check that shapes are equal
	if (t1.n() != t2.n()){
		std::cout << "Different shapes in sum!" << std::endl;
		exit(-1);
	}

	return t1 + (-1.0) * t2;
}
/*// TODO
Full& operator+=(Full& t, const Full& t1)
{
	// TODO
	if (this == &t1)
		return 2.0 * (*this);

	if (t.n() != t1.n()){
		std::cout << "Different shapes!" << std::endl;
		exit(-1);
	}



}
*/
Full operator *(const Full& t1, const Full& t2)
{
	// check that shapes are equal
	if (t1.n() != t2.n()){
		std::cout << "Different shapes in sum!" << std::endl;
		exit(-1);
	}
	Full result(t1.n1, t1.n2, t1.n3, t1.r1+t2.r1, t1.r2+t2.r2, t1.r3+t2.r3);

	int n1 = t1.n1;
	int n2 = t1.n2;
	int n3 = t1.n3;

	for (int i1 = 0; i1 < n1; ++i1) {
		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i3 = 0; i3 < n3; ++i3) {
				result.g[i1 * n2 * n3 + i2 * n3 + i3] = (t1.g[i1 * n2 * n3 + i2 * n3 + i3] * t2.g[i1 * n2 * n3 + i2 * n3 + i3]);
			}
		}
	}

	return result;
}

Full operator *(const double alpha, const Full& t)
{
	Full res(t.n1, t.n2, t.n3, t.r1, t.r2, t.r3);

	int n1 = t.n1;
	int n2 = t.n2;
	int n3 = t.n3;

	for (int i1 = 0; i1 < n1; ++i1) {
		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i3 = 0; i3 < n3; ++i3) {
				res.g[i1 * n2 * n3 + i2 * n3 + i3] = (alpha * t.g[i1 * n2 * n3 + i2 * n3 + i3]);
			}
		}
	}

	return res;
}

Full operator *(const Full& t, const double alpha)
{
	return alpha * t;
}

Full operator /(const Full& t1, const Full& t2)
{
	Full result(t1.n1, t1.n2, t1.n3, 1, 1, 1);

	int n1 = t1.n1;
	int n2 = t1.n2;
	int n3 = t1.n3;

	for (int i1 = 0; i1 < n1; ++i1) {
		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i3 = 0; i3 < n3; ++i3) {
				result.g[i1 * n2 * n3 + i2 * n3 + i3] = (t1.g[i1 * n2 * n3 + i2 * n3 + i3] / t2.g[i1 * n2 * n3 + i2 * n3 + i3]);
			}
		}
	}

	return result;
}

Full reflect(const Full& t, char axis)
{
	Full res(t);

	int n1 = t.n1;
	int n2 = t.n2;
	int n3 = t.n3;

	switch (axis) // TODO use enum
	{
	case 'X':
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i2 = 0; i2 < n2; ++i2) {
				for (int i3 = 0; i3 < n3; ++i3) {
					res.g[i1 * n2 * n3 + i2 * n3 + i3] = t.g[(n1 - i1 - 1) * n2 * n3 + i2 * n3 + i3];
				}
			}
		}
		return res;
	case 'Y':
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i2 = 0; i2 < n2; ++i2) {
				for (int i3 = 0; i3 < n3; ++i3) {
					res.g[i1 * n2 * n3 + i2 * n3 + i3] = t.g[i1 * n2 * n3 + (n2 - i2 - 1) * n3 + i3];
				}
			}
		}
		return res;
	case 'Z':
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i2 = 0; i2 < n2; ++i2) {
				for (int i3 = 0; i3 < n3; ++i3) {
					res.g[i1 * n2 * n3 + i2 * n3 + i3] = t.g[i1 * n2 * n3 + i2 * n3 + (n3 - i3 - 1)];
				}
			}
		}
		return res;
	default:
		std::cout << "Wrong axis, try X, Y, Z." << std::endl;
		exit(-1);
	}
}

Full round_t(const Full& t, double tol = 1e-14, int rmax = 1000000)
{
	Full res(t);

	res.round(tol, rmax);

	return res;
}

int Full::I(int i1, int i2, int i3)
{
	return i1 * n2 * n3 + i2 * n3 + i3;
}

std::vector<int> Full::multiI(int I)
{
	// TODO: implement
	return {0, 0, 0};
}
