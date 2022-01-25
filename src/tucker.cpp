#include "header.h"
#include "tucker.h"

void transpose(int n1, int n2, int n3, double* a, int dim)
{
	switch (dim) // TODO enum
	{
		case 120: {
			MKL_Dimatcopy ('R', 'T', n1, n2 * n3, 1.0, a, n2 * n3, n1);
			break;
		}
		case 201: {
			MKL_Dimatcopy ('R', 'T', n1 * n2, n3, 1.0, a, n3, n1 * n2);
			break;
		}
		case 210: {
			double* a_tmp = new double [n1 * n2 * n3];
			LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1 * n2 * n3, 1, a, 1, a_tmp, 1);
			for (int i = 0; i < n1; ++i) {
				for (int j = 0; j < n2; ++j) {
					for (int k = 0; k < n3; ++k) {
						a[k * n2 * n1 + j * n1 + i] = a_tmp[i * n2 * n3 + j * n3 + k];
					}
				}
			}
			delete [] a_tmp;
			break;
		}
		default: {
			std::cout << "Wrong dim, try 120, 201, 210." << std::endl;
			exit(-1);
		}
	}
}
// Constructors
// Default constructor
Tucker::Tucker()
: n1(1), n2(1), n3(1)
{
//	std::cout << "Default constructor called" << std::endl;
	r1 = 1;
	r2 = 1;
	r3 = 1;
	g = new double[1]();
	u1 = new double[1];
	u2 = new double[1];
	u3 = new double[1];
}
// Zero tensor with given ranks
Tucker::Tucker(int n1_, int n2_, int n3_, int r1_, int r2_, int r3_)
: n1(n1_), n2(n2_), n3(n3_)
{
//	std::cout << "Zero constructor called" << std::endl;
	// Create zero tensor
	r1 = r1_;
	r2 = r2_;
	r3 = r3_;
	g = new double[r1 * r2 * r3](); // () IMPORTANT
	u1 = new double[n1 * r1];
	u2 = new double[n2 * r2];
	u3 = new double[n3 * r3];
}
// Compress a tensor with given accuracy
Tucker::Tucker(int n1_, int n2_, int n3_, double *a, double eps)
: n1(n1_), n2(n2_), n3(n3_)
{
//	std::cout << "Compress constructor called" << std::endl;
	double **A;
	A = compress(n1, n2, n3, a, eps, r1, r2, r3);
	g = new double[r1 * r2 * r3];
	u1 = new double[n1 * r1];
	u2 = new double[n2 * r2];
	u3 = new double[n3 * r3];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r1*r2*r3, 1, A[0], 1, g, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*r1, 1, A[1], 1, u1, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2*r2, 1, A[2], 1, u2, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3*r3, 1, A[3], 1, u3, 1);

	for (int i = 0; i < 4; ++i) {
		delete [] A[i];
	}
	delete [] A;
}
// Create a rank-1 tensor from given factors
Tucker::Tucker(int n1_, int n2_, int n3_, double *u1_, double *u2_, double *u3_)
: n1(n1_), n2(n2_), n3(n3_)
{
	r1 = 1;
	r2 = 1;
	r3 = 1;

	g = new double[r1 * r2 * r3];
	u1 = new double[n1 * r1];
	u2 = new double[n2 * r2];
	u3 = new double[n3 * r3];

	g[0] = 1.0;
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*r1, 1, u1_, 1, u1, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2*r2, 1, u2_, 1, u2, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3*r3, 1, u3_, 1, u3, 1);
}
// Copy constructor
Tucker::Tucker(const Tucker& t)
: n1(t.n1), n2(t.n2), n3(t.n3)
{
//	std::cout << "Copy constructor called" << std::endl;
	r1 = t.r1;
	r2 = t.r2;
	r3 = t.r3;

	g = new double[r1 * r2 * r3];
	u1 = new double[n1 * r1];
	u2 = new double[n2 * r2];
	u3 = new double[n3 * r3];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r1*r2*r3, 1, t.g, 1, g, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*r1, 1, t.u1, 1, u1, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2*r2, 1, t.u2, 1, u2, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3*r3, 1, t.u3, 1, u3, 1);
}
// Copy assignment operator
Tucker& Tucker::operator =(const Tucker& t)
{
//	std::cout << "Copy assignment operator called" << std::endl;
	if (this == &t)
		return *this;

	double* g_tmp = new double[t.r1 * t.r2 * t.r3];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t.r1 * t.r2 * t.r3, 1, t.g, 1, g_tmp, 1);
	delete [] g;
	g = g_tmp;

	double* u1_tmp = new double[t.n1 * t.r1];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t.n1 * t.r1, 1, t.u1, 1, u1_tmp, 1);
	delete [] u1;
	u1 = u1_tmp;

	double* u2_tmp = new double[t.n2 * t.r2];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t.n2 * t.r2, 1, t.u2, 1, u2_tmp, 1);
	delete [] u2;
	u2 = u2_tmp;

	double* u3_tmp = new double[t.n3 * t.r3];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t.n3 * t.r3, 1, t.u3, 1, u3_tmp, 1);
	delete [] u3;
	u3 = u3_tmp;

	n1 = t.n1;
	n2 = t.n2;
	n3 = t.n3;

	r1 = t.r1;
	r2 = t.r2;
	r3 = t.r3;

	return *this;
}
// Destructor
Tucker::~Tucker()
{
//	std::cout << "Destructor called" << std::endl;
	delete [] g;
	delete [] u1;
	delete [] u2;
	delete [] u3;
}

std::vector<int> Tucker::n() const
{
	std::vector<int> n = {n1, n2, n3};
	return n;
}

std::vector<int> Tucker::r() const
{
	std::vector<int> ranks = {r1, r2, r3};
	return ranks;
}
// Print tensor
std::ostream& operator << (std::ostream &out, const Tucker& t)
{
	out << "This is a 3D tensor in the Tucker format with \n";
	out << "r1 = " << t.r1 << ", n1 = " << t.n1 << "\n";
	out << "r2 = " << t.r2 << ", n2 = " << t.n2 << "\n";
	out << "r3 = " << t.r3 << ", n3 = " << t.n3;

	return out;
}
// Get element
double Tucker::At(int i1, int i2, int i3) const
{
	double a = 0.0;
	for(int j1 = 0; j1 < r1; j1++) {
		for(int j2 = 0; j2 < r2; j2++) {
			for(int j3 = 0; j3 < r3; j3++) {
				a += g[j1 * r2 * r3 + j2 * r3 + j3] * u1[i1 * r1 + j1] * u2[i2 * r2 + j2] * u3[i3 * r3 + j3];
			}
		}
	}
	return a;
}
// Orthogonalize factors with QR
void Tucker::orthogonalize()
{
	const double alpha = 1.0;
	const double beta = 0.0;

	double **QR1 = qr(n1, r1, u1);
	double **QR2 = qr(n2, r2, u2);
	double **QR3 = qr(n3, r3, u3);
	// TODO maybe just change pointers
	// u1 = q1
	delete [] u1;
	u1 = QR1[0];
	// u2 = q2
	delete [] u2;
	u2 = QR2[0];
	// u3 = q3
	delete [] u3;
	u3 = QR3[0];
	// z1.shape = (r1, r2, r3)

	double* g_tmp = new double[std::min(n1, r1)*std::min(n2, r2)*std::min(n3, r3)];

	// z1.shape = (r3, r2, r1)
	transpose(r1, r2, r3, g, 210);
	// R1.shape = (r1, std::min(n1, r1))
	MKL_Dimatcopy ('R', 'T', std::min(n1, r1), r1, alpha, QR1[1], r1, std::min(n1, r1));
	double *z2 = new double[r3*r2*std::min(n1, r1)];
	// z2.shape = (r3, r2, std::min(n1, r1))
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			r3*r2, std::min(n1, r1), r1, alpha, g, r1, QR1[1], std::min(n1, r1), beta, z2, std::min(n1, r1));

	delete [] g;
	delete [] QR1[1];
	delete [] QR1;

	// z2.shape = (std::min(n1, r1), r2, r3)
	transpose(r3, r2, std::min(n1, r1), z2, 210); // TODO
	// z2.shape = (r3, std::min(n1, r1), r2)
	transpose(std::min(n1, r1), r2, r3, z2, 201);
	// R2.shape = (r2, std::min(n2, r2))
	MKL_Dimatcopy ('R', 'T', std::min(n2, r2), r2, alpha, QR2[1], r2, std::min(n2, r2));
	double *z3 = new double[r3*std::min(n1, r1)*std::min(n2, r2)];
	// z3.shape = (r3, std::min(n1, r1), std::min(n2, r2))
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			r3*std::min(n1, r1), std::min(n2, r2), r2, alpha, z2, r2, QR2[1], std::min(n2, r2), beta, z3, std::min(n2, r2));

	delete [] z2;
	delete [] QR2[1];
	delete [] QR2;

	// z3.shape = (std::min(n1, r1), std::min(n2, r2), r3)
	transpose(r3, std::min(n1, r1), std::min(n2, r2), z3, 120);
	// R3.shape = (r3, std::min(n3, r3))
	MKL_Dimatcopy ('R', 'T', std::min(n3, r3), r3, alpha, QR3[1], r3, std::min(n3, r3));
	// g_tmp.shape = (std::min(n1, r1), std::min(n2, r2), std::min(n3, r3))
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			std::min(n1, r1)*std::min(n2, r2), std::min(n3, r3), r3, alpha, z3, r3, QR3[1], std::min(n3, r3), beta, g_tmp, std::min(n3, r3));

	delete [] z3;
	delete [] QR3[1];
	delete [] QR3;

	g = g_tmp;

	r1 = std::min(n1, r1);
	r2 = std::min(n2, r2);
	r3 = std::min(n3, r3);
}
// Recompress tensor
void Tucker::round(double eps, int rmax)
{
	const double alpha = 1.0;
	const double beta = 0.0;

	int r1_tmp, r2_tmp, r3_tmp;
	orthogonalize();
	// TODO: replace with tuple
	double **A;
	A = compress(r1, r2, r3, g, eps, r1_tmp, r2_tmp, r3_tmp, rmax);

	delete [] g;
	g = A[0];

	double *u1_tmp = new double[n1 * r1_tmp];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
									n1, r1_tmp, r1, alpha, u1, r1, A[1], r1_tmp, beta, u1_tmp, r1_tmp);
	delete [] u1;
	u1 = u1_tmp;

	delete [] A[1];

	double *u2_tmp = new double[n2 * r2_tmp];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
									n2, r2_tmp, r2, alpha, u2, r2, A[2], r2_tmp, beta, u2_tmp, r2_tmp);
	delete [] u2;
	u2 = u2_tmp;

	delete [] A[2];

	double *u3_tmp = new double[n3 * r3_tmp];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
									n3, r3_tmp, r3, alpha, u3, r3, A[3], r3_tmp, beta, u3_tmp, r3_tmp);
	delete [] u3;
	u3 = u3_tmp;

	delete [] A[3];

	r1 = r1_tmp;
	r2 = r2_tmp;
	r3 = r3_tmp;

	delete [] A;
}

double *Tucker::full() const
{
	const double alpha = 1.0;
	const double beta = 0.0;

	double *z1 = new double[n1*r2*r3];
	double *z2 = new double[n1*n2*r3];
	double *res = new double[n1*n2*n3];

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n1, r2*r3, r1, alpha, u1, r1, g, r2*r3, beta, z1, r2*r3);

	MKL_Dimatcopy ('R', 'T', n1, r2*r3, alpha, z1, r2*r3, n1);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n2, r3*n1, r2, alpha, u2, r2, z1, r3*n1, beta, z2, r3*n1);

	MKL_Dimatcopy ('R', 'T', n2, r3*n1, alpha, z2, r3*n1, n2);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n3, n1*n2, r3, alpha, u3, r3, z2, n1*n2, beta, res, n1*n2);

	MKL_Dimatcopy ('R', 'T', n3, n1*n2, alpha, res, n1*n2, n3);

	delete [] z1;
	delete [] z2;

	return res;
}
// Compute sum of all elements
double Tucker::sum() const
{
	const double alpha = 1.0;
	const double beta = 0.0;

	double *S;
	Tucker tmp(1, 1, 1, r1, r2, r3);
	MKL_Domatcopy ('R', 'N', r1*r2, r3, alpha, g, r3, tmp.g, r3);
	int n = std::max({n1, n2, n3});
	double* ones = new double[n];
	for (int i = 0; i < n; ++i) {
		ones[i] = 1.0;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
											1, r1, n1, alpha, ones, n, u1, r1, beta, tmp.u1, r1);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
											1, r2, n2, alpha, ones, n, u2, r2, beta, tmp.u2, r2);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
											1, r3, n3, alpha, ones, n, u3, r3, beta, tmp.u3, r3);

	delete [] ones;

	S = tmp.full();
	return S[0];
}

double Tucker::norm()
{
	orthogonalize();

	double norm = LAPACKE_dlange (LAPACK_ROW_MAJOR, 'F',
			r1 * r2 * r3, 1, g, 1);

	return norm;
}

Tucker operator -(const Tucker& t)
{
	return (-1.0) * t;
}

Tucker operator +(const Tucker& t1, const Tucker& t2) // TODO return reference (or not)
{
	// check that shapes are equal
	if (t1.n() != t2.n()){
		std::cout << "Different shapes in sum!" << std::endl;
		exit(-1);
	}
	Tucker result(t1.n1, t1.n2, t1.n3, t1.r1+t2.r1, t1.r2+t2.r2, t1.r3+t2.r3);

	for (int i = 0; i < t1.r1; ++i) {
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.r2, t1.r3, t1.g+i*t1.r2*t1.r3, t1.r3, result.g+
				i*(t1.r2+t2.r2)*(t1.r3+t2.r3),(t1.r3+t2.r3));
	}
	for (int i = 0; i < t2.r1; ++i) {
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.r2, t2.r3, t2.g+i*t2.r2*t2.r3, t2.r3, result.g+
				t1.r1*(t1.r2+t2.r2)*(t1.r3+t2.r3) + t1.r2*(t1.r3+t2.r3) + t1.r3 + i*(t1.r2+t2.r2)*(t1.r3+t2.r3), (t1.r3+t2.r3));
	}

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.n1, t1.r1, t1.u1, t1.r1, result.u1, t1.r1+t2.r1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.n1, t2.r1, t2.u1, t2.r1, result.u1+t1.r1, t1.r1+t2.r1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.n2, t1.r2, t1.u2, t1.r2, result.u2, t1.r2+t2.r2);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.n2, t2.r2, t2.u2, t2.r2, result.u2+t1.r2, t1.r2+t2.r2);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.n3, t1.r3, t1.u3, t1.r3, result.u3, t1.r3+t2.r3);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.n3, t2.r3, t2.u3, t2.r3, result.u3+t1.r3, t1.r3+t2.r3);

	return result;
}

Tucker operator -(const Tucker& t1, const Tucker& t2) // TODO return reference (or not)
{
	// check that shapes are equal
	if (t1.n() != t2.n()){
		std::cout << "Different shapes in sum!" << std::endl;
		exit(-1);
	}

	return t1 + (-1.0) * t2;
}
/*// TODO
Tucker& operator+=(Tucker& t, const Tucker& t1)
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
Tucker operator *(const Tucker& t1, const Tucker& t2)
{
	// check that shapes are equal;
	if (t1.n() != t2.n()){
		std::cout << "Different shapes in mult!" << std::endl;
		exit(-1);
	}
	Tucker result(t1.n1, t1.n2, t1.n3, t1.r1*t2.r1, t1.r2*t2.r2, t1.r3*t2.r3);

	for (int i = 0; i < t1.r1; ++i) {
		for (int j = 0; j < t1.r2; ++j) {
			for (int k = 0; k < t1.r3; ++k) {
				for (int l = 0; l < t2.r1; ++l) {
					MKL_Domatcopy ('R', 'N', t2.r2, t2.r3, t1.g[i*t1.r3*t1.r2 + j*t1.r3 + k],
							t2.g+l*t2.r2*t2.r3, t2.r3, result.g+result.r2*result.r3*i*t2.r1+result.r3*j*t2.r2+k*t2.r3+
							result.r2*result.r3*l, result.r3);
				}
			}
		}
	}

	for (int i = 0; i < t1.n1; ++i) {
		for (int j = 0; j < t1.r1; ++j) {
			MKL_Domatcopy ('R', 'N', 1, t2.r1, t1.u1[i*t1.r1 + j],
					t2.u1+i*t2.r1, t2.r1, result.u1+result.r1*i+t2.r1*j, result.r1);
		}
	}

	for (int i = 0; i < t1.n2; ++i) {
		for (int j = 0; j < t1.r2; ++j) {
			MKL_Domatcopy ('R', 'N', 1, t2.r2, t1.u2[i*t1.r2 + j],
					t2.u2+i*t2.r2, t2.r2, result.u2+result.r2*i+t2.r2*j, result.r2);
		}
	}

	for (int i = 0; i < t1.n3; ++i) {
		for (int j = 0; j < t1.r3; ++j) {
			MKL_Domatcopy ('R', 'N', 1, t2.r3, t1.u3[i*t1.r3 + j],
					t2.u3+i*t2.r3, t2.r3, result.u3+result.r3*i+t2.r3*j, result.r3);
		}
	}

	return result;
}

Tucker operator *(const double alpha, const Tucker& t)
{
	Tucker res(t);

	MKL_Dimatcopy('R', 'N', 1, t.r1 * t.r2 * t.r3, alpha, res.g, 1, t.r1 * t.r2 * t.r3);

	return res;
}

Tucker operator *(const Tucker& t, const double alpha)
{
	return alpha * t;
}

Tucker operator /(const Tucker& t1, const Tucker& t2)
{
	if (t2.r() != std::vector<int>{1, 1, 1}) {
		std::cout << "Invalid ranks in divide!" << std::endl;
		exit(-1);
	}

	Tucker t = t1;

	MKL_Dimatcopy('R', 'N', 1, t1.r1 * t1.r2 * t1.r3, 1.0 / t2.g[0], t.g, 1, t1.r1 * t1.r2 * t1.r3);
	for (int i = 0; i < t1.n1; ++i){
		MKL_Dimatcopy('R', 'N', 1, t1.r1, 1.0 / t2.u1[i], t.u1 + i * t1.r1, 1, t1.r1);
	}
	for (int i = 0; i < t1.n2; ++i){
		MKL_Dimatcopy('R', 'N', 1, t1.r2, 1.0 / t2.u2[i], t.u2 + i * t1.r2, 1, t1.r2);
	}
	for (int i = 0; i < t1.n3; ++i){
		MKL_Dimatcopy('R', 'N', 1, t1.r3, 1.0 / t2.u3[i], t.u3 + i * t1.r3, 1, t1.r3);
	}

	return t;
}

Tucker reflect(const Tucker& t, char axis)
{
	Tucker res(t);
	switch (axis) // TODO use enum
	{
	case 'X':
		for (int i = 0; i < t.n1; ++i) {
			LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', 1, t.r1,
					t.u1 + i * t.r1, t.r1, res.u1 + (t.n1 - 1 - i) * t.r1, t.r1);
		}
		return res;
	case 'Y':
		for (int i = 0; i < t.n2; ++i) {
			LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', 1, t.r2,
					t.u2 + i * t.r2, t.r2, res.u2 + (t.n2 - 1 - i) * t.r2, t.r2);
		}
		return res;
	case 'Z':
		for (int i = 0; i < t.n3; ++i) {
			LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', 1, t.r3,
					t.u3 + i * t.r3, t.r3, res.u3 + (t.n3 - 1 - i) * t.r3, t.r3);
		}
		return res;
	default:
		std::cout << "Wrong axis, try X, Y, Z." << std::endl;
		exit(-1);
	}
}

Tucker round_t(const Tucker& t, double tol = 1e-14, int rmax = 1000000)
{
	Tucker res(t);

	res.round(tol, rmax);

	return res;
}

int Tucker::I(int i1, int i2, int i3)
{
	return i1 * n2 * n3 + i2 * n3 + i3;
}

std::vector<int> Tucker::multiI(int I)
{
	// TODO: implement
	return {0, 0, 0};
}

double *svd_trunc(int m, int n, double *a, double eps, int &r)
{
	int info;

	double *U = new double[m*m];
	double *S = new double[std::min(m,n)];
	double *VT = new double[n*n];
	double *superb = new double[std::min(m,n)-1];

	double *a_copy = new double[m*n];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', m*n, 1, a, 1, a_copy, 1);
	info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'N', m, n, a_copy, n,
							S, U, m, VT, n, superb );

	if( info > 0 ) {
			printf( "The algorithm computing SVD failed to converge.\n" );
			exit( 1 );
	}

	eps = eps * S[0] / sqrt(3);

	r = std::min(m, n);

	for( int i = 0; i < std::min(m, n); ++i ) {
		if (S[i] <= eps) {
			r = i;
			break;
		}
	}

	double *u = new double[m*r]; // TODO can be optimized

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', m, r, U, m, u, r);

	delete[] a_copy;
	delete[] U;
	delete[] S;
	delete[] VT;
	delete[] superb;

	return u;
}

double *svd_trunc_rmax(int m, int n, double *a, int rmax)
{
	int info;

	double *U = new double[m*m];
	double *S = new double[std::min(m,n)];
	double *VT = new double[n*n];
	double *superb = new double[std::min(m,n)-1];

	double *a_copy = new double[m*n];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', m*n, 1, a, 1, a_copy, 1);
	info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'N', m, n, a, n,
							S, U, m, VT, n, superb );

	if( info > 0 ) {
			printf( "The algorithm computing SVD failed to converge.\n" );
			exit( 1 );
	}

	double *u = new double[m*rmax]; // TODO can be optimized

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', m, rmax, U, m, u, rmax);

	delete[] a_copy;
	delete[] U;
	delete[] S;
	delete[] VT;
	delete[] superb;

	return u;
}

double **compress(int n1, int n2, int n3, double *a, double eps, int &r1, int &r2, int &r3, int rmax)
{
	double **result = new double *[4]; //TODO: tuple

	const double alpha = 1.0;
	const double beta = 0.0;

	double *u1, *u2, *u3;

	double *z1 = new double[n1*n2*n3];
	double *z2 = new double[n1*n2*n3];
	double *z3 = new double[n1*n2*n3];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, z1, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, z2, 1);
	MKL_Dimatcopy ('R', 'T', n1, n2*n3, alpha, z2, n2*n3, n1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, z3, 1);
	MKL_Dimatcopy ('R', 'T', n1*n2, n3, alpha, z3, n3, n1*n2);

	u1 = svd_trunc(n1, (n2 * n3), z1, eps, r1);
	u2 = svd_trunc(n2, (n1 * n3), z2, eps, r2);
	u3 = svd_trunc(n3, (n1 * n2), z3, eps, r3);

	// TODO can just copy part
	if (r1 > rmax) {
		delete [] u1;
		u1 = svd_trunc_rmax(n1, (n2 * n3), z1, rmax);
		r1 = rmax;
	}
	if (r2 > rmax) {
		delete [] u2;
		u2 = svd_trunc_rmax(n2, (n1 * n3), z2, rmax);
		r2 = rmax;
	}
	if (r3 > rmax) {
		delete [] u3;
		u3 = svd_trunc_rmax(n3, (n1 * n2), z3, rmax);
		r3 = rmax;
	}

	double *g = new double[r1*r2*r3];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n1*n2, r3, n3, alpha, a, n3, u3, r3, beta, z3, r3); // TODO trans

	MKL_Dimatcopy ('R', 'T', n1*n2, r3, alpha, z3, r3, n1*n2);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							r3*n1, r2, n2, alpha, z3, n2, u2, r2, beta, z2, r2);

	MKL_Dimatcopy ('R', 'T', r3*n1, r2, alpha, z2, r2, r3*n1);
	// g = z2 x u1
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							r2*r3, r1, n1, alpha, z2, n1, u1, r1, beta, g, r1);

	MKL_Dimatcopy ('R', 'T', r2*r3, r1, alpha, g, r1, r2*r3);

	result[0] = g;
	result[1] = u1;
	result[2] = u2;
	result[3] = u3;

	delete[] z1;
	delete[] z2;
	delete[] z3;

	return result;
}
// QR-decomposition of matrix a
double **qr(int n1, int n2, const double *a)
{
	double **result = new double *[2];

	int size = std::min(n1, n2);

	double *tau = new double[size];

	double *a_tmp = new double[n1 * n2];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1, n2, a, n2, a_tmp, n2);
	// QR in-place
	LAPACKE_dgeqrf (LAPACK_ROW_MAJOR, n1, n2, a_tmp, n2, tau);

	double *q = new double[n1 * size];
	double *r = new double[size * n2]();
	// R (only upper triangle)
//	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'U', n1, n2, a_tmp, n2, r, n2);
	for (int i = 0; i < size; ++i) {
		for (int j = i; j < n2; ++j) {
			r[i*n2 + j] = a_tmp[i*n2 + j];
		}
	}
	// get Q
	LAPACKE_dorgqr (LAPACK_ROW_MAJOR, n1, size, size, a_tmp, n2, tau);
	// Q
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1, size, a_tmp, n2, q, size);

	result[0] = q;
	result[1] = r;

	delete [] a_tmp;
	delete [] tau;

	return result;
}
