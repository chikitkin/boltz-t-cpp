#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>
using namespace std;
#include "tensor_class.h"

// Constructor
Tensor::Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_) :
	n1(n1_), n2(n2_), n3(n3_) {
	r1 = 0;
	r2 = 0;
	r3 = 0;
	g = new double[0];
	u1 = new double[0];
	u2 = new double[0];
	u3 = new double[0];
}
// Zero tensor with given ranks
Tensor::Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, MKL_INT r1_, MKL_INT r2_, MKL_INT r3_) :
	n1(n1_), n2(n2_), n3(n3_) {
	// Create zero tensor
	r1 = r1_;
	r2 = r2_;
	r3 = r3_;
	g = new double[r1 * r2 * r3]();
	u1 = new double[n1 * r1];
	u2 = new double[n2 * r2];
	u3 = new double[n3 * r3];
}
// Compress a tensor with given accuracy
Tensor::Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, double *a, double eps) :
	n1(n1_), n2(n2_), n3(n3_) {
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
Tensor::Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, double *u1_, double *u2_, double *u3_) :
	n1(n1_), n2(n2_), n3(n3_) {
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
// Destructor
Tensor::~Tensor(){
	delete [] g;
	delete [] u1;
	delete [] u2;
	delete [] u3;
}
// Print ranks
void Tensor::get_r() {
	cout << r1 << " " << r2 << " " << r3 << endl;
}
vector<int> Tensor::shape() const{
	auto dim = {n1, n2, n3};
	return dim;
}
// Get element
double Tensor::At(MKL_INT i1, MKL_INT i2, MKL_INT i3){
	double a = 0.0;
	for(MKL_INT j1 = 0; j1 < r1; j1++) {
		for(MKL_INT j2 = 0; j2 < r2; j2++) {
			for(MKL_INT j3 = 0; j3 < r3; j3++) {
				a += g[j1 * r2 * r3 + j2 * r3 + j3] * u1[i1 * r1 + j1] * u2[i2 * r2 + j2] * u3[i3 * r3 + j3];
			}
		}
	}
	return a;
}
// Orthogonalize factors with QR
void Tensor::orthogonalize() {

	const double alpha = 1.0;
	const double beta = 0.0;

	double **QR1, **QR2, **QR3;
	QR1 = qr(n1, r1, u1);
	QR2 = qr(n2, r2, u2);
	QR3 = qr(n3, r3, u3);

	delete [] u1;
	u1 = new double[n1*r1];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1, min(n1, r1), QR1[0], min(n1, r1), u1, min(n1, r1));
	delete [] u2;
	u2 = new double[n2*r2];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2, min(n2, r2), QR2[0], min(n2, r2), u2, min(n2, r2));
	delete [] u3;
	u3 = new double[n3*r3];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3, min(n3, r3), QR3[0], min(n3, r3), u3, min(n3, r3));

	double *z1 = new double[r1*r2*r3];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r1*r2*r3, 1, g, 1, z1, 1);

	delete [] g;
	g = new double[min(n1, r1)*min(n2, r2)*min(n3, r3)];

	mkl_dimatcopy ('R', 'T', r1, r2*r3, alpha, z1, r2*r3, r1);
	mkl_dimatcopy ('R', 'T', min(n1, r1), r1, alpha, QR1[1], r1, min(n1, r1));
	double *z2 = new double[min(n1, r1)*r2*r3];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
									r2*r3, min(n1, r1), r1, alpha, z1, r1, QR1[1], min(n1, r1), beta, z2, min(n1, r1));
	mkl_dimatcopy ('R', 'T', r2, min(n1, r1)*r3, alpha, z2, min(n1, r1)*r3, r2);
	mkl_dimatcopy ('R', 'T', min(n2, r2), r2, alpha, QR2[1], r2, min(n2, r2));
	double *z3 = new double[min(n1, r1)*min(n2, r2)*r3];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			min(n1, r1)*r3, min(n2, r2), r2, alpha, z2, r2, QR2[1], min(n2, r2), beta, z3, min(n2, r2));
	mkl_dimatcopy ('R', 'T', r3, min(n1, r1)*min(n2, r2), alpha, z3, min(n1, r1)*min(n2, r2), r3);
	mkl_dimatcopy ('R', 'T', min(n3, r3), r3, alpha, QR3[1], r3, min(n3, r3));
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			min(n1, r1)*min(n2, r2), min(n3, r3), r3, alpha, z3, r3, QR3[1], min(n3, r3), beta, g, min(n3, r3));

	delete [] z1;
	delete [] z2;
	delete [] z3;
	for (int i = 0; i < 2; ++i) {
		delete [] QR1[i];
		delete [] QR2[i];
		delete [] QR3[i];
	}
	delete [] QR1;
	delete [] QR2;
	delete [] QR3;

	r1 = min(n1, r1);
	r2 = min(n2, r2);
	r3 = min(n3, r3);

}
// Recompress tensor
void Tensor::round(double eps) {

	const double alpha = 1.0;
	const double beta = 0.0;

	MKL_INT r1_, r2_, r3_;

	orthogonalize();

	double **A;
	A = compress(r1, r2, r3, g, eps, r1_, r2_, r3_);

	delete [] g;
	g = new double[r1_ * r2_ * r3_];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r1_ * r2_ * r3_, 1, A[0], 1, g, 1);

	double *u1_ = new double[n1 * r1_];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
									n1, r1_, r1, alpha, u1, r1, A[1], r1_, beta, u1_, r1_);
	delete [] u1;
	u1 = new double[n1 * r1_];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1 * r1_, 1, u1_, 1, u1, 1);
	delete [] u1_;

	double *u2_ = new double[n2 * r2_];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
									n2, r2_, r2, alpha, u2, r2, A[2], r2_, beta, u2_, r2_);
	delete [] u2;
	u2 = new double[n2 * r2_];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2 * r2_, 1, u2_, 1, u2, 1);
	delete [] u2_;

	double *u3_ = new double[n3 * r3_];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
									n3, r3_, r3, alpha, u3, r3, A[3], r3_, beta, u3_, r3_);
	delete [] u3;
	u3 = new double[n3 * r3_];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3 * r3_, 1, u3_, 1, u3, 1);
	delete [] u3_;

	r1 = r1_;
	r2 = r2_;
	r3 = r3_;

	for (int i = 0; i < 4; ++i) {
		delete [] A[i];
	}
	delete [] A;

}

double *Tensor::full() {

	const double alpha = 1.0;
	const double beta = 0.0;

	double *z1 = new double[n1*r2*r3];
	double *z2 = new double[n1*n2*r3];
	double *res = new double[n1*n2*n3];

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n1, r2*r3, r1, alpha, u1, r1, g, r2*r3, beta, z1, r2*r3);

	mkl_dimatcopy ('R', 'T', n1, r2*r3, alpha, z1, r2*r3, n1);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n2, r3*n1, r2, alpha, u2, r2, z1, r3*n1, beta, z2, r3*n1);

	mkl_dimatcopy ('R', 'T', n2, r3*n1, alpha, z2, r3*n1, n2);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n3, n1*n2, r3, alpha, u3, r3, z2, n1*n2, beta, res, n1*n2);

	mkl_dimatcopy ('R', 'T', n3, n1*n2, alpha, res, n1*n2, n3);

	delete [] z1;
	delete [] z2;

	return res;
}
// Compute sum of all elements
double Tensor::sum() {

	const double alpha = 1.0;
	const double beta = 0.0;

	double *S;
	Tensor tmp(1, 1, 1, r1, r2, r3);
	mkl_domatcopy ('R', 'N', r1*r2, r3, alpha, g, r3, tmp.g, r3);
	MKL_INT n = max(n1, max(n2, n3));
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
	S = tmp.full();
	return S[0];
}

// Element-wise summation
Tensor add(const Tensor& t1, const Tensor& t2){
	// check that shapes are equal;
	if (t1.shape() != t2.shape()){
		cout << "Different shapes in sum!" << endl;
		exit(-1);
	}
	Tensor result(t1.n1, t1.n2, t1.n3, t1.r1+t2.r1, t1.r2+t2.r2, t1.r3+t2.r3);

	for (int i = 0; i < t1.r1; ++i) {
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.r2, t1.r3, t1.g+i*t1.r2*t1.r3, t1.r3, result.g+
				i*(t1.r2+t2.r2)*(t1.r3+t2.r3),(t1.r3+t2.r3));
	}
	for (int i = 0; i < t2.r1; ++i) {
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.r2, t2.r3, t2.g+i*t2.r2*t2.r3, t2.r3, result.g+
				t1.r1*(t1.r2+t2.r2)*(t1.r3+t2.r3) + t1.r2*(t1.r3+t2.r3) + t1.r3 + i*(t1.r2+t2.r2)*(t1.r3+t2.r3), (t1.r3+t2.r3));
	}

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.n1, t1.r1, t1.u1, t1.r1, result.u1, (t1.r1+t2.r1));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.n1, t2.r1, t2.u1, t2.r1, result.u1+t1.r1, (t1.r1+t2.r1));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.n2, t1.r2, t1.u2, t1.r2, result.u2, (t1.r2+t2.r2));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.n2, t2.r2, t2.u2, t2.r2, result.u2+t1.r2, (t1.r2+t2.r2));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.n3, t1.r3, t1.u3, t1.r3, result.u3, (t1.r3+t2.r3));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.n3, t2.r3, t2.u3, t2.r3, result.u3+t1.r3, (t1.r3+t2.r3));

	return result;
}

// Element-wise multiplication
Tensor mult(const Tensor& t1, const Tensor& t2){
	// check that shapes are equal;
	if (t1.shape() != t2.shape()){
		cout << "Different shapes in mult!" << endl;
		exit(-1);
	}
	Tensor result(t1.n1, t1.n2, t1.n3, t1.r1*t2.r1, t1.r2*t2.r2, t1.r3*t2.r3);

	for (int i = 0; i < t1.r1; ++i) {
		for (int j = 0; j < t1.r2; ++j) {
			for (int k = 0; k < t1.r3; ++k) {
				for (int l = 0; l < t2.r1; ++l) {
					mkl_domatcopy ('R', 'N', t2.r2, t2.r3, t1.g[i*t1.r3*t1.r2 + j*t1.r3 + k],
							t2.g+l*t2.r2*t2.r3, t2.r3, result.g+result.r2*result.r3*i*t2.r1+result.r3*j*t2.r2+k*t2.r3+
							result.r2*result.r3*l, result.r3);
				}
			}
		}
	}

	for (int i = 0; i < t1.n1; ++i) {
		for (int j = 0; j < t1.r1; ++j) {
			mkl_domatcopy ('R', 'N', 1, t2.r1, t1.u1[i*t1.r1 + j],
					t2.u1+i*t2.r1, t2.r1, result.u1+result.r1*i+t2.r1*j, result.r1);
		}
	}

	for (int i = 0; i < t1.n2; ++i) {
		for (int j = 0; j < t1.r2; ++j) {
			mkl_domatcopy ('R', 'N', 1, t2.r2, t1.u2[i*t1.r2 + j],
					t2.u2+i*t2.r2, t2.r2, result.u2+result.r2*i+t2.r2*j, result.r2);
		}
	}

	for (int i = 0; i < t1.n3; ++i) {
		for (int j = 0; j < t1.r3; ++j) {
			mkl_domatcopy ('R', 'N', 1, t2.r3, t1.u3[i*t1.r3 + j],
					t2.u3+i*t2.r3, t2.r3, result.u3+result.r3*i+t2.r3*j, result.r3);
		}
	}

	return result;
}

Tensor& Tensor::rmult(const double alpha){
	mkl_dimatcopy('R', 'N', 1, r1 * r2 * r3, alpha, g, 1, r1 * r2 * r3);

	return *this;
}

Tensor& Tensor::divide(const Tensor& t){
	if (vector<int>{t.r1, t.r2, t.r3} != vector<int>{1, 1, 1}){
		cout << "Invalid ranks in divide!" << endl;
		exit(-1);
	}
	mkl_dimatcopy('R', 'N', 1, r1 * r2 * r3, 1.0 / t.g[0], g, 1, r1 * r2 * r3);
	for (int i = 0; i < n1; ++i){
		mkl_dimatcopy('R', 'N', 1, r1, 1.0 / t.u1[i], u1+i*r1, 1, r1);
	}
	for (int i = 0; i < n2; ++i){
		mkl_dimatcopy('R', 'N', 1, r2, 1.0 / t.u2[i], u2+i*r2, 1, r2);
	}
	for (int i = 0; i < n3; ++i){
		mkl_dimatcopy('R', 'N', 1, r3, 1.0 / t.u3[i], u3+i*r3, 1, r3);
	}

	return *this;
}

Tensor& Tensor::operator =(const Tensor& t){

	if (this == &t)
			return *this;

	if (shape() != t.shape()){
		cout << "Different shapes!" << endl;
		exit(-1);
	}
	delete [] g;
	g = new double[t.r1 * t.r2 * t.r3];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t.r1 * t.r2 * t.r3, 1, t.g, 1, g, 1);

	delete [] u1;
	u1 = new double[t.n1 * t.r1];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t.n1 * t.r1, 1, t.u1, 1, u1, 1);

	delete [] u2;
	u2 = new double[t.n2 * t.r2];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t.n2 * t.r2, 1, t.u2, 1, u2, 1);

	delete [] u3;
	u3 = new double[t.n3 * t.r3];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t.n3 * t.r3, 1, t.u3, 1, u3, 1);

	r1 = t.r1;
	r2 = t.r2;
	r3 = t.r3;

	return *this;
}

// TODO: overload operators "+, *";
Tensor operator +(const Tensor& t1, const Tensor& t2){
	return add(t1, t2);
}

Tensor operator *(const Tensor& t1, const Tensor& t2){
	return mult(t1, t2);
}

int Tensor::I(int i1, int i2, int i3){
	return i1 * n2 * n3 + i2 * n3 + i3;
}
vector<int> Tensor::multiI(int I){
	// TODO: implement
}

double *svd_trunc(MKL_INT m, MKL_INT n, double *a, double eps, MKL_INT &r) {

	MKL_INT info;

	double *U = new double[m*m];
    double *S = new double[min(m,n)];
    double *VT = new double[n*n];
    double *superb = new double[min(m,n)-1];

	info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'N', m, n, a, n,
	                        S, U, m, VT, n, superb );

    if( info > 0 ) {
            printf( "The algorithm computing SVD failed to converge.\n" );
            exit( 1 );
    }

    eps = eps * S[0] / sqrt(3);

    r = min(m, n);

    for( int i = 0; i < min(m, n); ++i ) {
        if (S[i] <= eps) {
            r = i;
            break;
        }
    }

    double *u = new double[m*r];

    LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', m, r, U, m, u, r);

    delete[] U;
    delete[] S;
    delete[] VT;
    delete[] superb;

    return u;
}

double **compress(MKL_INT n1, MKL_INT n2, MKL_INT n3, double *a, double eps, MKL_INT &r1, MKL_INT &r2, MKL_INT &r3) {

	double **result = new double *[4]; //TODO: tuple

    const double alpha = 1.0;
    const double beta = 0.0;

    double *u1, *u2, *u3;

    double *z1 = new double[n1*n2*n3];
    double *z2 = new double[n1*n2*n3];
    double *z3 = new double[n1*n2*n3];

    LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, z1, 1);
    LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, z2, 1);
    mkl_dimatcopy ('R', 'T', n1, n2*n3, alpha, z2, n2*n3, n1);
    LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, z3, 1);
    mkl_dimatcopy ('R', 'T', n1*n2, n3, alpha, z3, n3, n1*n2);

    u1 = svd_trunc(n1, (n2 * n3), z1, eps, r1);
    u2 = svd_trunc(n2, (n1 * n3), z2, eps, r2);
    u3 = svd_trunc(n3, (n1 * n2), z3, eps, r3);

    double *g = new double[r1*r2*r3]();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            n1*n2, r3, n3, alpha, a, n3, u3, r3, beta, z3, r3);
    mkl_dimatcopy ('R', 'T', n1*n2, r3, alpha, z3, r3, n1*n2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            r3*n1, r2, n2, alpha, z3, n2, u2, r2, beta, z2, r2);
    mkl_dimatcopy ('R', 'T', r3*n1, r2, alpha, z2, r2, r3*n1);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            r2*r3, r1, n1, alpha, z2, n1, u1, r1, beta, g, r1);
    mkl_dimatcopy ('R', 'T', r2*r3, r1, alpha, g, r1, r2*r3);

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
double **qr(MKL_INT m, MKL_INT n, double *a) {

	double **result = new double *[2];

    const double alpha = 1.0;
    const double beta = 0.0;

	MKL_INT size = min(m, n);

	double *tau = new double[size];

	double *a_tmp = new double[m*n];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', m, n, a, n, a_tmp, n);

	LAPACKE_dgeqrf (LAPACK_ROW_MAJOR, m, n, a_tmp, n, tau);

	double *q = new double[m*size];
	double *r = new double[size*n]();

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'U', size, n, a_tmp, n, r, n);

	LAPACKE_dorgqr (LAPACK_ROW_MAJOR, m, size, size, a_tmp, n, tau);

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', m, size, a_tmp, n, q, size);

	result[0] = q;
	result[1] = r;

	delete [] a_tmp;
	delete [] tau;

	return result;
}

void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
		printf( "\n" );
	}
}

double f(double x, double y, double z) {
	return exp(-(x*x + y*y + z*z) / 1000.0) + (x / 100.0) + (x*y / 100.0) + 10.0 * sin(z) + exp(-(x*(x - 10.) + y*y + z*z) / 1000.0);
}

double g(double x, double y, double z) {
	return exp(-(x*x + y*y) / 100.0) + (x * z / 100.0) + (x*y / 100.0) + 10.0 * cos(x * y) + exp(-(x*(x - 5.) + y*y + z) / 1000.0);
}

int main(){
	MKL_INT n1, n2, n3;
	n1 = 10;
	n2 = 10;
	n3 = 20;

	double eps = 0.0005;
/*************************************************************************************/
	double *a1 = new double[n1*n2*n3];
	double *a2 = new double[n1*n2*n3];
	double *a3 = new double[n1*n2*n3];

	double s;
	double s_res;

	double a_norm;
	double diff_norm;
	double diff;
/*************************************************************************************/
	s = 0.0;
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a2[i*n2*n3 + j*n3 + k] = g(double(i), double(j), double(k));
				a3[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k)) + g(double(i), double(j), double(k));
				s += a3[i*n2*n3 + j*n3 + k];
			}
		}
	}
	Tensor t1_1(n1, n2, n3, a1, eps);
	Tensor t2_1(n1, n2, n3, a2, eps);
	Tensor t3_1(n1, n2, n3);
	t3_1 = t1_1 + t2_1; //TODO
	t3_1.round(eps);

	s_res = t3_1.sum();

	cout << "s_diff sum = " << (s_res - s) / s << endl;

	double *t3_1_full;
	t3_1_full = t3_1.full();

//	print_matrix("t3_full", n1*n2, n3, t3_full, n3);
//	print_matrix("a3", n1*n2, n3, a3, n3);
	a_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	cblas_daxpy(n1*n2*n3, -1.0, t3_1_full, 1, a3, 1);
	diff_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	diff = diff_norm / a_norm;
	cout << "diff sum = " << diff << endl;

	delete [] t3_1_full;
/*************************************************************************************/
	s = 0.0;
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a2[i*n2*n3 + j*n3 + k] = g(double(i), double(j), double(k));
				a3[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k)) * g(double(i), double(j), double(k));
				s += a3[i*n2*n3 + j*n3 + k];
			}
		}
	}
	Tensor t1_2(n1, n2, n3, a1, eps);
	Tensor t2_2(n1, n2, n3, a2, eps);
	Tensor t3_2(n1, n2, n3);
	t3_2 = t1_2 * t2_2;
	t3_2.round(eps);

	s_res = t3_2.sum();

	cout << "s_diff mult = " << (s_res - s) / s << endl;

	double *t3_2_full;
	t3_2_full = t3_2.full();

//	print_matrix("t3_full", n1*n2, n3, t3_full, n3);
//	print_matrix("a3", n1*n2, n3, a3, n3);
	a_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	cblas_daxpy(n1*n2*n3, -1.0, t3_2_full, 1, a3, 1);
	diff_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	diff = diff_norm / a_norm;
	cout << "diff mult = " << diff << endl;

	delete [] t3_2_full;
/*************************************************************************************/
	s = 0.0;
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a2[i*n2*n3 + j*n3 + k] = double(i + 1.0);
				a3[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k)) / double(i + 1.0);
				s += a3[i*n2*n3 + j*n3 + k];
			}
		}
	}
	Tensor t1_3(n1, n2, n3, a1, eps);
	Tensor t2_3(n1, n2, n3, a2, eps);
	Tensor t3_3(n1, n2, n3);
	t3_3 = t1_3.divide(t2_3);
	t3_3.round(eps);

	s_res = t3_3.sum();

	cout << "s_diff div = " << (s_res - s) / s << endl;

	double *t3_3_full;
	t3_3_full = t3_3.full();

//	print_matrix("t3_full", n1*n2, n3, t3_full, n3);
//	print_matrix("a3", n1*n2, n3, a3, n3);
	a_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	cblas_daxpy(n1*n2*n3, -1.0, t3_3_full, 1, a3, 1);
	diff_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	diff = diff_norm / a_norm;
	cout << "diff div = " << diff << endl;

	delete [] t3_3_full;
/*************************************************************************************/
	s = 0.0;
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a2[i*n2*n3 + j*n3 + k] = double(i + 1.0);
				a3[i*n2*n3 + j*n3 + k] = 2 * f(double(i), double(j), double(k));
				s += a3[i*n2*n3 + j*n3 + k];
			}
		}
	}
	Tensor t1_4(n1, n2, n3, a1, eps);
	Tensor t2_4(n1, n2, n3, a2, eps);
	Tensor t3_4(n1, n2, n3);
	t3_4 = t1_4.rmult(2.0);
	t3_4.round(eps);

	s_res = t3_4.sum();

	cout << "s_diff rmult = " << (s_res - s) / s << endl;

	double *t3_4_full;
	t3_4_full = t3_4.full();

//	print_matrix("t3_full", n1*n2, n3, t3_full, n3);
//	print_matrix("a3", n1*n2, n3, a3, n3);
	a_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	cblas_daxpy(n1*n2*n3, -1.0, t3_4_full, 1, a3, 1);
	diff_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	diff = diff_norm / a_norm;
	cout << "diff rmult = " << diff << endl;

	delete [] t3_4_full;
/*************************************************************************************/
	s = 0.0;
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a2[i*n2*n3 + j*n3 + k] = g(double(i), double(j), double(k));
				a3[i*n2*n3 + j*n3 + k] = 2 * f(double(i), double(j), double(k)) + g(double(i), double(j), double(k));
				s += a3[i*n2*n3 + j*n3 + k];
			}
		}
	}
	Tensor t1_5(n1, n2, n3, a1, eps);
	Tensor t2_5(n1, n2, n3, a2, eps);
	Tensor t3_5(n1, n2, n3);
	t3_5 = (t1_5 + t2_5) + t1_5; // TODO: temp?
	t3_5.round(eps);

	s_res = t3_5.sum();

	cout << "s_diff ++ = " << (s_res - s) / s << endl;

	double *t3_5_full;
	t3_5_full = t3_5.full();

//	print_matrix("t3_full", n1*n2, n3, t3_full, n3);
//	print_matrix("a3", n1*n2, n3, a3, n3);
	a_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	cblas_daxpy(n1*n2*n3, -1.0, t3_5_full, 1, a3, 1);
	diff_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	diff = diff_norm / a_norm;
	cout << "diff ++ = " << diff << endl;

	delete [] t3_5_full;
/*************************************************************************************/
	delete [] a1;
	delete [] a2;
	delete [] a3;

	return 0;
}
