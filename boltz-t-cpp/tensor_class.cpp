#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>
using namespace std;

extern double *svd_trunc(MKL_INT m, MKL_INT n, double *a, double eps, MKL_INT &r);
extern double **compress(MKL_INT n1, MKL_INT n2, MKL_INT n3,
		double *a, double eps, MKL_INT &r1, MKL_INT &r2, MKL_INT &r3);
extern double **qr(MKL_INT m, MKL_INT n, double *a);

class Tensor {
public:
	// Constructor
	Tensor(){

	}
	// Zero tensor with given ranks
	Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, MKL_INT r1_, MKL_INT r2_, MKL_INT r3_){
		// Create zero tensor
		n1 = n1_;
		n2 = n2_;
		n3 = n3_;
		r1 = r1_;
		r2 = r2_;
		r3 = r3_;
		g = new double[r1 * r2 * r3]();
		u1 = new double[n1 * r1];
		u2 = new double[n2 * r2];
		u3 = new double[n3 * r3];
	}
	// Compress a tensor with given accuracy
	Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, double *a, double eps){
		n1 = n1_;
		n2 = n2_;
		n3 = n3_;
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
	// Print ranks
	void get_r() {
		cout << r1 << " " << r2 << " " << r3 << endl;
	}
	// Destructor
	~Tensor(){
		cout << "destructed" << endl;
		delete [] g;
		delete [] u1;
		delete [] u2;
		delete [] u3;
	}

	vector<int> shape() const{
		auto dim = {n1, n2, n3};
		return dim;
	}
	// Get element
	double At(MKL_INT i1, MKL_INT i2, MKL_INT i3){
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
	void orthogonalize() {

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
	void round(double eps) {

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

	double *full() {

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
	double sum() {

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
	//TODO: const, Tensor&, new
	friend Tensor add(const Tensor& t1, const Tensor& t2){
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
	friend Tensor mult(const Tensor& t1, const Tensor& t2){
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

private:
	int I(int i1, int i2, int i3){
		return i1 * n2 * n3 + i2 * n3 + i3;
	}
	vector<int> multiI(int I){
		// TODO: implement
	}
	// dynamic array to store parameters
	double* g;
	// sizes along each dimension;
	MKL_INT n1, n2, n3;
	// u ranks;
	MKL_INT r1, r2, r3;
	// us;
	double* u1;
	double* u2;
	double* u3;

};
// TODO: overload operators "+, *";
Tensor operator +(Tensor& t1, Tensor& t2){
	return add(t1, t2);
}

Tensor operator *(Tensor& t1, Tensor& t2){
	return mult(t1, t2);
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

	double **result = new double *[4];

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

int main(){
	MKL_INT n1, n2, n3;
	n1 = 10;
	n2 = 10;
	n3 = 20;

	double eps = 0.0001;

	double *a1 = new double[n1*n2*n3];
	double *a2 = new double[n1*n2*n3];
	double *a3 = new double[n1*n2*n3];
	double s = 0.0;

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a2[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a3[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k)) * f(double(i), double(j), double(k));
				s += a3[i*n2*n3 + j*n3 + k];
			}
		}
	}
	/*
	Tensor *t1p = new Tensor(n1, n2, n3, a1, eps);
	Tensor t1 = *t1p;
	Tensor *t2p = new Tensor(n1, n2, n3, a1, eps);
	Tensor t2 = *t2p;
	Tensor *t3p = new Tensor;
	Tensor t3 = *t3p;
	delete [] t3p;
	t3 = t1 * t2;
	t3.round(eps);

	delete [] t1p;
	delete [] t2p;
	delete [] t3p;
	*/
	Tensor t1 = Tensor(n1, n2, n3, a1, eps);
	Tensor t2 = Tensor(n1, n2, n3, a1, eps);
	Tensor* t3p = new Tensor;
	Tensor t3 = *t3p;
	delete t3p;
	Tensor t3 = t1 * t2;
//	t3.round(eps);

	double s_res = t3.sum();

	cout << "s_diff = " << s_res - s << endl;

	double *t3_full;
	t3_full = t3.full();

//	print_matrix("t3_full", n1*n2, n3, t3_full, n3);
//	print_matrix("a3", n1*n2, n3, a3, n3);
	double a_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	cblas_daxpy(n1*n2*n3, -1.0, t3_full, 1, a3, 1);
	double diff_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	double diff = diff_norm / a_norm;
	cout << "diff = " << diff << endl;

	delete [] a1;
	delete [] a2;
	delete [] a3;
	delete [] t3_full;

	return 0;
}
