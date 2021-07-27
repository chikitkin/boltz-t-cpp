#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>

using namespace std;

#define N1 5
#define N2 5
#define N3 5

double f(double x, double y, double z) {
	return exp(-(x*x + y*y + z*z) / 1000.0);// + (x / 100.0) + (x*y / 100.0) + 10.0 * sin(z) + exp(-(x*(x - 10.) + y*y + z*z) / 1000.0);
}

extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern double *svd_trunc(MKL_INT m, MKL_INT n, double *a, double eps, MKL_INT &r);
extern double **compress(MKL_INT n1, MKL_INT n2, MKL_INT n3, double *a, double eps, MKL_INT &r1, MKL_INT &r2, MKL_INT &r3);
extern void full(MKL_INT n1, MKL_INT n2, MKL_INT n3, MKL_INT r1, MKL_INT r2, MKL_INT r3, double **tuc, double *res);
extern void sum(MKL_INT n1, MKL_INT n2, MKL_INT n3,
		MKL_INT r11, MKL_INT r12, MKL_INT r13, double **tuc1,
		MKL_INT r21, MKL_INT r22, MKL_INT r23, double **tuc2,
		double **tuc3);

int main() {

	MKL_INT n1 = N1, n2 = N2, n3 = N3;

	double eps = -0.01;

	double *a = new double[n1*n2*n3];
	double *a_ = new double[n1*n2*n3];
	double *b = new double[n1*n2*n3];
	double *c = new double[n1*n2*n3];
	double *c_ = new double[n1*n2*n3];

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a[i*n2*n3 + j*n3 + k] = double(k);//double(i);//f(double(i), double(j), double(k));
				b[i*n2*n3 + j*n3 + k] = double(k);//double(i);//2*f(double(i), double(j), double(k));
				c[i*n2*n3 + j*n3 + k] = 2*double(k);//double(i + i);//3*f(double(i), double(j), double(k));
			}
		}
	}

	double **A, **B;

	MKL_INT r11, r12, r13, r21, r22, r23;

	A = compress(n1, n2, n3, a, eps, r11, r12, r13);
	B = compress(n1, n2, n3, b, eps, r21, r22, r23);

	double *C[4];
	C[0] = new double[(r11+r21) * (r21+r22) * (r13+r23)];
	C[1] = new double[n1 * (r11+r21)];
	C[2] = new double[n2 * (r12+r22)];
	C[3] = new double[n3 * (r13+r23)];

	sum(n1, n2, n3,
			r11, r12, r13, A,
			r21, r22, r23, B,
			C);

//	full(n1, n2, n3, (r11+r21), (r12+r22), (r13+r23), C, c_);

	cout << "r11 = " << r11 << ", r12 = " << r12 << ", r13 = " << r13 << endl;
	cout << "r21 = " << r21 << ", r22 = " << r22 << ", r23 = " << r23 << endl;

	print_matrix("exact", n1*n2, n3, c, n3);
	print_matrix("sum", n1*n2, n3, c_, n3);

	delete [] a;
	delete [] a_;
	delete [] b;
	delete [] c;
	delete [] c_;
	for (int i = 0; i < 4; ++i) {
		delete [] A[i];
		delete [] B[i];
		delete [] C[i];
	}
	delete [] A;
	delete [] B;

	return 0;
}

void sum(MKL_INT n1, MKL_INT n2, MKL_INT n3,
		MKL_INT r11, MKL_INT r12, MKL_INT r13, double **tuc1,
		MKL_INT r21, MKL_INT r22, MKL_INT r23, double **tuc2,
		double **tuc3) {

	for (int i = 0; i < r11; ++i) {
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r12, r13, tuc1[0]+i*r12*r13, r13, tuc3[0]+
				i*(r12+r22)*(r13+r23), (r13+r23));
	}
	for (int i = 0; i < r21; ++i) {
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r22, r23, tuc2[0]+i*r22*r23, r23, tuc3[0]+
				r11*(r12+r22)*(r13+r23) + r12*(r13+r23) + r13 + i*(r12+r22)*(r13+r23), (r13+r23));
	}

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1, r11, tuc1[1], r11, tuc3[1], (r11+r21));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1, r21, tuc2[1], r21, tuc3[1]+r11, (r11+r21));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2, r12, tuc1[2], r12, tuc3[2], (r12+r22));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2, r22, tuc2[2], r22, tuc3[2]+r12, (r12+r22));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3, r13, tuc1[3], r13, tuc3[3], (r13+r23));
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3, r23, tuc2[3], r23, tuc3[3]+r13, (r13+r23));

}

void full(MKL_INT n1, MKL_INT n2, MKL_INT n3, MKL_INT r1, MKL_INT r2, MKL_INT r3, double **tuc, double *res) {

	const double alpha = 1.0;
	const double beta = 0.0;
	//TODO: if r > n
	double *z1 = new double[n1*n2*n3];
	double *z2 = new double[n1*n2*n3];

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n1, r2*r3, r1, alpha, tuc[1], r1, tuc[0], r2*r3, beta, z1, r2*r3);

	mkl_dimatcopy ('R', 'T', n1, r2*r3, alpha, z1, r2*r3, n1);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n2, r3*n1, r2, alpha, tuc[2], r2, z1, r3*n1, beta, z2, r3*n1);

	mkl_dimatcopy ('R', 'T', n2, r3*n1, alpha, z2, r3*n1, n2);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n3, n1*n2, r3, alpha, tuc[3], r3, z2, n1*n2, beta, res, n1*n2);

	mkl_dimatcopy ('R', 'T', n3, n1*n2, alpha, res, n1*n2, n3);

	delete [] z1;
	delete [] z2;

}

double **compress(MKL_INT n1, MKL_INT n2, MKL_INT n3, double *a, double eps, MKL_INT &r1, MKL_INT &r2, MKL_INT &r3) {

	double **result = new double *[4];

	const double alpha = 1.0;
	const double beta = 0.0;

	double *q1, *q2, *q3;

	double *z1 = new double[n1*n2*n3];
	double *z2 = new double[n1*n2*n3];
	double *z3 = new double[n1*n2*n3];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, z1, 1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, z2, 1);
	mkl_dimatcopy ('R', 'T', n1, n2*n3, alpha, z2, n2*n3, n1);
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*n2*n3, 1, a, 1, z3, 1);
	mkl_dimatcopy ('R', 'T', n1*n2, n3, alpha, z3, n3, n1*n2);

	q1 = svd_trunc(n1, (n2 * n3), z1, eps, r1);
	q2 = svd_trunc(n2, (n1 * n3), z2, eps, r2);
	q3 = svd_trunc(n3, (n1 * n2), z3, eps, r3);

	double *h = new double[r1*r2*r3]();
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n1*n2, r3, n3, alpha, a, n3, q3, r3, beta, z3, r3);
	mkl_dimatcopy ('R', 'T', n1*n2, r3, alpha, z3, r3, n1*n2);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							r3*n1, r2, n2, alpha, z3, n2, q2, r2, beta, z2, r2);
	mkl_dimatcopy ('R', 'T', r3*n1, r2, alpha, z2, r2, r3*n1);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							r2*r3, r1, n1, alpha, z2, n1, q1, r1, beta, h, r1);
	mkl_dimatcopy ('R', 'T', r2*r3, r1, alpha, h, r1, r2*r3);

	result[0] = h;
	result[1] = q1;
	result[2] = q2;
	result[3] = q3;

	delete[] z1;
	delete[] z2;
	delete[] z3;

	return result;
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

void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
		MKL_INT i, j;
		printf( "\n %s\n", desc );
		for( i = 0; i < m; i++ ) {
				for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
				printf( "\n" );
		}
}
