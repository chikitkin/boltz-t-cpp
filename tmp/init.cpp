#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>
#include <bits/stdc++.h>
#include <chrono>

using namespace std;

#define N1 10
#define N2 10
#define N3 20

double f(double x, double y, double z) {
	return exp(-(x*x + y*y + z*z) / 1000.0) + (x / 100.0) + (x*y / 100.0) + 10.0 * sin(z) + exp(-(x*(x - 10.) + y*y + z*z) / 1000.0);
}

extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern double *svd_trunc(MKL_INT m, MKL_INT n, double *a, double eps, MKL_INT &r);
extern void compress(double eps);

int main() {

	ofstream myfile;
	myfile.open ("plot.txt");
	myfile.close();
	compress(0.01);
	/*
	for (double eps = 0.0001; eps < 0.01; eps = eps + 0.0001) {
		compress(eps);
	}

	for (double eps = 0.01; eps < 1.0; eps = eps + 0.01) {
		compress(eps);
	}
	*/
	return 0;
}

void compress(double eps) {

	MKL_INT n1 = N1, n2 = N2, n3 = N3;

	MKL_INT r1, r2, r3;

	const double alpha = 1.0;
	const double beta = 0.0;

	double *a = new double[n1*n2*n3];

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a[i*n2*n3 + j*n3 + k] = double(k);//f(double(i), double(j), double(k));
			}
		}
	}

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

	auto start = chrono::high_resolution_clock::now();
/*
	for (int I = 0; I < r1; ++I) {
		for (int J = 0; J < r2; ++J) {
			for (int K = 0; K < r3; ++K) {
				for (int i = 0; i < n1; ++i) {
					for (int j = 0; j < n2; ++j) {
						for (int k = 0; k < n3; ++k) {
							h[I*r2*r3 + J*r3 + K] += q1[i*r1 + I] * q2[j*r2 + J] * q3[k*r3 + K] * a[i*n2*n3 + j*n3 + k];
						}
					}
				}
			}
		}
	}
*/

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n1*n2, r3, n3, alpha, a, n3, q3, r3, beta, z3, r3);

	print_matrix( "z3", n1*n2, r3, z3, r3 );

	mkl_dimatcopy ('R', 'T', n1*n2, r3, alpha, z3, r3, n1*n2);

	print_matrix( "z3", r3, n1*n2, z3, n1*n2 );

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							r3*n1, r2, n2, alpha, z3, n2, q2, r2, beta, z2, r2);

	print_matrix( "z2", r3*n1, r2, z2, r2 );

	mkl_dimatcopy ('R', 'T', r3*n1, r2, alpha, z2, r2, r3*n1);

	print_matrix( "z2", r2, r3*n1, z2, r3*n1 );

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							r2*r3, r1, n1, alpha, z2, n1, q1, r1, beta, h, r1);

	print_matrix( "h", r2*r3, r1, h, r1 );

	mkl_dimatcopy ('R', 'T', r2*r3, r1, alpha, h, r1, r2*r3);

	print_matrix( "h", r1, r2*r3, h, r2*r3 );

	auto end = chrono::high_resolution_clock::now();

	double time_taken =
		  chrono::duration_cast<chrono::nanoseconds>(end - start).count();

		time_taken *= 1e-9;

		cout << "Time taken by program is : " << fixed
			 << time_taken << setprecision(9);
		cout << " sec" << endl;

	double *a_appr = new double[n1*n2*n3]();
/*
	for (int I = 0; I < r1; ++I) {
		for (int J = 0; J < r2; ++J) {
			for (int K = 0; K < r3; ++K) {
				for (int i = 0; i < n1; ++i) {
					for (int j = 0; j < n2; ++j) {
						for (int k = 0; k < n3; ++k) {
							a_appr[i*n2*n3 + j*n3 + k] += h[I*r2*r3 + J*r3 + K] * q1[i*r1 + I] * q2[j*r2 + J] * q3[k*r3 + K];
						}
					}
				}
			}
		}
	}
*/

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n1, r2*r3, r1, alpha, q1, r1, h, r2*r3, beta, z1, r2*r3);

	mkl_dimatcopy ('R', 'T', n1, r2*r3, alpha, z1, r2*r3, n1);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n2, r3*n1, r2, alpha, q2, r2, z1, r3*n1, beta, z2, r3*n1);

	mkl_dimatcopy ('R', 'T', n2, r3*n1, alpha, z2, r3*n1, n2);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							n3, n1*n2, r3, alpha, q3, r3, z2, n1*n2, beta, a_appr, n1*n2);

	mkl_dimatcopy ('R', 'T', n3, n1*n2, alpha, a_appr, n1*n2, n3);

	print_matrix( "a", n1*n2, n3, a, n3 );
	print_matrix( "a_appr", n1*n2, n3, a_appr, n3 );

	cblas_daxpy(n1*n2*n3, -1.0, a, 1, a_appr, 1);

	double a_norm = cblas_dnrm2(n1*n2*n3, a, 1);

	double diff_norm = cblas_dnrm2(n1*n2*n3, a_appr, 1);

	double diff = diff_norm / a_norm;

	print_matrix( "factor1", n1, r1, q1, r1 );
	print_matrix( "factor2", n2, r2, q2, r2 );
	print_matrix( "factor3", n3, r3, q3, r3 );

	print_matrix( "h", r2*r3, r1, h, r1 );

	printf("r1 = %d, r2 = %d, r3 = %d\n", r1, r2, r3);

	printf("eps = %f\n", eps);
	printf("diff = %f\n", diff);

	delete[] a;
	delete[] a_appr;
	delete[] h;
	delete[] q1;
	delete[] q2;
	delete[] q3;
	delete[] z1;
	delete[] z2;
	delete[] z3;

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
