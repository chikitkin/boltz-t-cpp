#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"

#define min(a,b) ((a)>(b)?(b):(a))

#define M 6
#define N 5

extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern double *svd_trunc(MKL_INT m, MKL_INT n, double *a, double eps, MKL_INT &r);

int main() {

	MKL_INT m = M, n = N;

	MKL_INT r;

	double eps = 0.1;

	/*
	double *a = new double[m*n]{
		8.79,  9.93,  9.83, 5.45,  3.16,
		6.11,  6.91,  5.04, -0.27,  7.98,
	   -9.15, -7.93,  4.86, 4.85,  3.01,
		9.57,  1.64,  8.83, 0.74,  5.80,
	   -3.49,  4.02,  9.80, 10.00,  4.27,
		9.84,  0.15, -8.99, -6.02, -5.31
	};
	*/
	double *a = new double[m*n]{
		1,  2,  3, 4,  5,
		1,  2,  3, 4,  5,
		1,  2,  3, 4,  5,
		1,  2,  3, 4,  5,
		1,  2,  3, 4,  5,
		1,  2,  3, 4,  5,
	};

	double *u;

	print_matrix( "initial", m, n, a, n );

	u = svd_trunc(m, n, a, eps, r);

	printf("\n r = %d \n", r);

	print_matrix( "u trunc", m, r, u, r );

	delete[] a;
	delete[] u;

	return 0;

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

