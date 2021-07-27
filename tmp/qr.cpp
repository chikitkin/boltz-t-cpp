#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>

using namespace std;

#define N1 4
#define N2 10

extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

int main() {

	MKL_INT n1 = N1, n2 = N2;

	const double alpha = 1.0;
	const double beta = 0.0;

	double eps = 0.01;

	double *a = new double[n1*n2];

	delete [] a;

	a = new double[n1*n2];

	int size = min(n1, n2);

	double *tau = new double[size];

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			a[i*n2 + j] = double(j+1);
		}
	}

	print_matrix("a", n1, n2, a, n2);

	LAPACKE_dgeqrf (LAPACK_ROW_MAJOR, n1, n2, a, n2, tau);

	print_matrix("qr", n1, n2, a, n2);
	print_matrix("tau", size, 1, tau, 1);

	double *r = new double[size*n2]();

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'U', size, n2, a, n2, r, n2);

	print_matrix("r", size, n2, r, n2);

	LAPACKE_dorgqr (LAPACK_ROW_MAJOR, n1, size, size, a, n2, tau);

	print_matrix("q", n1, size, a, n2);

	double *a_ = new double[n1*n2];

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
								n1, n2, size, alpha, a, n2, r, n2, beta, a_, n2);

	print_matrix("a_", n1, n2, a_, n2);

	delete [] a;
	delete [] tau;
	delete [] r;
	delete [] a_;

}

void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
		MKL_INT i, j;
		printf( "\n %s\n", desc );
		for( i = 0; i < m; i++ ) {
				for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
				printf( "\n" );
		}
}
