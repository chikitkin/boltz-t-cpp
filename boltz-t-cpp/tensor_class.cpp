#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>
using namespace std;

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

class Tensor {
public:
	// Constructor
	Tensor(){

	}

	Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, MKL_INT r1_, MKL_INT r2_, MKL_INT r3_){
		// Create zero tensor
		n1 = n1_;
		n2 = n2_;
		n3 = n3_;
		r1 = r1_;
		r2 = r2_;
		r3 = r3_;
		// full tensor
		core = new double[r1 * r2 * r3]();
		mode1 = new double[n1 * r1];
		mode2 = new double[n2 * r2];
		mode3 = new double[n3 * r3];
	}

	Tensor(MKL_INT n1_, MKL_INT n2_, MKL_INT n3_, double *a, double eps){
		n1 = n1_;
		n2 = n2_;
		n3 = n3_;
		double **A;
		A = compress(n1, n2, n3, a, eps, r1, r2, r3);
		core = new double[r1 * r2 * r3];
		mode1 = new double[n1 * r1];
		mode2 = new double[n2 * r2];
		mode3 = new double[n3 * r3];
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r1*r2*r3, 1, A[0], 1, core, 1);
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*r1, 1, A[1], 1, mode1, 1);
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2*r2, 1, A[2], 1, mode2, 1);
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3*r3, 1, A[3], 1, mode3, 1);

		for (int i = 0; i < 4; ++i) {
				delete [] A[i];
			}
		delete [] A;
	}

	// Destructor
	~Tensor(){
		delete [] core;
		delete [] mode1;
		delete [] mode2;
		delete [] mode3;
	}
	vector<int> shape(){
		auto dim = {n1, n2, n3};
		return dim;
	}
	double At(int i1, int i2, int i3){
		// get element
		double a = 0.0;
		for(int j1 = 0; j1 < r1; j1++) {
			for(int j2 = 0; j2 < r2; j2++) {
				for(int j3 = 0; j3 < r3; j3++) {
					a += core[j1 * r2 * r3 + j2 * r3 + j3] * mode1[i1 * r1 + j1] * mode2[i2 * r2 + j2] * mode3[i3 * r3 + j3];
				}
			}
		}
		return a;
	}

	// rounding
	void round(){
		// do nothing for full tensor
	}
	// compute sum of all elements
	double sum() {
		double s = 0.0;
		for (int i = 0; i < n1 * n2 * n3; ++i){
			s += core[i];
		}
		return s;
	}

	// element-wise summation
	friend Tensor add(Tensor& t1, Tensor& t2){
		// check, that shapes are equal;
		if (t1.shape() != t2.shape()){
			cout << "Different shapes in sum!" << endl;
			exit(-1);
		}
		Tensor result(t1.n1, t1.n2, t1.n3, t1.r1+t2.r1, t1.r2+t2.r2, t1.r3+t2.r3);

		for (int i = 0; i < t1.r1; ++i) {
			LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.r2, t1.r3, t1.core+i*t1.r2*t1.r3, t1.r3, result.core+
					i*(t1.r2+t2.r2)*(t1.r3+t2.r3),(t1.r3+t2.r3));
		}
		for (int i = 0; i < t2.r1; ++i) {
			LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.r2, t2.r3, t2.core+i*t2.r2*t2.r3, t2.r3, result.core+
					t1.r1*(t1.r2+t2.r2)*(t1.r3+t2.r3) + t1.r2*(t1.r3+t2.r3) + t1.r3 + i*(t1.r2+t2.r2)*(t1.r3+t2.r3), (t1.r3+t2.r3));
		}

		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.n1, t1.r1, t1.mode1, t1.r1, result.mode1, (t1.r1+t2.r1));
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.n1, t2.r1, t2.mode1, t2.r1, result.mode1+t1.r1, (t1.r1+t2.r1));
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.n2, t1.r2, t1.mode2, t1.r2, result.mode2, (t1.r2+t2.r2));
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.n2, t2.r2, t2.mode2, t2.r2, result.mode2+t1.r2, (t1.r2+t2.r2));
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t1.n3, t1.r3, t1.mode3, t1.r3, result.mode3, (t1.r3+t2.r3));
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', t2.n3, t2.r3, t2.mode3, t2.r3, result.mode3+t1.r3, (t1.r3+t2.r3));

		return result;
	}

	//element-wise multiplication
	/*
	friend Tensor mult(Tensor& t1, Tensor& t2){
		// check, that shapes are equal;
		if (t1.shape() != t2.shape()){
			cout << "Different shapes in sum!" << endl;
			exit(-1);
		}
		vector<int> dim = t1.shape();
		Tensor result(dim[0], dim[1], dim[2], );

		for (int i = 0; i < dim[0] * dim[1] * dim[2]; ++i){
			result.core[i] = t1.core[i] * t2.core[i];
		}
		return result;
	}
	*/

	double *full() {

	    const double alpha = 1.0;
	    const double beta = 0.0;

	    double *z1 = new double[n1*n2*n3];
	    double *z2 = new double[n1*n2*n3];
	    double *res = new double[n1*n2*n3];

	    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	                            n1, r2*r3, r1, alpha, mode1, r1, core, r2*r3, beta, z1, r2*r3);

	    mkl_dimatcopy ('R', 'T', n1, r2*r3, alpha, z1, r2*r3, n1);

	    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	                            n2, r3*n1, r2, alpha, mode2, r2, z1, r3*n1, beta, z2, r3*n1);

	    mkl_dimatcopy ('R', 'T', n2, r3*n1, alpha, z2, r3*n1, n2);

	    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	                            n3, n1*n2, r3, alpha, mode3, r3, z2, n1*n2, beta, res, n1*n2);

	    mkl_dimatcopy ('R', 'T', n3, n1*n2, alpha, res, n1*n2, n3);

	    delete [] z1;
	    delete [] z2;

	    return res;

	}

private:
	int I(int i1, int i2, int i3){
		return i1 * n2 * n3 + i2 * n3 + i3;
	}
	vector<int> multiI(int I){
		// TODO: implement
	}
	// dynamic array to store parameters
	double* core;
	// sizes along each dimension;
	MKL_INT n1, n2, n3;
	// mode ranks;
	MKL_INT r1, r2, r3;
	// modes;
	double* mode1;
	double* mode2;
	double* mode3;

};

// TODO: overload operators "+, *";
Tensor operator +(Tensor& t1, Tensor& t2){
	return add(t1, t2);
}
/*
Tensor operator *(Tensor& t1, Tensor& t2){
	return mult(t1, t2);
}
*/
double f(double x, double y, double z) {
	return exp(-(x*x + y*y + z*z) / 1000.0);
}

void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}

int main(){
	MKL_INT n1, n2, n3;
	n1 = 10;
	n2 = 10;
	n3 = 10;

	double eps = 0.01;

	double *a1 = new double[n1*n2*n3];
	double *a2 = new double[n1*n2*n3];

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = double(k);
				a2[i*n2*n3 + j*n3 + k] = double(k);
			}
		}
	}

	Tensor t1 = Tensor(n1, n2, n3, a1, eps);
	Tensor t2 = Tensor(n1, n2, n3, a2, eps);
	Tensor t3 = t1 + t2;
//	t3 = t1 + t2;//add(t1, t2);

	double *t3_full;

	t3_full = t3.full();

	print_matrix("t3.full()", n1*n2, n3, t3_full, n3);

//	delete t1;
//	delete t2;
//	delete t3;
	delete [] t3_full;

	cout << "End" << endl;
	return 0;
}
