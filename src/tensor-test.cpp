#include "header.h"
#include "full.h"
#include "tucker.h"
#include "solver.h"

void print_matrix( char* desc, int m, int n, double* a, int lda )
{
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i * lda + j] );
		printf( "\n" );
	}
}

double f(double x, double y, double z) {
	return 2.0 + x + 4*y;
}

double g(double x, double y, double z) {
	return 3.0 + exp(x);
}

int main(){
	typedef Tucker Tensor;
	
// vn_abs_r1 test
	
	Tensor vn_abs_r1;

	int nvx = 44;
	int nvy = 44;
	int nvz = 44;
	int nv = nvx * nvy * nvz;
	
	double Ru = 8.3144598; // Universal gas constant
	double Mol = 40e-3; // = Mol
	double Rg = Ru / Mol;
	
	double T_s = 200.0;
	
	double v_s = pow(2. * Rg * T_s, 0.5);
	double vmax = 22.0 * v_s;

	double hv = 2.0 * vmax / nvx;

	double *vx__ = new double[nvx];

	for (int i = 0; i < nvx; ++i) {
		vx__[i] = - vmax + (hv / 2.0) + i * hv;
	}
	std::cout << vmax << std::endl;
	std::cout << hv << std::endl;
	
	double *vx = new double[nv];
	double *vy = new double[nv];
	double *vz = new double[nv];

	double *vx_ = new double[nvx];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', nvx, 1, vx__, 1, vx_, 1);
	double *vy_ = new double[nvy];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', nvy, 1, vx__, 1, vy_, 1);
	double *vz_ = new double[nvz];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', nvz, 1, vx__, 1, vz_, 1);

	for (int i = 0; i < nvx; ++i) {
		for (int j = 0; j < nvy; ++j) {
			for (int k = 0; k < nvz; ++k) {
				vx[i * nvy * nvz + j * nvz + k] = vx_[i];
				vy[i * nvy * nvz + j * nvz + k] = vx_[j];
				vz[i * nvy * nvz + j * nvz + k] = vx_[k];
			}
		}
	}
	
	double* vn_abs_r1_tmp = new double[nv];
	double sum = 0;
	for (int i = 0; i < nv; ++i) {
		vn_abs_r1_tmp[i] = pow(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i], 0.5);
		sum += vn_abs_r1_tmp[i] * vn_abs_r1_tmp[i];
	}
	sum = pow(sum, 0.5);
	std::cout << "sum " << sum << std::endl;

	vn_abs_r1 = Tensor(nvx, nvy, nvz, vn_abs_r1_tmp);
	delete [] vn_abs_r1_tmp;

	std::cout << "vn_abs_r1" << std::endl;
	std::cout << vn_abs_r1 << std::endl;
	std::cout << vn_abs_r1.norm() << std::endl;

	vn_abs_r1.round(1e-7, 1);

	std::cout << vn_abs_r1 << std::endl;
	std::cout << vn_abs_r1.norm() << std::endl;
	
/*************************************************************************************/
// Reflect test
/*
	int n1, n2, n3;
	n1 = 44;
	n2 = 44;
	n3 = 44;

	double eps = 1e-5;
	
	double *a1 = new double[n1*n2*n3];
	double *a2 = new double[n1*n2*n3];
	double *a3 = new double[n1*n2*n3];

	double s;
	double s_res;

	double a_norm;
	double diff_norm;
	double diff;

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = double(rand());
				a2[(n1 - i - 1)*n2*n3 + j*n3 + k] = a1[i*n2*n3 + j*n3 + k];
			}
		}
	}

	Tensor t_t(n1, n2, n3, a1);

	Tensor t_r_t(n1, n2, n3, a2);

	Tensor t_t_r = reflect(t_t, 'X');

	cout << t_r_t;
	cout << t_t_r;

	Tensor error = t_r_t - t_t_r;

	double *error_full = error.full();

//	for (int i = 0; i < n1*n2*n3; ++i) {
//		cout << a1[i] << " ";
//	}
//	cout << endl;


	cout << "reflect " << (t_r_t - t_t_r).norm() / t_r_t.norm() << endl;




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
//	cout << t1_1.r()[0] << t1_1.r()[1] << t1_1.r()[2] << '\n';
	Tensor t2_1(n1, n2, n3, a2, eps);
	Tensor t3_1;
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
/*
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
	Tensor t3_2;
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
/*
	s = 0.0;
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = 2.0;
				a2[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a3[i*n2*n3 + j*n3 + k] = a1[i*n2*n3 + j*n3 + k] / a2[i*n2*n3 + j*n3 + k];
				s += a3[i*n2*n3 + j*n3 + k];
			}
		}
	}
	Tensor t1_3(n1, n2, n3, a1, eps);
	Tensor t2_3(n1, n2, n3, a2, eps);
	t2_3.round(1e-14, 1);
	Tensor t3_3;
	t3_3 = t1_3 / t2_3;
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
/*
	s = 0.0;
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a2[i*n2*n3 + j*n3 + k] = g(double(i), double(j), double(k));
				a3[i*n2*n3 + j*n3 + k] = 2 * f(double(i), double(j), double(k));
				s += a3[i*n2*n3 + j*n3 + k];
			}
		}
	}
	Tensor t1_4(n1, n2, n3, a1, eps);
	Tensor t2_4(n1, n2, n3, a2, eps);
	Tensor t3_4;
	t3_4 = 2.0 * t1_4;
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
/*
	s = 0.0;
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a2[i*n2*n3 + j*n3 + k] = g(double(i), double(j), double(k));
				a3[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k)) + g(double(i), double(j), double(k)) + f(double(i), double(j), double(k));
				s += a3[i*n2*n3 + j*n3 + k];
			}
		}
	}
	Tensor t1_5(n1, n2, n3, a1, eps);
	Tensor t2_5(n1, n2, n3, a2, eps);
	Tensor t3_5;
	t3_5 = t1_5;
	t3_5 = t3_5 + t3_5;
	t3_5 = t3_5 + t2_5;
	t3_5.round(eps);

	s_res = t3_5.sum();

	cout << "s_diff sum = " << (s_res - s) / s << endl;

	double *t3_5_full;
	t3_5_full = t3_5.full();

//	print_matrix("t3_full", n1*n2, n3, t3_full, n3);
//	print_matrix("a3", n1*n2, n3, a3, n3);
	a_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	cblas_daxpy(n1*n2*n3, -1.0, t3_5_full, 1, a3, 1);
	diff_norm = cblas_dnrm2(n1*n2*n3, a3, 1);
	diff = diff_norm / a_norm;
	cout << "diff sum = " << diff << endl;

	delete [] t3_5_full;
/*************************************************************************************/
/*
	delete [] a1;
	delete [] a2;
	delete [] a3;
/*************************************************************************************/
	return 0;
}

/*************************************************************************************/
// QR test
/*
	double* A = new double[n1*n2];
	double* B = new double[n1*n2]();

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			A[i*n2 + j] = double(i*n2 + j + 1);
		}
	}

	int size = min(n1, n2);

	double *tau = new double[size];

	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'U', n1, n2, A, n2, B, n2);

	print_matrix("B", n1, n2, B, n2);

	print_matrix("A", n1, n2, A, n2);

	double** QR = qr(n1, n2, A);

	print_matrix("Q", n1, size, QR[0], size);
	print_matrix("R", size, n2, QR[1], n2);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			n1, n2, size, 1.0, QR[0], size, QR[1], n2, 0.0, A, n2);

	print_matrix("A result", n1, n2, A, n2);

	cout << endl;

	delete [] A;
	delete [] B;

/*************************************************************************************/
/*************************************************************************************/
// Constructor test
/*
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a1[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
			}
		}
	}

	Tensor t(n1, n2, n3, a1, 1e-7);

	cout << t.r()[0] << " " << t.r()[1] << " " << t.r()[2] << endl;

	Tensor t_trunc = round_t(t, 1e-14, 1);

	cout << t_trunc.r()[0] << " " << t_trunc.r()[1] << " " << t_trunc.r()[2] << endl;

	cout << (t - t_trunc).norm() / t.norm() << endl;

	cout << t_trunc << endl;

/*************************************************************************************/
// Transpose test
/*
	double *a_0 = new double[n1*n2*n3];
	double *a = new double[n1*n2*n3];
	double *a_true = new double[n1*n2*n3];

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k) {
				a_0[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a[i*n2*n3 + j*n3 + k] = f(double(i), double(j), double(k));
				a_true[k*n2*n1 + j*n1 + i] = a[i*n2*n3 + j*n3 + k];
			}
		}
	}

	print_matrix("A", 1, n1*n2*n3, a, n1*n2*n3);

	transpose(n1, n2, n3, a, 210);

	print_matrix("A 0", 1, n1*n2*n3, a_0, n1*n2*n3);
	print_matrix("A", 1, n1*n2*n3, a, n1*n2*n3);
	print_matrix("A true", 1, n1*n2*n3, a_true, n1*n2*n3);

	delete [] a_0;
	delete [] a;
	delete [] a_true;
/*************************************************************************************/

