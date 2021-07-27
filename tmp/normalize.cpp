Tensor normalize() {

		const double alpha = 1.0;
		const double beta = 0.0;

		Tensor t = Tensor(n1, n2, n3, min(n1, r1), min(n2, r2), min(n3, r3));

		double **QR1, **QR2, **QR3;
		QR1 = qr(n1, r1, mode1);
		QR2 = qr(n2, r2, mode2);
		QR3 = qr(n3, r3, mode3);

		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1, min(n1, r1), QR1[0], min(n1, r1), t.mode1, min(n1, r1));
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2, min(n2, r2), QR2[0], min(n2, r2), t.mode2, min(n2, r2));
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3, min(n3, r3), QR3[0], min(n3, r3), t.mode3, min(n3, r3));

		double *z1 = new double[r1*r2*r3];
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r1*r2*r3, 1, core, 1, z1, 1);
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
				min(n1, r1)*min(n2, r2), min(n3, r3), r3, alpha, z3, r3, QR3[1], min(n3, r3), beta, t.core, min(n3, r3));

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

		return t;

	}




void round(double eps) {

		const double alpha = 1.0;
		const double beta = 0.0;

		MKL_INT r1_, r2_, r3_;

		Tensor t = Tensor(n1, n2, n3, r1, r2, r3);
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n1*r1, 1, mode1, 1, t.mode1, 1);
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n2*r2, 1, mode2, 1, t.mode2, 1);
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', n3*r3, 1, mode3, 1, t.mode3, 1);
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r1*r2*r3, 1, core, 1, t.core, 1);

		t.normalize();

		double **A;
		A = compress(t.r1, t.r2, t.r3, t.core, eps, r1_, r2_, r3_);

		delete [] core;
		core = new double[r1_ * r2_ * r3_];
		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', r1_ * r2_ * r3_, 1, A[0], 1, core, 1);

		delete [] mode1;
		mode1 = new double[n1 * r1_];
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
										n1, r1_, t.r1, alpha, t.mode1, t.r1, A[1], r1_, beta, mode1, r1_);
		delete [] mode2;
		mode2 = new double[n2 * r2_];
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
										n2, r2_, t.r2, alpha, t.mode2, t.r2, A[2], r2_, beta, mode2, r2_);
		delete [] mode3;
		mode3 = new double[n3 * r3_];
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
										n3, r3_, t.r3, alpha, t.mode3, t.r3, A[3], r3_, beta, mode3, r3_);

		r1 = r1_;
		r2 = r2_;
		r3 = r3_;

		for (int i = 0; i < 4; ++i) {
			delete [] A[i];
		}
		delete [] A;

	}
