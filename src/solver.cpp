#include <iostream>
#include <fstream>
#include <sstream>
// #include <filesystem>
#include <string>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>
#include <numeric>
#include "mesh.h"
#include "solver.h"

#include "tensor.h"
using namespace std;

VelocityGrid::VelocityGrid(int nvx_, int nvy_, int nvz_, double *vx__, double *vy__, double *vz__)
: nvx(nvx_), nvy(nvy_), nvz(nvz_)
{
	nv = nvx * nvy * nvz;

	vx_ = new double[nvx];
	vy_ = new double[nvy];
	vz_ = new double[nvz];

	hvx = vx_[1] - vx_[0];
	hvy = vy_[1] - vy_[0];
	hvz = vz_[1] - vz_[0];
	hv3 = hvx * hvy * hvz; // TODO better designations

	// TODO vx vy vz
	vx = new double[nv];
	vy = new double[nv];
	vz = new double[nv];

	zerox = new double[nvx];
	onesx = new double[nvx];
	for (int i = 0; i < nvx; ++i) {
		zerox[i] = 0.0;
		onesx[i] = 1.0;
	}

	zeroy = new double[nvy];
	onesy = new double[nvy];
	for (int i = 0; i < nvy; ++i) {
		zeroy[i] = 0.0;
		onesy[i] = 1.0;
	}

	zeroz = new double[nvz];
	onesz = new double[nvz];
	for (int i = 0; i < nvz; ++i) {
		zeroz[i] = 0.0;
		onesz[i] = 1.0;
	}

	vx_t = Tensor(nvx, nvy, nvz, vx_, onesy, onesz);
	vy_t = Tensor(nvx, nvy, nvz, onesx, vy_, onesz);
	vz_t = Tensor(nvx, nvy, nvz, onesx, onesy, vz_);

	v2 = vx_t * vx_t + vy_t * vy_t + vz_t * vz_t;
	v2.round(1e-7);

	zero = Tensor(nvx, nvy, nvz, zerox, zeroy, zeroz);
	ones = Tensor(nvx, nvy, nvz, onesx, onesy, onesz);
}

VelocityGrid::~VelocityGrid()
{
	delete [] zerox;
	delete [] zeroy;
	delete [] zeroz;

	delete [] onesx;
	delete [] onesy;
	delete [] onesz;

	delete [] vx_;
	delete [] vy_;
	delete [] vz_;

	delete [] vx;
	delete [] vy;
	delete [] vz;
}
// TODO copy or create??
double *f_maxwell(VelocityGrid v,
		double n, double ux, double uy, double uz,
		double T, double Rg)
{
	double *fmax = new double[v.nv];

	double C = n * pow(1.0 / (2.0 * PI * Rg * T), 1.5); // TODO pi
	double s = 2.0 * Rg * T;

	for (int i = 0; i < v.nv; ++i) {
		fmax[i] = C * exp(-((pow(v.vx[i] - ux, 2) + pow(v.vy[i] - uy, 2) + pow(v.vz[i] - uz, 2)) / s));
	}

	return fmax;
}

Tensor f_maxwell_t(VelocityGrid v,
		double n, double ux, double uy, double uz,
		double T, double Rg)
{
	double C = n * pow(1.0 / (2.0 * PI * Rg * T), 1.5);
	double s = 2.0 * Rg * T;

	double *u1 = new double[v.nvx];
	for (int i = 0; i < v.nvx; ++i) {
		u1[i] = exp(-(pow(v.vx[i] - ux, 2.0) / s));
	}

	double *u2 = new double[v.nvy];
	for (int i = 0; i < v.nvy; ++i) {
		u2[i] = exp(-(pow(v.vy[i] - uy, 2.0) / s));
	}

	double *u3 = new double[v.nvz];
	for (int i = 0; i < v.nvz; ++i) {
		u3[i] = exp(-(pow(v.vz[i] - uz, 2.0) / s));
	}

	Tensor fmax(v.nvx, v.nvy, v.nvz, u1, u2, u3);

	delete [] u1;
	delete [] u2;
	delete [] u3;

	return C * fmax;
}

GasParams::GasParams(double Mol_, double Pr_, double g_, double d_,
		double C_, double T_0_, double mu_0_)
: Mol(Mol_), Pr(Pr_), g(g_), d(d_), C(C_), T_0(T_0_), mu_0(mu_0_)
{
	Rg = Ru / Mol; // = self.Ru  / self.Mol  # J / (kg * K)
	m = Mol / Na; // # kg
}

double GasParams::mu_suth(double T) const {
	return mu_0 * ((T_0 + C) / (T + C)) * ((T / T_0) ** (3.0 / 2.0));
}

double GasParams::mu(double T) const {
	return mu_suth(200.0) * ((T / 200.0) ** (0.734));
}

Tensor set_bc(const GasParams& gas_params,
		const string& bc_type, const Tensor& bc_data, const Tensor& f, const VelocityGrid& v,
		const Tensor& vn, const Tensor& vnp, const Tensor& vnm, double tol)
{
	switch (bc_type) // TODO use enum
	{
	case "SYM_X":
		return reflect(f, 'X');
	case "SYM_Y":
		return reflect(f, 'Y');
	case "SYM_Z":
		return reflect(f, 'Z');
	case "SYM":
		return Tensor(f);
	case "IN":
		return Tensor(bc_data);
	case "OUT":
		return Tensor(bc_data);
	case "WALL":
	{
		return 0; // TODO
	}
	default:
		cout << "Wrong BC type." << endl;
		exit(-1);
	}
}

vector <double> comp_macro_params(const Tensor& f, const VelocityGrid& v, const GasParams& gas_params)
{
	double n = v.hv3 * f.sum();

	if (n <= 0.0) {
		n = 1e+10;
	}

	double ux = (1.0 / n) * v.hv3 * (v.vx_t * f).sum();
	double uy = (1.0 / n) * v.hv3 * (v.vy_t * f).sum();
	double uz = (1.0 / n) * v.hv3 * (v.vz_t * f).sum();

	double u2 = pow(ux, 2.0) + pow(uy, 2.0) + pow(uz, 2.0);

	double T = (1.0 / (3.0 * n * gas_params.Rg)) * (v.hv3 * (v.v2 * f).sum() - n * u2);

	if (T <= 0.0) {
		T = 1.0;
	}

	double rho = gas_params.m * n;
	double p = rho * gas_params.Rg * T;
	double mu = gas_params.mu(T);
	double nu = p / mu;

	return vector <double> {n, ux, uy, uz, T, rho, p, nu};
}

Tensor comp_j(const vector <double> params, const Tensor& f, const VelocityGrid& v, const GasParams& gas_params)
// TODO const v
{
	double n = params[0];
	double ux = params[1];
	double uy = params[2];
	double uz = params[3];
	double T = params[4];
	double rho = params[5];
	double p = params[6];
	double nu = params[7];

	double *tmp = new double[max(v.nvx, max(v.nvy, v.nvz))];

	for (int i = 0; i < v.nvx; ++i) {
		tmp[i] = (1.0 / pow(2.0 * gas_params.Rg * T, 0.5)) * (v.vx_[i] - ux);
	}
	Tensor cx(v.nvx, v.nvy, v.nvz, tmp, v.onesy, v.onesz);

	for (int i = 0; i < v.nvy; ++i) {
		tmp[i] = (1.0 / pow(2.0 * gas_params.Rg * T, 0.5)) * (v.vy_[i] - uy);
	}
	Tensor cy(v.nvx, v.nvy, v.nvz, v.onesx, tmp, v.onesz);

	for (int i = 0; i < v.nvz; ++i) {
		tmp[i] = (1.0 / pow(2.0 * gas_params.Rg * T, 0.5)) * (v.vz_[i] - uz);
	}
	Tensor cz(v.nvx, v.nvy, v.nvz, v.onesx, v.onesy, tmp);

	Tensor c2 = ((cx * cx) + (cy * cy) + (cz * cz));
	c2.round(1e-7);

	double Sx = (1.0 / n) * v.hv3 * (cx * c2 * f).sum();
	double Sy = (1.0 / n) * v.hv3 * (cy * c2 * f).sum();
	double Sz = (1.0 / n) * v.hv3 * (cz * c2 * f).sum();

	Tensor fmax = f_maxwell_t(v, n, ux, uy, uz, T, gas_params.Rg);

	Tensor f_plus = fmax * (v.ones + ((4.0 / 5.0) * (1.0 - gas_params.Pr) * (Sx*cx + Sy*cy + Sz*cz) * ((c2 - (5.0 / 2.0) * v.ones))));
	Tensor J = nu * (f_plus - f); // TODO
	J.round(1e-7);

	delete [] tmp;

	return J;
}

Solution::Solution()
{
	path = "../job-"; // dummy
}

int main()
{
	return 0;
}
