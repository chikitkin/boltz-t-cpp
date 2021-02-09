#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
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
#include "tensor_class.h"
#include "mesh.h"
using namespace std;

const double PI = acos(-1.0); // TODO

class VelocityGrid {
public:
	// Constructor
	VelocityGrid(int nvx_, int nvy_, int nvz_, double *vx__, double *vy__, double *vz__); // TODO
//	~VelocityGrid(); // TODO

	double *vx_;
	double *vy_;
	double *vz_;

	int nvx, nvy, nvz;
	int nv;

	double hvx, hvy, hvz;
	double hv3;

	double *vx;
	double *vy;
	double *vz;

	double *zerox;
	double *zeroy;
	double *zeroz;

	double *onesx;
	double *onesy;
	double *onesz;

	Tensor vx_t;
	Tensor vy_t;
	Tensor vz_t;

	Tensor v2;

	Tensor zero;
	Tensor ones;
};

double *f_maxwell(VelocityGrid v,
		double n, double ux, double uy, double uz,
		double T, double Rg);

Tensor f_maxwell_t(VelocityGrid v,
		double n, double ux, double uy, double uz,
		double T, double Rg);

class GasParams {
public:
	// Constructor
	GasParams(double Mol_ = 40e-3, double Pr_ = 2.0 / 3.0, double g_ = 5.0 / 3.0, double d_ = 3.418e-10,
			double C_ = 144.4, double T_0_ = 273.11, double mu_0_ = 2.125e-5);

    double Na = 6.02214129e+23; // Avogadro constant
    double kB = 1.381e-23; // Boltzmann constant, J / K
    double Ru = 8.3144598; // Universal gas constant

    double Mol; // = Mol
    double Rg; // = self.Ru  / self.Mol  # J / (kg * K)
    double m; // # kg

	double g; // # specific heat ratio
	double d; // # diameter of molecule

    double Pr;

    double C;
	double T_0;
	double mu_0;

	double mu_suth(double T) const;

	double mu(double T) const;
};

class Problem {
public:
	// Constructor
	Problem();
	vector < string > bc_type_list;
	vector < Tensor > bc_data;
	Tensor f_init(double x, double y, double z, VelocityGrid v);
};

Tensor set_bc(const GasParams& gas_params,
		const string& bc_type, const Tensor& bc_data, const Tensor& f, const VelocityGrid& v,
		const Tensor& vn, const Tensor& vnp, const Tensor& vnm, double tol);

vector <double> comp_macro_params(const Tensor& f, const VelocityGrid& v, const GasParams& gas_params);
Tensor comp_j(const vector <double> params, const Tensor& f, const VelocityGrid& v, const GasParams& gas_params);

class Solution {
public:
	GasParams gas_params;

	Mesh mesh;
	// velocity mesh parameters
	int nv;
	double vmax;
	double tol;
	// path to the folder
	string path;
	// number of iteration
	int it;

	double tau;
	double t;

	double *vx_;

	double *vx, *vy, *vz;

	vector < double * > vn;

	Tensor vx_t, vy_t, vz_t;

	vector < Tensor > vn_t;
	vector < Tensor > vnp_t;
	vector < Tensor > vnm_t;
	vector < Tensor > vn_abs_t;

	vector < Tensor > f;

	vector < Tensor > fp;
	vector < Tensor > fm;

	vector < Tensor > flux;

	vector < Tensor > rhs;

	double *n;
	double *ux, *uy, *uz;
	double *T;


	// Constructor
	Solution();
	// Destructor
	~Solution();

	void load_t();
	void save_t();

	void maketimesteps(double CFL, int nt, ...);
};

#endif /* SOLVER_H_ */
