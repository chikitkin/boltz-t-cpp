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

class GasParams {
public:
	// Constructor
	GasParams();
	~GasParams();

    double Na; // Avogadro constant
    double kB; // Boltzmann constant, J / K
    double Ru; // Universal gas constant

    double Mol; // = Mol
    double Rg; // = self.Ru  / self.Mol  # J / (kg * K)
    double m; // # kg

    double Pr;

    double C;
	double T_0;
	double mu_0;

	double mu_suth(double T);

	double mu(double T);

	double g; // # specific heat ratio
	double d; // # diameter of molecule
};

class Problem {
public:
	vector < string > bc_type_list;
	vector < vector < double > > bc_data;
	vector < Tensor > f_init;
};

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
