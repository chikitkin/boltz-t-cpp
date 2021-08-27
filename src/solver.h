#ifndef SOLVER_H_
#define SOLVER_H_

#include "header.h"

#include "full.h"
#include "tucker.h"
#include "mesh.h"
using namespace std;

const double PI = acos(-1.0); // TODO

class GasParams {
public:
	// Default constructor

	double Na = 6.02214129e+23; // Avogadro constant
	double kB = 1.381e-23; // Boltzmann constant, J / K
	double Ru = 8.3144598; // Universal gas constant

	double Mol = 40e-3; // = Mol
	double Rg = Ru / Mol; // = self.Ru  / self.Mol  # J / (kg * K)
	double m = Mol / Na; // # kg

	double g = 5.0 / 3.0; // # specific heat ratio
	double d = 3.418e-10; // # diameter of molecule

	double Pr = 2.0 / 3.0;

	double C = 144.4;
	double T_0 = 273.11;
	double mu_0 = 2.125e-5;

	double mu_suth(double T) const;
	double mu(double T) const;
};

template <class Tensor>
class VelocityGrid {
public:
	// Constructor
	VelocityGrid(int nvx_, int nvy_, int nvz_, double *vx__, double *vy__, double *vz__);
	// Copy constructor
	VelocityGrid(const VelocityGrid<Tensor>& v);
	// Destructor
	~VelocityGrid();

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

/*
 * Для каждого г.у. должна быть своя функция
 * И в солвере задаём г.у. на всех гранях с одним г.у., потом с другим и т.д.
 */

template <class Tensor>
class Problem {
public:
	vector < Tensor > init_tensor_list;
	Tensor f_init(double x, double y, double z);

	vector < char > bc_types;
	vector < Tensor > bc_data;
	Tensor set_bc(const GasParams& gas_params, const VelocityGrid<Tensor>& v,
			char bc_type, const Tensor& bc_data,
			const Tensor& f,
			const Tensor& vn, const Tensor& vnp, const Tensor& vnm,
			double tol);
};

struct Config {
	string solver_type = "explicit";

	double CFL = 0.5;
	double tol = 1e-7;

	string init_type = "default";
	string init_filename = "";

	int save_tec_step = 1e+5;
	int save_macro_step = 1e+5;
};

template <class Tensor> vector <double> comp_macro_params(const Tensor& f, const VelocityGrid<Tensor>& v, const GasParams& gas_params);
template <class Tensor> Tensor comp_j(const vector <double>& params, const Tensor& f, const VelocityGrid<Tensor>& v, const GasParams& gas_params);

template <class Tensor>
class Solution {
public:
	// Constructor
	Solution(
			const GasParams & gas_params_,
			const Mesh & mesh_,
			const VelocityGrid<Tensor> & v_,
			const Problem<Tensor> & problem_,
			const Config & config_
			);
	// Destructor
//	~Solution();

	GasParams gas_params;
	Mesh mesh;
	VelocityGrid<Tensor> v;

	Problem<Tensor> problem;

	Config config;

	// path to folder with output
	string path;

	vector < Tensor > vn;
	vector < Tensor > vnm;
	vector < Tensor > vnp;
	vector < Tensor > vn_abs;
	Tensor vn_abs_r1;

	double h;
	double tau;

	vector < Tensor > diag;
	vector < Tensor > diag_r1;

	vector < Tensor > f;

	vector < Tensor > fp;
	vector < Tensor > fm;
	vector < Tensor > flux;
	vector < Tensor > rhs;
	vector < Tensor > df;

	Tensor J;
	Tensor vnm_loc;
	Tensor div_tmp;
	Tensor incr;

	// Arrays for macroparameters
	vector < double > n;
	vector < double > rho;
	vector < double > ux, uy, uz;
	vector < double > p;
	vector < double > T;
	vector < double > nu;
	vector < double > rank;
	vector < vector < double > > data;

	int it;
	vector <double> frob_norm_iter;

	void create_res();
	void update_res();

	void save_tec();
	void save_macro();

	void load_restart();
	void save_restart();

	void plot_macro();

	void make_time_steps(const Config& config_, int nt);
};

template <class Tensor>
double *f_maxwell(const VelocityGrid<Tensor> & v,
		double n, double ux, double uy, double uz,
		double T, double Rg);

template <class Tensor>
Tensor f_maxwell_t(const VelocityGrid<Tensor> & v,
		double n, double ux, double uy, double uz,
		double T, double Rg);

#endif /* SOLVER_H_ */
