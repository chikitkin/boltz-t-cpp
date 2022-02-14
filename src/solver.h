#ifndef SOLVER_H_
#define SOLVER_H_

#include "read_mesh.h"
#include "header.h"

#include "full.h"
#include "tucker.h"

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
	Tensor vn_abs_r1;

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
	std::vector < Tensor > init_tensor_list;
	std::vector < int > bc_types;
	Tensor f_init(double x, double y, double z);

	std::vector < Tensor > bc_data;
	Tensor set_bc(const GasParams& gas_params, const VelocityGrid<Tensor>& v,
			int bc_type, const Tensor& bc_data,
			const Tensor& f,
			const Tensor& vn, const Tensor& vnp, const Tensor& vnm,
			double tol);
};

struct Config {
	std::string solver_type = "explicit";

	double CFL = 0.5;
	double tol = 1e-7;

	std::string init_type = "default";
	std::string init_filename = "";

	int save_tec_step = 1e+5;
	int save_macro_step = 1e+5;
};

template <class Tensor> std::vector <double> comp_macro_params(const Tensor& f, const VelocityGrid<Tensor>& v, const GasParams& gas_params);
template <class Tensor> Tensor comp_j(const std::vector <double>& params, const Tensor& f, const VelocityGrid<Tensor>& v, const GasParams& gas_params);

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
	std::string path;

	std::vector < Tensor > vn;
	std::vector < Tensor > vnm;
	std::vector < Tensor > vnp;
	std::vector < Tensor > vn_abs;
	Tensor vn_abs_r1;

	double h;
	double tau;

	std::vector < Tensor > diag;
	std::vector < Tensor > diag_r1;

	std::vector < Tensor > f;

	std::vector < Tensor > fp;
	std::vector < Tensor > fm;
	std::vector < Tensor > flux;
	std::vector < Tensor > rhs;
	std::vector < Tensor > df;

	// Arrays for macroparameters
	std::vector < double > n;
	std::vector < double > rho;
	std::vector < double > ux, uy, uz;
	std::vector < double > p;
	std::vector < double > T;
	std::vector < double > nu;
	std::vector < double > rank;
	std::vector < std::vector < double > > data;

	int it;
	std::vector <double> frob_norm_iter;

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
