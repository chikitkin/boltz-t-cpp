#ifndef SOLVER_H_
#define SOLVER_H_

#include "mesh.h"
#include "header.h"

#include "full.h"
#include "tucker.h"

const double PI = acos(-1.0); // TODO

enum AlgoritmParts {
	RECONSTRUCTION,
	BOUNDARY_CONDITIONS,
	FLUXES,
	RHS,
	UPDATE
};

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
class BoundaryCondition {
public:
	int jf; // TODO maybe fix?
	std::shared_ptr < GasParams > gas_params;
	std::shared_ptr < VelocityGrid<Tensor> > v;
	Tensor bcData;

	virtual Tensor applyBC(double x, double y, double z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			double tol) = 0;
};

template <class Tensor>
class BCSYMMETRYX : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(double x, double y, double z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			double tol) override {
				return reflect(f, 'X');
	}
};

template <class Tensor>
class BCSYMMETRYY : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(double x, double y, double z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			double tol) override {
				return reflect(f, 'Y');
	}
};

template <class Tensor>
class BCSYMMETRYZ : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(double x, double y, double z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			double tol) override {
				return reflect(f, 'Z');
	}
};

template <class Tensor>
class BCINLET : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(double x, double y, double z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			double tol) override {
				return Tensor(BoundaryCondition<Tensor>::bcData);
	}
};

template <class Tensor>
class BCOUTLET : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(double x, double y, double z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			double tol) override {
				return Tensor(BoundaryCondition<Tensor>::bcData);
	}
};

template <class Tensor>
class BCWALL : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(double x, double y, double z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			double tol) override {
				double Ni = BoundaryCondition<Tensor>::v->hv3 * (round_t(0.5 * f * (vn + vn_abs), tol, 1e+6)).sum();
				double Nr = BoundaryCondition<Tensor>::v->hv3 * (round_t(0.5 * BoundaryCondition<Tensor>::bcData * (vn - vn_abs), tol, 1e+6)).sum();
				return -(Ni / Nr) * BoundaryCondition<Tensor>::bcData;
	}
};

template <class Tensor>
class Problem {
public:
	std::shared_ptr < GasParams > gas_params;
	std::shared_ptr < VelocityGrid<Tensor> > v;

	std::vector < Tensor > initData;
	Tensor getInit(double x, double y, double z,
			const std::vector<Tensor>& initData);

	std::vector < int > bcTags;
	std::vector < bcType > bcTypes;
	std::vector < Tensor > bcData;
};

struct Config {
	bool isImplicit = false;

	double CFL = 0.5;
	double tol = 1e-7;

	std::string initType = "default";
	std::string initFilename = "";

	int saveTecStep = 1e+5;
	int saveMacroStep = 1e+5;
};

template <class Tensor>
class Solution {
public:
	// Constructor
	Solution(
			std::shared_ptr < GasParams > gas_params,
			std::shared_ptr < Mesh > mesh,
			std::shared_ptr < VelocityGrid<Tensor> > v,
			std::shared_ptr < Problem<Tensor> > problem,
			std::shared_ptr < Config > config
			);
	// Destructor
//	~Solution();

	std::shared_ptr < GasParams > gas_params;
	std::shared_ptr < Mesh > mesh;
	std::shared_ptr < VelocityGrid<Tensor> > v;

	std::shared_ptr < Problem<Tensor> > problem;

	std::shared_ptr < Config > config;

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
	std::vector<std::vector<Tensor>> fLeftRight;
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
	std::vector < double > comp;
	std::vector < std::vector < double > > data;
	std::map<AlgoritmParts, std::vector<double>> timings;

	std::vector < BoundaryCondition<Tensor> *> bcList; // TODO make shared

	std::vector <double> comp_macro_params(const Tensor& f);
	Tensor comp_j(const std::vector <double>& params, const Tensor& f);

	int it;
	std::vector <double> frob_norm_iter;

	void create_res();
	void update_res();

	void save_tec();
	void save_macro();

	void load_restart();
	void save_restart();

	void plot_macro();

	void make_time_steps(std::shared_ptr<Config> config, int nt);


};

template <class Tensor>
double *f_maxwell(std::shared_ptr < VelocityGrid<Tensor> > v,
		double n, double ux, double uy, double uz,
		double T, double Rg);

template <class Tensor>
Tensor f_maxwell_t(std::shared_ptr < VelocityGrid<Tensor> > v,
		double n, double ux, double uy, double uz,
		double T, double Rg);

#endif /* SOLVER_H_ */
