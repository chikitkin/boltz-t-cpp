#ifndef SOLVER_H_
#define SOLVER_H_

#include "mesh.h"
#include "header.h"

#include "full.h"
#include "tucker.h"

const REAL PI = acos(-1.0); // TODO

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

	REAL Na = 6.02214129e+23; // Avogadro constant
	REAL kB = 1.381e-23; // Boltzmann constant, J / K
	REAL Ru = 8.3144598; // Universal gas constant

	REAL Mol = 40e-3; // = Mol
	REAL Rg = Ru / Mol; // = self.Ru  / self.Mol  # J / (kg * K)
	REAL m = Mol / Na; // # kg

	REAL g = 5.0 / 3.0; // # specific heat ratio
	REAL d = 3.418e-10; // # diameter of molecule

	REAL Pr = 2.0 / 3.0;

	REAL C = 144.4;
	REAL T_0 = 273.11;
	REAL mu_0 = 2.125e-5;

	REAL mu_suth(REAL T) const;
	REAL mu(REAL T) const;
};

template <class Tensor>
class VelocityGrid {
public:
	// Constructor
	VelocityGrid(int nvx_, int nvy_, int nvz_, REAL *vx__, REAL *vy__, REAL *vz__);
	// Destructor
	~VelocityGrid();

	REAL *vx_;
	REAL *vy_;
	REAL *vz_;

	int nvx, nvy, nvz;
	int nv;

	REAL hvx, hvy, hvz;
	REAL hv3;

	REAL *vx;
	REAL *vy;
	REAL *vz;

	REAL *zerox;
	REAL *zeroy;
	REAL *zeroz;

	REAL *onesx;
	REAL *onesy;
	REAL *onesz;

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

	bcType type;

	virtual Tensor applyBC(REAL x, REAL y, REAL z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			REAL tol) = 0;
};

template <class Tensor>
class BCSYMMETRYX : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(REAL x, REAL y, REAL z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			REAL tol) override {
				return reflect(f, 'X');
	}
};

template <class Tensor>
class BCSYMMETRYY : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(REAL x, REAL y, REAL z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			REAL tol) override {
				return reflect(f, 'Y');
	}
};

template <class Tensor>
class BCSYMMETRYZ : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(REAL x, REAL y, REAL z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			REAL tol) override {
				return reflect(f, 'Z');
	}
};

template <class Tensor>
class BCINLET : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(REAL x, REAL y, REAL z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			REAL tol) override {
				return Tensor(BoundaryCondition<Tensor>::bcData);
	}
};

template <class Tensor>
class BCOUTLET : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(REAL x, REAL y, REAL z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			REAL tol) override {
				return Tensor(BoundaryCondition<Tensor>::bcData);
	}
};

template <class Tensor>
class BCWALL : public BoundaryCondition<Tensor> {
public:
	Tensor applyBC(REAL x, REAL y, REAL z, 
			const Tensor& f, 
			const Tensor& vn, const Tensor& vn_abs, 
			REAL tol) override {
				REAL Ni = BoundaryCondition<Tensor>::v->hv3 * (round_t(0.5 * f * (vn + vn_abs), tol, 1e+6)).sum();
				REAL Nr = BoundaryCondition<Tensor>::v->hv3 * (round_t(0.5 * BoundaryCondition<Tensor>::bcData * (vn - vn_abs), tol, 1e+6)).sum();
				return -(Ni / Nr) * BoundaryCondition<Tensor>::bcData;
	}
};

template <class Tensor>
class Problem {
public:
	std::shared_ptr < GasParams > gas_params;
	std::shared_ptr < VelocityGrid<Tensor> > v;

	std::vector < Tensor > initData;
	Tensor getInit(REAL x, REAL y, REAL z,
			const std::vector<Tensor>& initData);

	std::vector < int > bcTags;
	std::vector < bcType > bcTypes;
	std::vector < Tensor > bcData;
};

struct Config {
	bool isImplicit = false;

	REAL CFL = 0.5;
	REAL tol = 1e-7;

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

	REAL h;
	REAL tau;

	std::vector < Tensor > diag;
	std::vector < Tensor > diag_r1;

	std::vector < Tensor > f;
	std::vector < std::vector<Tensor> > fLeftRight;
	std::vector < Tensor > slope;
	std::vector < Tensor > flux;
	std::vector < Tensor > rhs;
	std::vector < Tensor > df;

	// Arrays for macroparameters
	std::vector < REAL > n;
	std::vector < REAL > rho;
	std::vector < REAL > ux, uy, uz;
	std::vector < REAL > p;
	std::vector < REAL > T;
	std::vector < REAL > nu;
	std::vector < REAL > rank;
	std::vector < REAL > compression;
	std::vector < std::vector < REAL > > data;
	std::map<AlgoritmParts, std::vector<double>> timings;

	std::vector < BoundaryCondition<Tensor> *> bcList; // TODO make shared

	std::vector <REAL> comp_macro_params(const Tensor& f);
	Tensor comp_j(const std::vector <REAL>& params, const Tensor& f);
	void write_wall_params();

	int it;
	std::vector <REAL> frob_norm_iter;

	void plot_residual();
	void create_res();
	void update_res(REAL frob_norm);

	void save_tec();
	void save_macro();

	void load_restart();
	void save_restart();

	void plot_macro();

	void make_time_steps(std::shared_ptr<Config> config, int nt);
    void reconstruction_2nd_order();

};

template <class Tensor>
REAL *f_maxwell(std::shared_ptr < VelocityGrid<Tensor> > v,
		REAL n, REAL ux, REAL uy, REAL uz,
		REAL T, REAL Rg);

template <class Tensor>
Tensor f_maxwell_t(std::shared_ptr < VelocityGrid<Tensor> > v,
		REAL n, REAL ux, REAL uy, REAL uz,
		REAL T, REAL Rg);

#endif /* SOLVER_H_ */
