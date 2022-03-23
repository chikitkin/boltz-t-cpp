#include "mesh.h"
#include "header.h"
#include "solver.h"

#include "full.h"
#include "tucker.h"

template <class Tensor>
double *f_maxwell(std::shared_ptr < VelocityGrid<Tensor> > v,
		double n, double ux, double uy, double uz,
		double T, double Rg)
{
	double *fmax = new double[v->nv];

	double C = n * pow(1.0 / (2.0 * PI * Rg * T), 1.5); // TODO pi
	double s = 2.0 * Rg * T;

	for (int i = 0; i < v->nv; ++i) {
		fmax[i] = C * exp(-((pow(v->vx[i] - ux, 2) + pow(v->vy[i] - uy, 2) + pow(v->vz[i] - uz, 2)) / s));
	}

	return fmax;
}

template <class Tensor>
Tensor f_maxwell_t(std::shared_ptr < VelocityGrid<Tensor> > v,
		double n, double ux, double uy, double uz,
		double T, double Rg)
{
	double C = n * pow(1.0 / (2.0 * PI * Rg * T), 1.5);
	double s = 2.0 * Rg * T;

	double *u1 = new double[v->nvx];
	for (int i = 0; i < v->nvx; ++i) {
		u1[i] = exp(-(pow(v->vx_[i] - ux, 2.0) / s));
	}

	double *u2 = new double[v->nvy];
	for (int i = 0; i < v->nvy; ++i) {
		u2[i] = exp(-(pow(v->vy_[i] - uy, 2.0) / s));
	}

	double *u3 = new double[v->nvz];
	for (int i = 0; i < v->nvz; ++i) {
		u3[i] = exp(-(pow(v->vz_[i] - uz, 2.0) / s));
	}

	Tensor fmax(v->nvx, v->nvy, v->nvz, u1, u2, u3);

	delete [] u1;
	delete [] u2;
	delete [] u3;

	return C * fmax;
}

double GasParams::mu_suth(double T) const {
	return mu_0 * ((T_0 + C) / (T + C)) * (pow(T / T_0, 3.0 / 2.0));
}

double GasParams::mu(double T) const {
	return mu_suth(200.0) * (pow(T / 200.0, 0.734));
}

template <class Tensor>
Tensor Problem<Tensor>::getInit(std::shared_ptr < GasParams > gas_params, std::shared_ptr < VelocityGrid<Tensor> > v,
		double x, double y, double z,
		const std::vector<Tensor>& initData) {
	if (x < 0.0) {
		return initData[0];
	}
	else {
		return initData[1];
	}
}

template <class Tensor>
Tensor Problem<Tensor>::getBC(std::shared_ptr < GasParams > gas_params, std::shared_ptr < VelocityGrid<Tensor> > v,
		double x, double y, double z,
		bcType bc_type, const Tensor& bc_data,
		const Tensor& f,
		const Tensor& vn, const Tensor& vnp, const Tensor& vnm,
		double tol)
{
	switch (bc_type) // TODO use enum
	{
	case SYMMETRYX:
		return reflect(f, 'X');
	case SYMMETRYY:
		return reflect(f, 'Y');
	case SYMMETRYZ:
		return reflect(f, 'Z');
	case INLET:
		return Tensor(bc_data);
	case OUTLET:
		return Tensor(bc_data);
	case WALL:
	{
		double Ni = v->hv3 * (round_t(f * vnp, tol, 1e+6)).sum();
		double Nr = v->hv3 * (round_t(bc_data * vnm, tol, 1e+6)).sum();
		return -(Ni / Nr) * bc_data;
	}
	default:
		std::cout << "Wrong BC type." << std::endl;
		exit(-1);
	}
}

template <class Tensor>
VelocityGrid<Tensor>::VelocityGrid(int nvx_, int nvy_, int nvz_, double *vx__, double *vy__, double *vz__)
: nvx(nvx_), nvy(nvy_), nvz(nvz_)
{
	nv = nvx * nvy * nvz;

	vx_ = new double[nvx];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', nvx, 1, vx__, 1, vx_, 1);
	vy_ = new double[nvy];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', nvy, 1, vy__, 1, vy_, 1);
	vz_ = new double[nvz];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', nvz, 1, vz__, 1, vz_, 1);

	hvx = vx_[1] - vx_[0];
	hvy = vy_[1] - vy_[0];
	hvz = vz_[1] - vz_[0];
	hv3 = hvx * hvy * hvz; // TODO better designations

	vx = new double[nv];
	vy = new double[nv];
	vz = new double[nv];

	for (int i = 0; i < nvx; ++i) {
		for (int j = 0; j < nvy; ++j) {
			for (int k = 0; k < nvz; ++k) {
				vx[i * nvy * nvz + j * nvz + k] = vx_[i];
				vy[i * nvy * nvz + j * nvz + k] = vy_[j];
				vz[i * nvy * nvz + j * nvz + k] = vz_[k];
			}
		}
	}

	zerox = new double[nvx]();
	onesx = new double[nvx];
	for (int i = 0; i < nvx; ++i) {
		onesx[i] = 1.0;
	}

	zeroy = new double[nvy]();
	onesy = new double[nvy];
	for (int i = 0; i < nvy; ++i) {
		onesy[i] = 1.0;
	}

	zeroz = new double[nvz]();
	onesz = new double[nvz];
	for (int i = 0; i < nvz; ++i) {
		onesz[i] = 1.0;
	}

	vx_t = Tensor(nvx, nvy, nvz, vx_, onesy, onesz);
	vy_t = Tensor(nvx, nvy, nvz, onesx, vy_, onesz);
	vz_t = Tensor(nvx, nvy, nvz, onesx, onesy, vz_);

	v2 = vx_t * vx_t + vy_t * vy_t + vz_t * vz_t;
	v2.round(1e-7);

	zero = Tensor(nvx, nvy, nvz, zerox, zeroy, zeroz);
	ones = Tensor(nvx, nvy, nvz, onesx, onesy, onesz);
	
	double* vn_abs_r1_tmp = new double[nv];
	for (int i = 0; i < nv; ++i) {
		vn_abs_r1_tmp[i] = pow(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i], 0.5);
	}
	vn_abs_r1 = Tensor(nvx, nvy, nvz, vn_abs_r1_tmp);
	delete [] vn_abs_r1_tmp;
	vn_abs_r1.round(1e-14, 1);
}

template <class Tensor>
VelocityGrid<Tensor>::~VelocityGrid()
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

template <class Tensor>
std::vector <double> comp_macro_params(const Tensor& f, std::shared_ptr < VelocityGrid<Tensor> > v, std::shared_ptr < GasParams > gas_params)
{
	double n = v->hv3 * f.sum();

	if (n <= 0.0) {
		n = 1e+10;
	}

	double ux = (1.0 / n) * v->hv3 * (v->vx_t * f).sum();
	double uy = (1.0 / n) * v->hv3 * (v->vy_t * f).sum();
	double uz = (1.0 / n) * v->hv3 * (v->vz_t * f).sum();

	double u2 = pow(ux, 2.0) + pow(uy, 2.0) + pow(uz, 2.0);

	double T = (1.0 / (3.0 * n * gas_params->Rg)) * (v->hv3 * (v->v2 * f).sum() - n * u2);

	if (T <= 0.0) {
		T = 1.0;
	}

	double rho = gas_params->m * n;
	double p = rho * gas_params->Rg * T;
	double mu = gas_params->mu(T);
	double nu = p / mu;

	return {n, ux, uy, uz, T, rho, p, nu};
}

template <class Tensor>
Tensor comp_j(const std::vector <double>& params, const Tensor& f, std::shared_ptr < VelocityGrid<Tensor> > v, std::shared_ptr < GasParams > gas_params)
{
	double n = params[0];
	double ux = params[1];
	double uy = params[2];
	double uz = params[3];
	double T = params[4];
	double rho = params[5];
	double p = params[6];
	double nu = params[7];

	double *tmp = new double[std::max({v->nvx, v->nvy, v->nvz})];

	for (int i = 0; i < v->nvx; ++i) {
		tmp[i] = (1.0 / pow(2.0 * gas_params->Rg * T, 0.5)) * (v->vx_[i] - ux);
	}
	Tensor cx(v->nvx, v->nvy, v->nvz, tmp, v->onesy, v->onesz);

	for (int i = 0; i < v->nvy; ++i) {
		tmp[i] = (1.0 / pow(2.0 * gas_params->Rg * T, 0.5)) * (v->vy_[i] - uy);
	}
	Tensor cy(v->nvx, v->nvy, v->nvz, v->onesx, tmp, v->onesz);

	for (int i = 0; i < v->nvz; ++i) {
		tmp[i] = (1.0 / pow(2.0 * gas_params->Rg * T, 0.5)) * (v->vz_[i] - uz);
	}
	Tensor cz(v->nvx, v->nvy, v->nvz, v->onesx, v->onesy, tmp);

	Tensor c2 = ((cx * cx) + (cy * cy) + (cz * cz));
	c2.round(1e-7);

	double Sx = (1.0 / n) * v->hv3 * (cx * c2 * f).sum();
	double Sy = (1.0 / n) * v->hv3 * (cy * c2 * f).sum();
	double Sz = (1.0 / n) * v->hv3 * (cz * c2 * f).sum();

	Tensor fmax = f_maxwell_t(v, n, ux, uy, uz, T, gas_params->Rg);

	Tensor f_plus = fmax * (v->ones + ((4.0 / 5.0) * (1.0 - gas_params->Pr) * (Sx*cx + Sy*cy + Sz*cz) * ((c2 + (- 5.0 / 2.0) * v->ones))));
	Tensor J = nu * (f_plus - f);
	J.round(1e-7);

	delete [] tmp;

	return J;
}

template <class Tensor>
Solution<Tensor>::Solution(
		std::shared_ptr < GasParams > gas_params,
		std::shared_ptr < Mesh > mesh,
		std::shared_ptr < VelocityGrid<Tensor> > v,
		std::shared_ptr < Problem<Tensor> > problem,
		std::shared_ptr < Config > config
		)
: gas_params(gas_params),
  mesh(mesh),
  v(v),
  problem(problem),
  config(config)
{
//	path = "./" + "job_tuck_" + config.solver_type + '_' + '/'; // TODO datetime.now().strftime("%Y.%m.%d_%H:%M:%S") + '/';
	path = "./";
	// TODO make directory

	vn.resize(mesh->nFaces, Tensor());
	double* vn_tmp = new double[v->nv];
	vnm.resize(mesh->nFaces, Tensor());
	double* vnm_tmp = new double[v->nv];
	vnp.resize(mesh->nFaces, Tensor());
	double* vnp_tmp = new double[v->nv];
	vn_abs.resize(mesh->nFaces, Tensor());
	double* vn_abs_tmp = new double[v->nv];

	for (int jf = 0; jf < mesh->nFaces; ++jf) {
		for (int i = 0; i < v->nv; ++i) {
			vn_tmp[i] = 
					mesh->faceNormals[jf][0] * v->vx[i] +
					mesh->faceNormals[jf][1] * v->vy[i] +
					mesh->faceNormals[jf][2] * v->vz[i];
			if (vn_tmp[i] <= 0.0) {
				vnm_tmp[i] = vn_tmp[i];
				vnp_tmp[i] = 0.0;
			}
			else {
				vnm_tmp[i] = 0.0;
				vnp_tmp[i] = vn_tmp[i];
			}
			vn_abs_tmp[i] = vnp_tmp[i] - vnm_tmp[i];
		}
		vn[jf] = Tensor(v->nvx, v->nvy, v->nvz, vn_tmp, 1e-3);
		vnm[jf] = Tensor(v->nvx, v->nvy, v->nvz, vnm_tmp, config->tol);
		vnp[jf] = Tensor(v->nvx, v->nvy, v->nvz, vnp_tmp, config->tol);
		vn_abs[jf] = Tensor(v->nvx, v->nvy, v->nvz, vn_abs_tmp);
		vn_abs[jf].round(1e-14, 6);
	}
	h = *std::min_element(mesh->cellDiameters.begin(), mesh->cellDiameters.end());
	tau = h * config->CFL / (*std::max_element(v->vx_, v->vx_ + v->nvx) * (pow(3.0, 0.5)));

	diag.resize(mesh->nCells, Tensor());
	diag_r1.resize(mesh->nCells, Tensor());
	
	double *diag_tmp = new double [v->nv];
	double diag_sc;
	double *diag_t_full;
	double *ratio = new double [v->nv];

	double *zero = new double [v->nv]();

	for (int ic = 0; ic < mesh->nCells; ++ic) {

		LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', v->nv, 1, zero, 1, diag_tmp, 1);
		diag_sc = 0.0;

		for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
			int jf = mesh->cellFaces[ic][j];

			for (int i = 0; i < v->nv; ++i) {
				vn_tmp[i] = mesh->getOutSign(ic, j) * (
				mesh->faceNormals[jf][0] * v->vx[i] +
				mesh->faceNormals[jf][1] * v->vy[i] +
				mesh->faceNormals[jf][2] * v->vz[i]);
				if (vn_tmp[i] <= 0.0) {
					vnm_tmp[i] = vn_tmp[i];
					vnp_tmp[i] = 0.0;
				}
				else {
					vnm_tmp[i] = 0.0;
					vnp_tmp[i] = vn_tmp[i];
				}
				vn_abs_tmp[i] = vnp_tmp[i] - vnm_tmp[i];

				diag_tmp[i] += (mesh->faceAreas[jf] / mesh->cellVolumes[ic]) * vnp_tmp[i];
			}
			diag_sc += 0.5 * (mesh->faceAreas[jf] / mesh->cellVolumes[ic]);
		}
		
		diag_r1[ic] = diag_sc * v->vn_abs_r1;
		diag_t_full = diag_r1[ic].full();

		for (int i = 0; i < v->nv; ++i) {
			ratio[i] = diag_t_full[i] / diag_tmp[i];
		}
		
		diag_r1[ic] = (1.0 / *std::min_element(ratio, ratio + v->nv)) * diag_r1[ic];
		
		delete [] diag_t_full;
	}

	f.resize(mesh->nCells, Tensor());

	if (config->initType == "default") {
		double x;
		double y;
		double z;
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			x = mesh->cellCenters[ic][0];
			y = mesh->cellCenters[ic][1];
			z = mesh->cellCenters[ic][2];
			f[ic] = problem->getInit(gas_params, v,
					x, y, z,
					problem->initData);
		}
	}
	// TODO: other inits

	fLeftRight.resize(mesh->nFaces, std::vector<Tensor>{Tensor(), Tensor()});
	flux.resize(mesh->nFaces, Tensor());
	rhs.resize(mesh->nCells, Tensor());
	df.resize(mesh->nCells, Tensor());

	n.resize(mesh->nCells, 0.0);
	rho.resize(mesh->nCells, 0.0);
	ux.resize(mesh->nCells, 0.0);
	uy.resize(mesh->nCells, 0.0);
	uz.resize(mesh->nCells, 0.0);
	p.resize(mesh->nCells, 0.0);
	T.resize(mesh->nCells, 0.0);
	nu.resize(mesh->nCells, 0.0);
	rank.resize(mesh->nCells, 0.0);
	data.resize(mesh->nCells, std::vector < double >());

	it = 0;
	//  TODO: create_res

	delete [] vn_tmp;
	delete [] vnm_tmp;
	delete [] vnp_tmp;
	delete [] vn_abs_tmp;

	delete [] zero;
	delete [] diag_tmp;

}

template <class Tensor>
void Solution<Tensor>::make_time_steps(std::shared_ptr<Config> config, int nt)
{
	tau = h * config->CFL / (*std::max_element(v->vx_, v->vx_ + v->nvx) * (pow(3.0, 0.5)));

	Tensor J;
	
	Tensor vnm_loc;
	Tensor div_tmp;
	Tensor incr;

	it = 0;

	while(it < nt) {
		std::cout << "Step " << it << std::endl;
		it += 1;
		// reconstruction for inner faces
		// 1st order
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			for (int j = 0; j < mesh->getNumFacesOfCell(ic); j++) {
				int jf = mesh->getFacesOfCell(ic, j);
				// 0 if outer, 1 else
				fLeftRight[jf][1 - mesh->getOutIndex(ic, j)] = f[ic];
			}
		}
		// boundary condition
		// loop over all boundary faces
		for (int ibc = 0; ibc < problem->bcTypes.size(); ++ibc) {
			int tag = mesh->boundaryFaceTagsSet[ibc];
			std::vector<int> BCFaces = mesh->getFacesForTag(tag);
			bcType type = problem->bcTypes[ibc];
			Tensor &data = problem->bcData[ibc];
			for (const int &jf: BCFaces) {
				double x = mesh->faceCenters[jf][0];
				double y = mesh->faceCenters[jf][1];
				double z = mesh->faceCenters[jf][2];
				fLeftRight[jf][mesh->getOutIndex(jf)] = problem->getBC(
						gas_params, v,
						x, y, z,
						type, data,
						fLeftRight[jf][1 - mesh->getOutIndex(jf)],
						mesh->getOutSign(jf) * vn[jf],
						(1 - mesh->getOutIndex(jf)) * vnp[jf] - mesh->getOutIndex(jf) * vnm[jf],
						(1 - mesh->getOutIndex(jf)) * vnm[jf] - mesh->getOutIndex(jf) * vnp[jf], // TODO fix
						config->tol
				);
			}
		}
		// Riemann solver - compute fluxes
		for (int jf = 0; jf < mesh->nFaces; ++jf) {
			flux[jf] = 0.5 * mesh->faceAreas[jf] *
					((fLeftRight[jf][0] + fLeftRight[jf][1]) * vn[jf] - (fLeftRight[jf][1] - fLeftRight[jf][0]) * vn_abs[jf]);
			flux[jf].round(config->tol);
		}
		// computation of the right hand side
		std::vector < double > params(8, 0.0); // for J
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			rhs[ic] = v->zero;
			// sum up fluxes from all faces of this cell
			for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
				int jf = mesh->cellFaces[ic][j];
				rhs[ic] = rhs[ic] - (mesh->getOutSign(ic, j) / mesh->cellVolumes[ic]) * flux[jf];
				rhs[ic].round(config->tol);
			}
			// compute macroparameters and collision integral
			params = comp_macro_params(f[ic], v, gas_params);

			n[ic] = params[0];
			ux[ic] = params[1];
			uy[ic] = params[2];
			uz[ic] = params[3];
			T[ic] = params[4];
			rho[ic] = params[5];
			p[ic] = params[6];
			nu[ic] = params[7];

			double compression =
					static_cast<double>(f[ic].r()[0] * f[ic].r()[1] * f[ic].r()[2] +
					f[ic].r()[0] * f[ic].n()[0] +
					f[ic].r()[1] * f[ic].n()[1] +
					f[ic].r()[2] * f[ic].n()[2]) / static_cast<double>(f[ic].n()[0] * f[ic].n()[1] * f[ic].n()[2]);

			data[ic] = {
					n[ic],
					ux[ic],
					uy[ic],
					uz[ic],
					T[ic],
					compression};

			J = comp_j(params, f[ic], v, gas_params);
			rhs[ic] = rhs[ic] + J;
			rhs[ic].round(config->tol);
		}
		double frob_norm = 0.0;
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			frob_norm += pow(rhs[ic].norm(), 2.0);
		}
		frob_norm = pow(frob_norm / mesh->nCells, 0.5);
		frob_norm_iter.push_back(frob_norm);

		std::cout << "Frob norm = " << frob_norm << std::endl;

		if (!config->isImplicit) {
			// Update values
			for (int ic = 0; ic < mesh->nCells; ++ic) {
				f[ic] = f[ic] + tau * rhs[ic];
				f[ic].round(config->tol);
			}
		}

		else {
			for (int ic = mesh->nCells - 1; ic >= 0; --ic) {
				df[ic] = rhs[ic];
			}
			// Backward sweep
			for (int ic = mesh->nCells - 1; ic >= 0; --ic) {
				// loop over neighbors of cell ic
				for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
					int jf = mesh->cellFaces[ic][j];
					int icn = mesh->cellNeighbors[ic][j]; // index of neighbor
					if ((icn >= 0) && (icn > ic)) {
						vnm_loc = 0.5 * (-v->vn_abs_r1 + mesh->getOutSign(ic, j) * vn[jf]); // vnm[jf] or -vnp[jf]
						df[ic] = df[ic] - (mesh->faceAreas[jf] / mesh->cellVolumes[ic]) * vnm_loc * df[icn];
						df[ic].round(config->tol);
					}
				}
				// divide by diagonal coefficient
				div_tmp = ((1.0 / tau + nu[ic]) * v->ones + diag_r1[ic]);
				div_tmp.round(1e-3, 1);
				df[ic] = df[ic] / div_tmp;
				df[ic].round(config->tol);
			}
			// Forward sweep
			for (int ic = 0; ic < mesh->nCells; ++ic) {
				// loop over neighbors of cell ic
				incr = v->zero;
				for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
					int jf = mesh->cellFaces[ic][j];
					int icn = mesh->cellNeighbors[ic][j]; // index of neighbor, -1 if no neighbor
					if ((icn >= 0) && (icn < ic)) {
						vnm_loc = 0.5 * (-v->vn_abs_r1 + mesh->getOutSign(ic, j) * vn[jf]); // vnm[jf] or -vnp[jf]
						incr = incr - (mesh->faceAreas[jf] / mesh->cellVolumes[ic]) * vnm_loc * df[icn];
						incr.round(config->tol);
					}
				}
				// divide by diagonal coefficient
				div_tmp = ((1.0 / tau + nu[ic]) * v->ones + diag_r1[ic]);
				div_tmp.round(1e-3, 1);
				df[ic] = df[ic] + (incr / div_tmp);
				df[ic].round(config->tol);
			}
			// Update values
			for (int ic = 0; ic < mesh->nCells; ++ic) {
				f[ic] = f[ic] + df[ic];
				f[ic].round(config->tol);
			}
		}
		// TODO save
	}
	mesh->write_tecplot(data, "tec.dat",
			{"n", "ux", "uy", "uz", "T", "comp"});
}

template class Full f_maxwell_t<Full>(std::shared_ptr< VelocityGrid<Full> >, double, double, double, double, double, double);
template class std::vector<double> comp_macro_params<Full>(Full const&, std::shared_ptr< VelocityGrid<Full> >, std::shared_ptr<GasParams>);
template class VelocityGrid<Full>;
template class Problem<Full>;
template class Solution<Full>;

template class Tucker f_maxwell_t<Tucker>(std::shared_ptr< VelocityGrid<Tucker> >, double, double, double, double, double, double);
template class std::vector<double> comp_macro_params<Tucker>(Tucker const&, std::shared_ptr< VelocityGrid<Tucker> >, std::shared_ptr<GasParams>);
template class VelocityGrid<Tucker>;
template class Problem<Tucker>;
template class Solution<Tucker>;
