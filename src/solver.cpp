#include "mesh.h"
#include "header.h"
#include "solver.h"

#include "full.h"
#include "tucker.h"

std::vector<int> getParallelRanges(const std::vector<double>& times, int numThreads) {

	int n = times.size();

	double totalTime = std::accumulate(times.begin(), times.end(), 0.0);
	std::vector<int> threadRanges(numThreads + 1);
	threadRanges[0] = 0;
	int thread = 0;
	double sum = 0.0;
	for (int i = 0; i < n; ++i) {
		sum += times[i];
		if (sum > totalTime / numThreads) {
			threadRanges[thread + 1] = i;
			++thread;
			sum = 0.0;
		}
	}
	threadRanges[numThreads] = n;

	return threadRanges;
}

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
Tensor Problem<Tensor>::getInit(double x, double y, double z,
		const std::vector<Tensor>& initData) {
	if (x < 0.0) {
		return initData[0];
	}
	else {
		return initData[1];
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
std::vector <double> Solution<Tensor>::comp_macro_params(const Tensor& f)
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
Tensor Solution<Tensor>::comp_j(const std::vector <double>& params, const Tensor& f)
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
void Solution<Tensor>::comp_M()
{
	std::ofstream file;
	file.open("M.txt", std::ofstream::trunc);

	for (int ibf = 0; ibf < bcList.size(); ++ibf) {
		if (bcList[ibf]->type == WALL) {
			int jf = bcList[ibf]->jf;
			double x = mesh->faceCenters[jf][0];
			double y = mesh->faceCenters[jf][1];
			double z = mesh->faceCenters[jf][2];
			Tensor fWall = fLeftRight[jf][mesh->getOutIndex(jf)];
			double Mx = 0.5 * v->hv3 * (v->vx_t * v->v2 * fWall).sum();
			double My = 0.5 * v->hv3 * (v->vy_t * v->v2 * fWall).sum();
			double Mz = 0.5 * v->hv3 * (v->vz_t * v->v2 * fWall).sum();
			file << x << " " << Mx << " " << My << " " << Mz << "\n";
		}
	}
	file.close();
}

template <class Tensor>
void Solution<Tensor>::plot_residual() {
	// mglGraph gr;// create canvas
	// mglData d; d.Set(frob_norm_iter);  // convert to internal format
	// gr.Plot(d);   // plot it
	// gr.Axis();    // draw axis if you need
	// gr.WritePNG("res.png"); // save it
}

template <class Tensor>
void Solution<Tensor>::create_res() {
	std::ofstream file;
	file.open("res.txt", std::ofstream::trunc);
	file.close();
}

template <class Tensor>
void Solution<Tensor>::update_res(double frob_norm) {
	std::ofstream file;
	file.open("res.txt", std::ofstream::app);
	file << frob_norm << "\n";
	file.close();
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
	vnm.resize(mesh->nFaces, Tensor());
	vnp.resize(mesh->nFaces, Tensor());
	vn_abs.resize(mesh->nFaces, Tensor());

	#pragma omp parallel for schedule(dynamic)
	for (int jf = 0; jf < mesh->nFaces; ++jf) {
		// TODO why is it so slow?
		double* vn_tmp = new double[v->nv];
		double* vnm_tmp = new double[v->nv];
		double* vnp_tmp = new double[v->nv];
		double* vn_abs_tmp = new double[v->nv];
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

		delete [] vn_tmp;
		delete [] vnm_tmp;
		delete [] vnp_tmp;
		delete [] vn_abs_tmp;
	}
	h = *std::min_element(mesh->cellDiameters.begin(), mesh->cellDiameters.end());
	tau = h * config->CFL / (*std::max_element(v->vx_, v->vx_ + v->nvx) * (pow(3.0, 0.5)));

	diag.resize(mesh->nCells, Tensor());
	diag_r1.resize(mesh->nCells, Tensor());

	#pragma omp parallel for schedule(dynamic)
	for (int ic = 0; ic < mesh->nCells; ++ic) {

		double *diag_tmp = new double [v->nv]();
		double diag_sc = 0.0;

		for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
			int jf = mesh->cellFaces[ic][j];

			double* vn_tmp = new double[v->nv];
			double* vnm_tmp = new double[v->nv];
			double* vnp_tmp = new double[v->nv];
			double* vn_abs_tmp = new double[v->nv];
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
			delete [] vn_tmp;
			delete [] vnm_tmp;
			delete [] vnp_tmp;
			delete [] vn_abs_tmp;
		}
		
		diag_r1[ic] = diag_sc * v->vn_abs_r1;
		double *diag_t_full = diag_r1[ic].full();

		double *ratio = new double [v->nv];
		for (int i = 0; i < v->nv; ++i) {
			ratio[i] = diag_t_full[i] / diag_tmp[i];
		}
		
		diag_r1[ic] = (1.0 / *std::min_element(ratio, ratio + v->nv)) * diag_r1[ic];
		
		delete [] diag_tmp;
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
			f[ic] = problem->getInit(x, y, z,
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
	comp.resize(mesh->nCells, 0.0);
	rank.resize(mesh->nCells, 0.0);
	data.resize(mesh->nCells, std::vector < double >());

	bcList.reserve(mesh->nBoundaryFaces);

	for (int ibc = 0; ibc < problem->bcTags.size(); ++ibc) {
		int tag = problem->bcTags[ibc];
		bcType type = problem->bcTypes[ibc];
		Tensor &data = problem->bcData[ibc]; // TODO make shared
		if (mesh->boundaryFacesForEachTag.count(tag)) {
			std::vector<int> bcFaces = mesh->boundaryFacesForEachTag[tag];
			for (const int &jf: bcFaces) {
				if (type == SYMMETRYX) {
					BCSYMMETRYX<Tensor> * pBoundaryCondition = new BCSYMMETRYX<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == SYMMETRYY) {
					BCSYMMETRYY<Tensor> * pBoundaryCondition = new BCSYMMETRYY<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == SYMMETRYZ) {
					BCSYMMETRYZ<Tensor> * pBoundaryCondition = new BCSYMMETRYZ<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == INLET) {
					BCINLET<Tensor> * pBoundaryCondition = new BCINLET<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == OUTLET) {
					BCOUTLET<Tensor> * pBoundaryCondition = new BCOUTLET<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == WALL) {
					BCWALL<Tensor> * pBoundaryCondition = new BCWALL<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
			}
		}
	}
				
	it = 0;
	create_res();

	mesh->divideMesh(omp_get_max_threads());
}

template <class Tensor>
void Solution<Tensor>::make_time_steps(std::shared_ptr<Config> config, int nt)
{
	tau = h * config->CFL / (*std::max_element(v->vx_, v->vx_ + v->nvx) * (pow(3.0, 0.5)));

	it = 0;


	int numThreads = omp_get_max_threads();
	std::cout << "Number of threads is: " << numThreads << std::endl;
	std::vector<double> timesForFaceFluxes(mesh->nFaces, 1.0);
	std::vector<double> timesForCellsRHS(mesh->nCells, 1.0);
	std::vector<double> timesForCellsUpdate(mesh->nCells, 1.0);

	while(it < nt) {
		it += 1;
		// reconstruction for inner faces
		// 1st order
		auto t0 = omp_get_wtime();
		#pragma omp parallel for schedule(dynamic)
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			for (int j = 0; j < mesh->cellFaces[ic].size(); j++) {
				int jf = mesh->cellFaces[ic][j];
				// 0 if outer, 1 else
				fLeftRight[jf][1 - mesh->getOutIndex(ic, j)] = f[ic];
			}
		}
		auto t1 = omp_get_wtime();
		// boundary condition
		// loop over all boundary faces
		#pragma omp parallel for schedule(dynamic)
		for (int ibf = 0; ibf < bcList.size(); ++ibf) {
			int jf = bcList[ibf]->jf;
			double x = mesh->faceCenters[jf][0];
			double y = mesh->faceCenters[jf][1];
			double z = mesh->faceCenters[jf][2];
			fLeftRight[jf][mesh->getOutIndex(jf)] = bcList[ibf]->applyBC(
					x, y, z,
					fLeftRight[jf][1 - mesh->getOutIndex(jf)],
					mesh->getOutSign(jf) * vn[jf],
					vn_abs[jf],
					config->tol
			);			
		}

		auto t2 = omp_get_wtime();

		// Compute ranges for each thread
		//std::vector<int> threadRanges = getParallelRanges(timesForFaceFluxes, numThreads);

		// Riemann solver - compute fluxes
		#pragma omp parallel for schedule(dynamic)
		for (int jf = 0; jf < mesh->nFaces; ++jf) {
			double begin = omp_get_wtime();
			flux[jf] = 0.5 * mesh->faceAreas[jf] *
					((fLeftRight[jf][0] + fLeftRight[jf][1]) * vn[jf] - (fLeftRight[jf][1] - fLeftRight[jf][0]) * vn_abs[jf]);
			flux[jf].round(config->tol);
			double end = omp_get_wtime();
			timesForFaceFluxes[jf] = end - begin;
		}
		auto t3 = omp_get_wtime();
		// Compute ranges for each thread
		//threadRanges = getParallelRanges(timesForCellsRHS, numThreads);
		// computation of the right hand side
		#pragma omp parallel for schedule(dynamic)
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			double begin = omp_get_wtime();
			rhs[ic] = v->zero;
			// sum up fluxes from all faces of this cell
			for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
				int jf = mesh->cellFaces[ic][j];
				rhs[ic] = rhs[ic] - (mesh->getOutSign(ic, j) / mesh->cellVolumes[ic]) * flux[jf];
				rhs[ic].round(config->tol);
			}
			// compute macroparameters and collision integral
			std::vector<double> params = comp_macro_params(f[ic]);
			Tensor J = comp_j(params, f[ic]);
			rhs[ic] = rhs[ic] + J;
			rhs[ic].round(config->tol);
			double end = omp_get_wtime();
			timesForCellsRHS[ic] = end - begin;

			n[ic] = params[0];
			ux[ic] = params[1];
			uy[ic] = params[2];
			uz[ic] = params[3];
			T[ic] = params[4];
			rho[ic] = params[5];
			p[ic] = params[6];
			nu[ic] = params[7];
			comp[ic] = f[ic].compression();

			data[ic] = {
					n[ic],
					ux[ic],
					uy[ic],
					uz[ic],
					T[ic],
					comp[ic]
			};
		}

		double frob_norm = 0.0;
		#pragma omp parallel for reduction(+:frob_norm)
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			frob_norm += pow(rhs[ic].norm(), 2.0);
		}
		frob_norm = pow(frob_norm / mesh->nCells, 0.5);
		frob_norm_iter.push_back(frob_norm);
		update_res(frob_norm);

		auto t4 = omp_get_wtime();
		if (!config->isImplicit) {
			// Update values
			#pragma omp parallel for schedule(dynamic)
			for (int ic = 0; ic < mesh->nCells; ++ic) {
				double begin = omp_get_wtime();
				f[ic] = f[ic] + tau * rhs[ic];
				f[ic].round(config->tol);
				double end = omp_get_wtime();
				timesForCellsUpdate[ic] = end - begin;
			}
		}
		else {
			std::cout << "implicit" << std::endl;
			#pragma omp parallel for schedule(dynamic)
			for (int ic = 0; ic < mesh->nCells; ++ic) {
				df[ic] = rhs[ic];
			}
			// Backward sweep
			#pragma omp parallel
			{
			int partition = omp_get_thread_num();
			for (int color = mesh->nColors - 1; color >= 0; --color) {
				for (int i = mesh->C[partition][color].size() - 1; i >= 0; --i) {
					int ic = mesh->C[partition][color][i];
					int ic_perm = mesh->iPerm[ic];
					Tensor vnm_loc;
					Tensor div_tmp;
					// loop over neighbors of cell ic
					for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
						int jf = mesh->cellFaces[ic][j];
						int icn = mesh->cellNeighbors[ic][j]; // index of neighbor
						int icn_perm = mesh->iPerm[icn];
						if ((icn >= 0) && (icn_perm > ic_perm)) {
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
				#pragma omp barrier
			}

			// Forward sweep
			for (int color = 0; color < mesh->nColors; ++color) {
				for (int i = 0; i < mesh->C[partition][color].size(); ++i) {
					int ic = mesh->C[partition][color][i];
					int ic_perm = mesh->iPerm[ic];
					Tensor vnm_loc;
					Tensor incr = v->zero;
					Tensor div_tmp;
					// loop over neighbors of cell ic
					for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
						int jf = mesh->cellFaces[ic][j];
						int icn = mesh->cellNeighbors[ic][j]; // index of neighbor, -1 if no neighbor
						int icn_perm = mesh->iPerm[icn];
						if ((icn >= 0) && (icn_perm < ic_perm)) {
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
				#pragma omp barrier
			}
			}
			// Update values
			#pragma omp parallel for schedule(dynamic)
			for (int ic = 0; ic < mesh->nCells; ++ic) {
				f[ic] = f[ic] + df[ic];
				f[ic].round(config->tol);
			}
		}
		auto t5 = omp_get_wtime();
		// TODO save
		timings[RECONSTRUCTION].push_back(t1 - t0);
		timings[BOUNDARY_CONDITIONS].push_back(t2 - t1);
		timings[FLUXES].push_back(t3 - t2);
		timings[RHS].push_back(t4 - t3);
		timings[UPDATE].push_back(t5 - t4);

		if (it % 10 == 0) {
			mesh->write_tecplot(data, "tec.dat",
					{"n", "ux", "uy", "uz", "T", "comp"});
		}
	}
	mesh->write_tecplot(data, "tec.dat",
			{"n", "ux", "uy", "uz", "T", "comp"});
}

template class VelocityGrid<Full>;
template class Problem<Full>;
template class Solution<Full>;
template Full f_maxwell_t(std::shared_ptr < VelocityGrid<Full> >, double, double, double, double, double, double);

template class VelocityGrid<Tucker>;
template class Problem<Tucker>;
template class Solution<Tucker>;
template Tucker f_maxwell_t(std::shared_ptr < VelocityGrid<Tucker> >, double, double, double, double, double, double);
