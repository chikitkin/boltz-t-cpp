#include "mesh.h"
#include "header.h"
#include "solver.h"

#include "full.h"
#include "tucker.h"

REAL distance_3d(std::vector<REAL> x, std::vector<REAL> y) {
    REAL squared = pow(x[0] - y[0], 2.0) + pow(x[1] - y[1], 2.0) + pow(x[2] - y[2], 2.0);
    return pow(squared, 0.5);
}

// Power law
REAL GasParams::mu(REAL T) const {
	return pow(T, omega);
}

// Baseline to compute mu_star
REAL GasParams::mu_suth(REAL T) const {
	return mu_0 * ((T_0 + C) / (T + C)) * (pow(T / T_0, 3.0 / 2.0));
}

std::vector<int> getParallelRanges(const std::vector<double>& times, int numThreads) {

	int n = times.size();

	auto totalTime = std::accumulate(times.begin(), times.end(), 0.0);
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
REAL *f_maxwell(std::shared_ptr < VelocityGrid<Tensor> > v,
		REAL n, REAL ux, REAL uy, REAL uz,
		REAL T, REAL Rg)
{
	REAL *fmax = new REAL[v->nv];

	REAL C = n / pow((PI * T), 1.5); // TODO pi

	for (int i = 0; i < v->nv; ++i) {
		fmax[i] = C * exp( -(pow(v->vx[i] - ux, 2) + pow(v->vy[i] - uy, 2) + pow(v->vz[i] - uz, 2)) / T);
	}

	return fmax;
}

template <class Tensor>
Tensor f_maxwell_t(std::shared_ptr < VelocityGrid<Tensor> > v,
		REAL n, REAL ux, REAL uy, REAL uz,
		REAL T, REAL Rg)
{
	REAL C = n / pow((PI * T), 1.5); // TODO pi

	REAL *u1 = new REAL[v->nvx];
	for (int i = 0; i < v->nvx; ++i) {
		u1[i] = exp( -(pow(v->vx_[i] - ux, 2.0)) / T);
	}

	REAL *u2 = new REAL[v->nvy];
	for (int i = 0; i < v->nvy; ++i) {
		u2[i] = exp( -(pow(v->vy_[i] - uy, 2.0)) / T);
	}

	REAL *u3 = new REAL[v->nvz];
	for (int i = 0; i < v->nvz; ++i) {
		u3[i] = exp( -(pow(v->vz_[i] - uz, 2.0)) / T);
	}

	Tensor fmax(v->nvx, v->nvy, v->nvz, u1, u2, u3);

	delete [] u1;
	delete [] u2;
	delete [] u3;

	return C * fmax;
}

template <class Tensor>
std::vector <REAL> comp_macro_params(const Tensor& f, std::shared_ptr < VelocityGrid<Tensor> > v, std::shared_ptr < GasParams > gas_params)
{
	REAL n = v->hv3 * f.sum();

	if (n < 0.0) {
		std::cout << "n < 0" << std::endl;
		n = 1e-8;
	}

	REAL ux = (1.0 / n) * v->hv3 * (v->vx_t * f).sum();
	REAL uy = (1.0 / n) * v->hv3 * (v->vy_t * f).sum();
	REAL uz = (1.0 / n) * v->hv3 * (v->vz_t * f).sum();

	REAL u2 = pow(ux, 2.0) + pow(uy, 2.0) + pow(uz, 2.0);

	REAL T = (2.0 / (3.0 * n)) * (v->hv3 * (v->v2 * f).sum() - n * u2);

	if (T <= 0.0) {
		std::cout << "T < 0" << std::endl;
		T = 1e-8;
	}

	REAL mu = gas_params->mu(T);
	REAL nu = (8.0 / (5.0 * pow(PI, 0.5))) * (n * T / mu) / gas_params->Kn;

	// REAL rho = n;
	// REAL p = rho * T;
	// REAL Mach = pow((ux*ux + uy*uy + uz*uz) / (gas_params->g * gas_params->Rg * T), 0.5)

	return {n, ux, uy, uz, T, nu};
}

template <class Tensor>
Tensor comp_j(const std::vector <REAL>& params, const Tensor& f, REAL tol, std::shared_ptr < VelocityGrid<Tensor> > v, std::shared_ptr < GasParams > gas_params)
{
	REAL n  = params[0];
	REAL ux = params[1];
	REAL uy = params[2];
	REAL uz = params[3];
	REAL T  = params[4];
	REAL nu = params[5];

	Tensor vx = v->vx_t + (-ux) * v->ones;
	Tensor vy = v->vy_t + (-uy) * v->ones;
	Tensor vz = v->vz_t + (-uz) * v->ones;

	Tensor v2 = ((vx * vx) + (vy * vy) + (vz * vz));
	v2.round(static_cast<REAL>(1e-14));
	
	REAL qx = 0.5 * v->hv3 * (vx * v2 * f).sum();
	REAL qy = 0.5 * v->hv3 * (vy * v2 * f).sum();
	REAL qz = 0.5 * v->hv3 * (vz * v2 * f).sum();

	Tensor fmax = f_maxwell_t(v, n, ux, uy, uz, T, gas_params->Rg);

	Tensor f_plus = fmax * (v->ones + ((8.0 / 5.0) * (1.0 - gas_params->Pr) * (1.0 / (n*T*T)) * (vx*qx + vy*qy + vz*qz) * (((1.0 / T) * v2 + (- 5.0 / 2.0) * v->ones)))); // TODO round
	Tensor J = nu * (f_plus - f);
	J.round(static_cast<REAL>(tol));

	return J;
}

// TODO fix HARDCODE
template <class Tensor>
Tensor Problem<Tensor>::getInit(REAL x, REAL y, REAL z,
		const std::vector<Tensor>& initData) {
	if (x <= 0.0) {
		return initData[0];
	}
	else {
		return initData[1];
	}
}

template <class Tensor>
VelocityGrid<Tensor>::VelocityGrid(int nvx_, int nvy_, int nvz_, REAL *vx__, REAL *vy__, REAL *vz__)
: nvx(nvx_), nvy(nvy_), nvz(nvz_)
{
	nv = nvx * nvy * nvz;

	vx_ = new REAL[nvx];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', nvx, 1, vx__, 1, vx_, 1);
	vy_ = new REAL[nvy];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', nvy, 1, vy__, 1, vy_, 1);
	vz_ = new REAL[nvz];
	LAPACKE_dlacpy (LAPACK_ROW_MAJOR, 'A', nvz, 1, vz__, 1, vz_, 1);

	hvx = vx_[1] - vx_[0];
	hvy = vy_[1] - vy_[0];
	hvz = vz_[1] - vz_[0];
	hv3 = hvx * hvy * hvz; // TODO better designations

	vx = new REAL[nv];
	vy = new REAL[nv];
	vz = new REAL[nv];

	for (int i = 0; i < nvx; ++i) {
		for (int j = 0; j < nvy; ++j) {
			for (int k = 0; k < nvz; ++k) {
				vx[i * nvy * nvz + j * nvz + k] = vx_[i];
				vy[i * nvy * nvz + j * nvz + k] = vy_[j];
				vz[i * nvy * nvz + j * nvz + k] = vz_[k];
			}
		}
	}

	zerox = new REAL[nvx]();
	onesx = new REAL[nvx];
	std::fill_n(onesx, nvx, 1.0);

	zeroy = new REAL[nvy]();
	onesy = new REAL[nvy];
	std::fill_n(onesy, nvy, 1.0);

	zeroz = new REAL[nvz]();
	onesz = new REAL[nvz];
	std::fill_n(onesz, nvy, 1.0);

	vx_t = Tensor(nvx, nvy, nvz, vx_, onesy, onesz);
	vy_t = Tensor(nvx, nvy, nvz, onesx, vy_, onesz);
	vz_t = Tensor(nvx, nvy, nvz, onesx, onesy, vz_);

	v2 = vx_t * vx_t + vy_t * vy_t + vz_t * vz_t;
	v2.round(static_cast<REAL>(1e-7));

	zero = Tensor(nvx, nvy, nvz, zerox, zeroy, zeroz);
	ones = Tensor(nvx, nvy, nvz, onesx, onesy, onesz);
	
	REAL* vn_abs_r1_tmp = new REAL[nv];
	for (int i = 0; i < nv; ++i) {
		vn_abs_r1_tmp[i] = pow(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i], 0.5);
	}
	vn_abs_r1 = Tensor(nvx, nvy, nvz, vn_abs_r1_tmp);
	delete [] vn_abs_r1_tmp;
	vn_abs_r1.round(static_cast<REAL>(1e-14), 1);
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
void Solution<Tensor>::write_boundary_params() // TODO FIX FOR DIMENSIONLESS
{
	std::ofstream file;
	file.precision(17); // TODO magic number
	file.open("boundary.txt", std::ofstream::trunc);

	file << "x" << " " << "y" << " " << "z" << " ";
	file << "n" << " " << "T" << " ";
	file << "Px" << " " << "Py" << " " << "Pz" << " ";
	file << "Mx" << " " << "My" << " " << "Mz" << " ";
	file << "type" << " " << "p_inf" << " " << "S_inf";
	file << "\n";

	for (int ibf = 0; ibf < bcList.size(); ++ibf) {
		{
			int jf = bcList[ibf]->jf;
			REAL x = mesh->faceCenters[jf][0];
			REAL y = mesh->faceCenters[jf][1];
			REAL z = mesh->faceCenters[jf][2];
			file << x << " " << y << " " << z << " ";
			Tensor f = fLeftRight[jf][1 - mesh->getOutIndex(jf)];

			std::vector<REAL> params = comp_macro_params(f, v, gas_params);
			REAL n = params[0];
			REAL T = params[4];
			file << n << " " << T << " ";

			REAL Px = 2.0 * v->hv3 * (vn[jf] * v->vx_t * f).sum();
			REAL Py = 2.0 * v->hv3 * (vn[jf] * v->vy_t * f).sum();
			REAL Pz = 2.0 * v->hv3 * (vn[jf] * v->vz_t * f).sum();
			file << Px << " " << Py << " " << Pz << " ";
			REAL Mx = 0.5 * v->hv3 * (v->vx_t * v->v2 * f).sum();
			REAL My = 0.5 * v->hv3 * (v->vy_t * v->v2 * f).sum();
			REAL Mz = 0.5 * v->hv3 * (v->vz_t * v->v2 * f).sum();
			file << Mx << " " << My << " " << Mz << " ";
			file << bcList[ibf]->type << " ";
			file << gas_params->p_s << " " << gas_params->S_inf;

			file << "\n";
		}
	}
	file.close();
}

template <class Tensor>
void Solution<Tensor>::write_rate() // TODO FIX FOR DIMENSIONLESS
{
	std::ofstream file;
	file.precision(17); // TODO magic number
	file.open("rate.txt", std::ofstream::trunc);

	file << "x y z n ux uy uz T A nx ny nz" << std::endl;

	for (int jf = 0; jf < mesh->nFaces; ++jf) {
		{
			REAL x = mesh->faceCenters[jf][0];
			REAL y = mesh->faceCenters[jf][1];
			REAL z = mesh->faceCenters[jf][2];
			file << x << " " << y << " " << z << " ";

			Tensor f = 0.5 * (fLeftRight[jf][0] + fLeftRight[jf][1]);

			std::vector<REAL> params = comp_macro_params(f, v, gas_params);
			REAL n  = params[0];
			REAL ux = params[1];
			REAL uy = params[2];
			REAL uz = params[3];
			REAL T  = params[4];
			file << n << " " << ux << " " << uy << " " << uz << " " << T << " ";

			file << mesh->faceAreas[jf] << " ";
			file << mesh->faceNormals[jf][0] << " " << mesh->faceNormals[jf][1] << " " << mesh->faceNormals[jf][2];

			file << "\n";
		}
	}
	file.close();
}

template <class Tensor>
void Solution<Tensor>::write_restart()
{
	std::ofstream file;
	file.precision(17); // TODO magic number
	file.open("restart.txt", std::ofstream::trunc);

	for (int ic = 0; ic < mesh->nCells; ++ic) {
	    file << f[ic].to_string() << "\n";
	}
	file.close();
}

template <class Tensor>
void Solution<Tensor>::write_macro_restart()
{
	std::ofstream file;
	file.precision(17); // TODO magic number
	file.open("macro_restart.txt", std::ofstream::trunc);

	for (int ic = 0; ic < mesh->nCells; ++ic) {
		file << mesh->cellCenters[ic][0] << " " << mesh->cellCenters[ic][1] << " " << mesh->cellCenters[ic][2] << " ";
		file << n[ic] << " " << ux[ic] << " " << uy[ic] << " " << uz[ic] << " " << T[ic] << " " << mesh->cellVolumes[ic] << " " << compression[ic] << "\n";
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
void Solution<Tensor>::update_res(REAL frob_norm) {
	std::ofstream file;
	file.open("res.txt", std::ofstream::app);
	file << frob_norm << "\n";
	file.close();
}

template <class Tensor>
REAL Solution<Tensor>::vector_norm(std::vector<Tensor> vec) {
	REAL res = 0.0;
	#pragma omp parallel for reduction(+:res)
	for (int ic = 0; ic < mesh->nCells; ++ic) {
		res += pow(vec[ic].norm(), 2.0);
	}
	res = pow(res, 0.5); // was / mesh->nCells
	return res;
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
	std::cout << "Init started." << std::endl;
	
	std::cout << "TOL   = " << config->tol << std::endl;
	std::cout << "ORDER = " << config->order << std::endl;
	std::cout << "NVX   = " << v->nvx << std::endl;

	int numThreads = omp_get_max_threads();
	std::cout << "Number of threads is: " << numThreads << std::endl;

	vn.resize(mesh->nFaces, Tensor());
	vn_abs.resize(mesh->nFaces, Tensor());
	vn_abs_max.resize(mesh->nFaces);

	auto TIME0 = omp_get_wtime();
	if (config->vnAbsRestart != 2) {
		#pragma omp parallel for schedule(dynamic)
		for (int jf = 0; jf < mesh->nFaces; ++jf) {
			// TODO why is it so slow?
			REAL* vn_tmp = new REAL[v->nv];
			REAL* vn_abs_tmp = new REAL[v->nv];
			REAL vn_abs_max_tmp = 0.0;
			for (int i = 0; i < v->nv; ++i) {
				vn_tmp[i] = 
						mesh->faceNormals[jf][0] * v->vx[i] +
						mesh->faceNormals[jf][1] * v->vy[i] +
						mesh->faceNormals[jf][2] * v->vz[i];
				vn_abs_tmp[i] = abs(vn_tmp[i]);
				if (vn_abs_tmp[i] > vn_abs_max_tmp) {
					vn_abs_max_tmp = vn_abs_tmp[i];
				}
			}
			vn[jf] = Tensor(v->nvx, v->nvy, v->nvz, vn_tmp);
			vn[jf].round(static_cast<REAL>(1e-3));
			vn_abs[jf] = Tensor(v->nvx, v->nvy, v->nvz, vn_abs_tmp);
			vn_abs[jf].round(static_cast<REAL>(1e-14), 6);
			vn_abs_max[jf] = vn_abs_max_tmp;

			delete [] vn_tmp;
			delete [] vn_abs_tmp;
		}
	}
	auto TIME1 = omp_get_wtime();
	std::cout << "vn abs time " << TIME1 - TIME0 << " s" << std::endl;

	h = *std::min_element(mesh->cellDiameters.begin(), mesh->cellDiameters.end());
	tau = h * config->CFL / pow(pow(v->hvx, 2) + pow(v->hvy, 2) + pow(v->hvz, 2), 0.5);

	diag.resize(mesh->nCells, Tensor());
	diag_r1.resize(mesh->nCells, Tensor());

	if (config->vnAbsRestart != 2) {
		#pragma omp parallel for schedule(dynamic)
		for (int ic = 0; ic < mesh->nCells; ++ic) {

			REAL *diag_tmp = new REAL [v->nv]();
			REAL diag_sc = 0.0;

			for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
				int jf = mesh->cellFaces[ic][j];

				REAL* vn_tmp = new REAL[v->nv];
				REAL* vnp_tmp = new REAL[v->nv];
				for (int i = 0; i < v->nv; ++i) {
					vn_tmp[i] = mesh->getOutSign(ic, j) * (
					mesh->faceNormals[jf][0] * v->vx[i] +
					mesh->faceNormals[jf][1] * v->vy[i] +
					mesh->faceNormals[jf][2] * v->vz[i]);
					if (vn_tmp[i] <= 0.0) {
						vnp_tmp[i] = 0.0;
					}
					else {
						vnp_tmp[i] = vn_tmp[i];
					}
					diag_tmp[i] += (mesh->faceAreas[jf] / mesh->cellVolumes[ic]) * vnp_tmp[i];
				}
				diag_sc += 0.5 * (mesh->faceAreas[jf] / mesh->cellVolumes[ic]);
				delete [] vn_tmp;
				delete [] vnp_tmp;
			}
			
			diag_r1[ic] = diag_sc * v->vn_abs_r1;
			REAL *diag_t_full = diag_r1[ic].full();

			REAL *ratio = new REAL [v->nv];
			for (int i = 0; i < v->nv; ++i) {
				ratio[i] = diag_t_full[i] / diag_tmp[i];
			}
			
			diag_r1[ic] = (1.0 / *std::min_element(ratio, ratio + v->nv)) * diag_r1[ic];
			
			delete [] diag_tmp;
			delete [] diag_t_full;
			delete [] ratio;
		}
	}
	auto TIME2 = omp_get_wtime();
	std::cout << "diag_r1 time " << TIME2 - TIME1 << " s" << std::endl;

	if (config->vnAbsRestart == 1) {
		// SAVE to file
		std::ofstream file;
		file.precision(17); // TODO magic number
		file.open("../vn_abs_restart.txt", std::ofstream::trunc);
		for (int jf = 0; jf < mesh->nFaces; ++jf) {
			file << vn[jf].to_string() << "\n";
			file << vn_abs[jf].to_string() << "\n";
			file << vn_abs_max[jf] << "\n";
		}
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			file << diag_r1[ic].to_string() << "\n";
		}
		file.close();
	}
	if (config->vnAbsRestart == 2) {
		// READ FROM FILE
		std::ifstream file("../vn_abs_restart.txt");
		std::string line;
		for (int jf = 0; jf < mesh->nFaces; ++jf) {
			getline(file, line); vn[jf] = from_string(line, v->zero);
			getline(file, line); vn_abs[jf] = from_string(line, v->zero);
			getline(file, line); std::istringstream ss(line); ss.precision(17); ss >> vn_abs_max[jf];
		}
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			getline(file, line); diag_r1[ic] = from_string(line, v->zero);
		}
		file.close();
	}
	std::cout << "READ/SAVE TO FILE " << omp_get_wtime() - TIME2 << "s" << std::endl;

	f.resize(mesh->nCells, Tensor());

	n.   resize(mesh->nCells, 0.0);
	ux.  resize(mesh->nCells, 0.0);
	uy.  resize(mesh->nCells, 0.0);
	uz.  resize(mesh->nCells, 0.0);
	T.   resize(mesh->nCells, 0.0);
	nu.  resize(mesh->nCells, 0.0);

	compression.resize(mesh->nCells, 0.0);
	rank_x.resize(mesh->nCells, 0.0);
	rank_y.resize(mesh->nCells, 0.0);
	rank_z.resize(mesh->nCells, 0.0);
	max_rank.resize(mesh->nCells, 0.0);
	data.resize(mesh->nCells, std::vector < REAL >());

	std::cout << "f init start, ";
	if (config->initType == 0) {
		std::cout << "initType=0" << std::endl;
		REAL x;
		REAL y;
		REAL z;
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			x = mesh->cellCenters[ic][0];
			y = mesh->cellCenters[ic][1];
			z = mesh->cellCenters[ic][2];
			f[ic] = problem->getInit(x, y, z,
					problem->initData);
			std::cout << ic << " ";
		}
	}
	else if (config->initType == 1) {
		std::cout << "initType=1" << std::endl;
        std::string init_path = config->initFilename;
	    std::ifstream init(init_path);
		init.precision(17); // TODO magic number
        std::string line;
        int ic = 0;
		while (getline(init, line)) {
			f[ic] = from_string(line, v->zero);
			++ic;
			std::cout << ic << " ";
		}
		init.close();
	}
	else if (config->initType == 2) {
		std::cout << "initType=2" << std::endl;
        std::string init_path = config->initFilename;
	    std::ifstream init(init_path);
		init.precision(17); // TODO magic number
        std::string line;
        int ic = 0;
		REAL x, y, z;
        while (getline(init, line)) {
            std::istringstream line_stream(line);
			line_stream.precision(17); // TODO magic number
            line_stream >> x >> y >> z >> n[ic] >> ux[ic] >> uy[ic] >> uz[ic] >> T[ic];
            f[ic] = f_maxwell_t(v, n[ic], ux[ic], uy[ic], uz[ic], T[ic], gas_params->Rg);
            ++ic;
			std::cout << ic << " ";
        }
		init.close();
	}
	else {
	    std::cout << "Incorrect init" << std::endl;
	}
	std::cout << "\n";
	auto TIME3 = omp_get_wtime();
	std::cout << "initial time " << TIME3 - TIME2 << " s" << std::endl;

	fLeftRight.resize(mesh->nFaces, std::vector<Tensor>{Tensor(), Tensor()});
	slope.resize(mesh->nFaces, Tensor());
	flux.resize(mesh->nFaces, Tensor());
	rhs.resize(mesh->nCells, Tensor());
	df.resize(mesh->nCells, Tensor());

	bcList.reserve(mesh->nBoundaryFaces);
	auto TIME4 = omp_get_wtime();
	std::cout << "resize  time " << TIME4 - TIME3 << " s" << std::endl;

	for (int ibc = 0; ibc < problem->bcTags.size(); ++ibc) {
		int tag = problem->bcTags[ibc];
		bcType type = problem->bcTypes[ibc];
		Tensor &data_ = problem->bcData[ibc]; // TODO make shared
		if (mesh->boundaryFacesForEachTag.count(tag)) {
			std::vector<int> bcFaces = mesh->boundaryFacesForEachTag[tag];
			for (const int &jf: bcFaces) {
				if (type == SYMMETRYX) {
					BCSYMMETRYX<Tensor> * pBoundaryCondition = new BCSYMMETRYX<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data_;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == SYMMETRYY) {
					BCSYMMETRYY<Tensor> * pBoundaryCondition = new BCSYMMETRYY<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data_;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == SYMMETRYZ) {
					BCSYMMETRYZ<Tensor> * pBoundaryCondition = new BCSYMMETRYZ<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data_;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == INLET) {
					BCINLET<Tensor> * pBoundaryCondition = new BCINLET<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data_;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == OUTLET) {
					BCOUTLET<Tensor> * pBoundaryCondition = new BCOUTLET<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data_;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
				else if (type == WALL) {
					BCWALL<Tensor> * pBoundaryCondition = new BCWALL<Tensor>();
					pBoundaryCondition->jf = jf;
					pBoundaryCondition->gas_params = gas_params;
					pBoundaryCondition->v = v;
					pBoundaryCondition->bcData = data_;
					pBoundaryCondition->type = type;
					bcList.push_back(pBoundaryCondition);
				}
			}
		}
	}
	auto TIME5 = omp_get_wtime();
	std::cout << "bounds  time " << TIME5 - TIME4 << " s" << std::endl;

	create_res();

	std::cout << "Init finished." << std::endl;
}

template <class Tensor>
void Solution<Tensor>::reconstruction_2nd_order() {
    // compute slopes
    #pragma omp parallel for schedule(dynamic)
    for (int jf = 0; jf < mesh->nFaces; ++jf) {
        std::vector < int > leftRightCell = mesh->leftRightCells[jf];
        
        if ((leftRightCell[0] == -1) || (leftRightCell[1] == -1)) {
            slope[jf] = v->zero;
            continue;
        }

        std::vector < REAL > leftCellCenter = mesh->cellCenters[leftRightCell[0]];
        std::vector < REAL > rightCellCenter = mesh->cellCenters[leftRightCell[1]];
        REAL delta = distance_3d(leftCellCenter, rightCellCenter);
        slope[jf] = (1.0 / delta) * (f[leftRightCell[1]] - f[leftRightCell[0]]);
        slope[jf].round(config->tol);
    }
    #pragma omp parallel for schedule(dynamic)
    for (int ic = 0; ic < mesh->nCells; ++ic) {
        std::vector < int > hexaFaces = mesh->cellFaces[ic];
        
        Tensor slope0 = minmod(-mesh->getOutSign(ic, 0) * slope[hexaFaces[0]], mesh->getOutSign(ic, 2) * slope[hexaFaces[2]], config->tol);
        Tensor slope1 = minmod(-mesh->getOutSign(ic, 1) * slope[hexaFaces[1]], mesh->getOutSign(ic, 3) * slope[hexaFaces[3]], config->tol);
        Tensor slope2 = minmod(-mesh->getOutSign(ic, 4) * slope[hexaFaces[4]], mesh->getOutSign(ic, 5) * slope[hexaFaces[5]], config->tol);
        
        std::vector<REAL> cellCenter = mesh->cellCenters[ic];
        
        fLeftRight[hexaFaces[0]][1 - mesh->getOutIndex(ic, 0)] = round_t(f[ic] - distance_3d(cellCenter, mesh->faceCenters[hexaFaces[0]]) * slope0, config->tol, 1000000);
        fLeftRight[hexaFaces[2]][1 - mesh->getOutIndex(ic, 2)] = round_t(f[ic] + distance_3d(cellCenter, mesh->faceCenters[hexaFaces[2]]) * slope0, config->tol, 1000000);
        
        fLeftRight[hexaFaces[1]][1 - mesh->getOutIndex(ic, 1)] = round_t(f[ic] - distance_3d(cellCenter, mesh->faceCenters[hexaFaces[1]]) * slope1, config->tol, 1000000);
        fLeftRight[hexaFaces[3]][1 - mesh->getOutIndex(ic, 3)] = round_t(f[ic] + distance_3d(cellCenter, mesh->faceCenters[hexaFaces[3]]) * slope1, config->tol, 1000000);
        
        fLeftRight[hexaFaces[4]][1 - mesh->getOutIndex(ic, 4)] = round_t(f[ic] - distance_3d(cellCenter, mesh->faceCenters[hexaFaces[4]]) * slope2, config->tol, 1000000);
        fLeftRight[hexaFaces[5]][1 - mesh->getOutIndex(ic, 5)] = round_t(f[ic] + distance_3d(cellCenter, mesh->faceCenters[hexaFaces[5]]) * slope2, config->tol, 1000000);
    }
}

template <class Tensor>
void Solution<Tensor>::make_time_steps(std::shared_ptr<Config> config, int nt)
{
	tau = h * config->CFL / pow(pow(v->hvx, 2) + pow(v->hvy, 2) + pow(v->hvz, 2), 0.5);

	int numThreads = omp_get_max_threads();
	std::cout << "Number of threads is: " << numThreads << std::endl;
	std::vector<double> timesForFaceFluxes(mesh->nFaces, 1.0);
	std::vector<double> timesForCellsRHS(mesh->nCells, 1.0);
	std::vector<double> timesForCellsUpdate(mesh->nCells, 1.0);

	std::vector < double > lusgs_timings(mesh->nCells, 1.0);
	mesh->divideMesh(numThreads, lusgs_timings);

	for (int it = 0; it < nt; ++it) {

		std::cout << "Step " << it << "." << std::endl;

		// reconstruction for inner faces
		auto TIME0 = omp_get_wtime();
		// 1st order
		if (config->order == 1) {
		    #pragma omp parallel for schedule(dynamic)
		    for (int ic = 0; ic < mesh->nCells; ++ic) {
			    for (int j = 0; j < mesh->cellFaces[ic].size(); j++) {
				    int jf = mesh->cellFaces[ic][j];
				    // 0 if outer, 1 else
				    fLeftRight[jf][1 - mesh->getOutIndex(ic, j)] = f[ic];
			    }
		    }
		}
		// 2nd order
        else if (config->order == 2) {
		    reconstruction_2nd_order();
        }
        else {
            std::cout << "WRONG RECONSTRUCTION ORDER " << config->order << std::endl;
        }
		auto TIME1 = omp_get_wtime();
		timings[RECONSTRUCTION].push_back(TIME1 - TIME0);
		std::cout << "RECONSTRUCTION      time " << TIME1 - TIME0 << " s" << std::endl;

		// boundary condition
		// loop over all boundary faces
		#pragma omp parallel for schedule(dynamic)
		for (int ibf = 0; ibf < bcList.size(); ++ibf) {
			int jf = bcList[ibf]->jf;
			REAL x = mesh->faceCenters[jf][0];
			REAL y = mesh->faceCenters[jf][1];
			REAL z = mesh->faceCenters[jf][2];
			fLeftRight[jf][mesh->getOutIndex(jf)] = bcList[ibf]->applyBC(
					x, y, z,
					fLeftRight[jf][1 - mesh->getOutIndex(jf)],
					mesh->getOutSign(jf) * vn[jf],
					vn_abs[jf],
					config->tol
			);			
		}
		auto TIME2 = omp_get_wtime();
		timings[BOUNDARY_CONDITIONS].push_back(TIME2 - TIME1);
		std::cout << "BOUNDARY_CONDITIONS time " << TIME2 - TIME1 << " s" << std::endl;

		// Compute ranges for each thread
		//std::vector<int> threadRanges = getParallelRanges(timesForFaceFluxes, numThreads);

		// Riemann solver - compute fluxes
		// loop over all faces
		#pragma omp parallel for schedule(dynamic)
		for (int jf = 0; jf < mesh->nFaces; ++jf) {
			auto begin = omp_get_wtime();
			if (!config->isRusanov) {
			    // flux[jf] = 0.5 * mesh->faceAreas[jf] *
				// 	    (round_t(round_t(fLeftRight[jf][0] + fLeftRight[jf][1], config->tol) * vn[jf], config->tol) - 
				// 	        round_t(round_t(fLeftRight[jf][1] - fLeftRight[jf][0], config->tol) * vn_abs[jf], config->tol));
				flux[jf] = 0.5 * mesh->faceAreas[jf] * ((fLeftRight[jf][0] + fLeftRight[jf][1]) * vn[jf] - (fLeftRight[jf][1] - fLeftRight[jf][0]) * vn_abs[jf]);
			}
			else {
			    // flux[jf] = 0.5 * mesh->faceAreas[jf] *
				// 	    (round_t(round_t(fLeftRight[jf][0] + fLeftRight[jf][1], config->tol) * vn[jf], config->tol) - 
				// 	        vn_abs_max[jf] * round_t(fLeftRight[jf][1] - fLeftRight[jf][0], config->tol));
			    flux[jf] = 0.5 * mesh->faceAreas[jf] * ((fLeftRight[jf][0] + fLeftRight[jf][1]) * vn[jf] - vn_abs_max[jf] * (fLeftRight[jf][1] - fLeftRight[jf][0]));
			}
			flux[jf].round(config->tol);
			auto end = omp_get_wtime();
			timesForFaceFluxes[jf] = end - begin;
		}
		auto TIME3 = omp_get_wtime();
		timings[FLUXES].push_back(TIME3 - TIME2);
		std::cout << "FLUXES              time " << TIME3 - TIME2 << " s" << std::endl;

		// Compute ranges for each thread
		//threadRanges = getParallelRanges(timesForCellsRHS, numThreads);

		// computation of the right hand side
		// loop over all cells
		#pragma omp parallel for schedule(dynamic)
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			auto begin = omp_get_wtime();
			// compute macroparameters and collision integral
			std::vector<REAL> params = comp_macro_params(f[ic], v, gas_params);
			rhs[ic] = comp_j(params, f[ic], config->tol, v, gas_params);
			// sum up fluxes from all faces of this cell
			for (int j = 0; j < mesh->cellFaces[ic].size(); ++j) {
				int jf = mesh->cellFaces[ic][j];
				rhs[ic] = rhs[ic] - (mesh->getOutSign(ic, j) / mesh->cellVolumes[ic]) * flux[jf];
				rhs[ic].round(config->tol);
			}
			auto end = omp_get_wtime();
			timesForCellsRHS[ic] = end - begin;

			n[ic]   = params[0];
			ux[ic]  = params[1];
			uy[ic]  = params[2];
			uz[ic]  = params[3];
			T[ic]   = params[4];
			nu[ic]  = params[5];
			compression[ic] = f[ic].compression();
			rank_x[ic] = f[ic].r()[0];
			rank_y[ic] = f[ic].r()[1];
			rank_z[ic] = f[ic].r()[2];
			max_rank[ic] = std::max({f[ic].r()[0], f[ic].r()[1], f[ic].r()[2]});
			
			REAL u = pow(ux[ic]*ux[ic] + uy[ic]*uy[ic] + uz[ic]*uz[ic], 0.5)*gas_params->v_s;
			REAL Mach = u / pow(gas_params->g * gas_params->Rg * T[ic]*gas_params->T_s, 0.5);

			data[ic] = {
					n[ic]       * gas_params->n_s,   // n
					ux[ic]      * gas_params->v_s,   // ux
					uy[ic]      * gas_params->v_s,   // uy
					uz[ic]      * gas_params->v_s,   // uz
					T[ic]       * gas_params->T_s,   // T
					n[ic]       * gas_params->rho_s, // rho
					n[ic]*T[ic] * gas_params->p_s,   // p
					Mach,                            // Mach
					compression[ic],
					rank_x[ic],
					rank_y[ic],
					rank_z[ic],
					max_rank[ic]
			};
		}

		REAL frob_norm = vector_norm(rhs);
		update_res(frob_norm);
		auto TIME4 = omp_get_wtime();
		timings[RHS].push_back(TIME4 - TIME3);
		std::cout << "RHS                 time " << TIME4 - TIME3 << " s" << std::endl;

		// Update values
		// 
		if (!config->isImplicit) {
			#pragma omp parallel for schedule(dynamic)
			for (int ic = 0; ic < mesh->nCells; ++ic) {
				auto begin = omp_get_wtime();
				f[ic] = f[ic] + tau * rhs[ic];
				f[ic].round(config->tol);
				auto end = omp_get_wtime();
				timesForCellsUpdate[ic] = end - begin;
			}
		}
		else {

			auto divide_start = omp_get_wtime();
			mesh->divideMesh(numThreads, lusgs_timings);

			// BEGIN PARTITION PLOT
			std::vector < std::vector <REAL> > data;
			for (int ic = 0; ic < mesh->nCells; ++ic) {
				data.push_back(std::vector <REAL> {static_cast<REAL>(mesh->cellPartitions[ic]), static_cast<REAL>(mesh->cellColors[ic])});
			}
            std::ostringstream it_ss;
            it_ss << it;
            std::string it_string = it_ss.str();
			mesh->write_tecplot(data, "partiton_" + it_string + ".dat", {"partition", "color"});
			// END PARTITION PLOT

			std::fill_n(lusgs_timings.begin(), mesh->nCells, 0.0);
			auto divide_end   = omp_get_wtime();
			std::cout << "Divide time " << divide_end - divide_start << " s." << std::endl;
			
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
					auto start = omp_get_wtime();

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
						    if (!config->isRusanov) {
								// TODO important
							    vnm_loc = 0.5 * (-v->vn_abs_r1 + mesh->getOutSign(ic, j) * vn[jf]); // vnm[jf] or -vnp[jf]
							    // vnm_loc = 0.5 * (-vn_abs[jf] + mesh->getOutSign(ic, j) * vn[jf]); // vnm[jf] or -vnp[jf]
                            }
                            else {
                                vnm_loc = 0.5 * (-vn_abs_max[jf] * v->ones + mesh->getOutSign(ic, j) * vn[jf]); // vnm[jf] or -vnp[jf]
                            }
						    // df[ic] = df[ic] - (mesh->faceAreas[jf] / mesh->cellVolumes[ic]) * round_t(vnm_loc * df[icn], config->tol);
							df[ic] = df[ic] - (mesh->faceAreas[jf] / mesh->cellVolumes[ic]) * (vnm_loc * df[icn]);
						    df[ic].round(config->tol);
						}
					}
					// divide by diagonal coefficient
					div_tmp = ((1.0 / tau + nu[ic]) * v->ones + diag_r1[ic]);
					div_tmp.round(static_cast<REAL>(1e-3), 1); // TODO magic number
					// df[ic].round(config->tol); // TODO less rounding?
					df[ic] = df[ic] / div_tmp;

					auto finish = omp_get_wtime();
					lusgs_timings[ic] += finish - start;
				}
				#pragma omp barrier
			}

			// Forward sweep
			for (int color = 0; color < mesh->nColors; ++color) {
				for (int i = 0; i < mesh->C[partition][color].size(); ++i) {
					auto start = omp_get_wtime();

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
						    if (!config->isRusanov) {
								// TODO important
							    vnm_loc = 0.5 * (-v->vn_abs_r1 + mesh->getOutSign(ic, j) * vn[jf]); // vnm[jf] or -vnp[jf]
							    // vnm_loc = 0.5 * (-vn_abs[jf] + mesh->getOutSign(ic, j) * vn[jf]); // vnm[jf] or -vnp[jf]
					        }
					        else {
					            vnm_loc = 0.5 * (-vn_abs_max[jf] * v->ones + mesh->getOutSign(ic, j) * vn[jf]); // vnm[jf] or -vnp[jf]
					        }
						    // incr = incr - (mesh->faceAreas[jf] / mesh->cellVolumes[ic]) * round_t(vnm_loc * df[icn], config->tol);
						    incr = incr - (mesh->faceAreas[jf] / mesh->cellVolumes[ic]) * (vnm_loc * df[icn]);
						    incr.round(config->tol);
						}
					}
					// divide by diagonal coefficient
					div_tmp = ((1.0 / tau + nu[ic]) * v->ones + diag_r1[ic]);
					div_tmp.round(static_cast<REAL>(1e-3), 1); // TODO magic number
					df[ic] = df[ic] + (incr / div_tmp);
					df[ic].round(config->tol);

					auto finish = omp_get_wtime();
					lusgs_timings[ic] += finish - start;
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
			int  f_max_rank = 0;
			int df_max_rank = 0;
			for (int ic = 0; ic < mesh->nCells; ++ic) {
				int  f_max_rank_tmp = std::max({ f[ic].r()[0],  f[ic].r()[1],   f[ic].r()[2]});
				int df_max_rank_tmp = std::max({df[ic].r()[0], df[ic].r()[1],  df[ic].r()[2]});
				if (f_max_rank_tmp > f_max_rank) {
					f_max_rank = f_max_rank_tmp;
				}
				if (df_max_rank_tmp > df_max_rank) {
					df_max_rank = df_max_rank_tmp;
				}
			}
			std::cout << " f max rank = " <<  f_max_rank << std::endl;
			std::cout << "df max rank = " << df_max_rank << std::endl;
			std::cout << "df norm     = " << vector_norm(df) << std::endl;
//			for (int ic = 0; ic < mesh->nCells; ++ic) {
//				std::cout << "Tensor, r=(" << df[ic].r()[0] << "," << df[ic].r()[1] << "," << df[ic].r()[2] << ")";
//				std::cout << ", n=(" << df[ic].n()[0] << "," << df[ic].n()[1] << "," << df[ic].n()[2] << ")" << std::endl;
//			}
		}
		auto TIME5 = omp_get_wtime();
		timings[UPDATE].push_back(TIME5 - TIME4);
		std::cout << "UPDATE              time " << TIME5 - TIME4 << " s" << std::endl;

		// TODO save timings

		if ((it > 0) && (it % config->saveTecStep == 0)) {
            std::ostringstream it_ss;
            it_ss << it;
            std::string it_string = it_ss.str();
			mesh->write_tecplot(data, "tec_" + it_string + ".dat",
					{"n", "ux", "uy", "uz", "T", "rho", "p", "Mach", "compression", "rank_x", "rank_y", "rank_z", "max_rank"});
		}
		if ((it > 0) && (it % config->saveRestartStep == 0)) {
			write_restart();
		}
		if ((it > 0) && (it % config->saveMacroStep == 0)) {
			write_macro_restart();
			write_boundary_params();
		}
	}
	mesh->write_tecplot(data, "tec_final.dat",
			{"n", "ux", "uy", "uz", "T", "rho", "p", "Mach", "compression", "rank_x", "rank_y", "rank_z", "max_rank"});
	write_restart();
	write_macro_restart();
	write_boundary_params();
	// write_rate();
}

template class VelocityGrid<Full>;
template class Problem<Full>;
template class Solution<Full>;
template Full f_maxwell_t(std::shared_ptr < VelocityGrid<Full> >, REAL, REAL, REAL, REAL, REAL, REAL);
// template std::vector <REAL> comp_macro_params(const Full&, std::shared_ptr < VelocityGrid<Full> >, std::shared_ptr < GasParams >);     
// template Full comp_j(const std::vector <REAL>&, const Full&, REAL, std::shared_ptr < VelocityGrid<Full> > v, std::shared_ptr < GasParams >);
// template void Solution<Full>::reconstruction_2nd_order();


template class VelocityGrid<Tucker>;
template class Problem<Tucker>;
template class Solution<Tucker>;
template Tucker f_maxwell_t(std::shared_ptr < VelocityGrid<Tucker> >, REAL, REAL, REAL, REAL, REAL, REAL);
// template std::vector <REAL> comp_macro_params(const Full&, std::shared_ptr < VelocityGrid<Full> >, std::shared_ptr < GasParams >);     
// template Tucker comp_j(const std::vector <REAL>&, const Tucker&, REAL, std::shared_ptr < VelocityGrid<Tucker> > v, std::shared_ptr < GasParams >);
// template void Solution<Tucker>::reconstruction_2nd_order();


