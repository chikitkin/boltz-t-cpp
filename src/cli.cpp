#include "read_mesh.h"
#include "header.h"

#include "full.h"
#include "tucker.h"
#include "solver.h"

#include <ctime>

int main(int argc, char *argv[])
{
	typedef Full Tensor;

	GasParams gas_params;
    
    double Mach;
	double Kn;
	double delta;
	
	double n_l;
	double u_l;
	double T_l;
	double T_w;
	
	double n_r;
	double u_r;
	double T_r;

	int nv;
	
	Config config;
	
	int steps;
		
	std::string mesh_path;
	
	
    std::string cfg_path = argv[1];
	std::ifstream cfg(cfg_path);
    std::string line;
    
    while (getline(cfg, line)) {
        std::istringstream sin(line.substr(line.find("=") + 1));
        if (line.find("mesh_path") != -1) sin >> mesh_path;
            
        else if (line.find("Mach") != -1) sin >> Mach;
        else if (line.find("Kn") != -1) sin >> Kn;
        else if (line.find("delta") != -1) sin >> delta;
        else if (line.find("n_l") != -1) sin >> n_l;
        else if (line.find("u_l") != -1) sin >> u_l;
        else if (line.find("T_l") != -1) sin >> T_l;
        else if (line.find("T_w") != -1) sin >> T_w;
        
        else if (line.find("nv") != -1) sin >> nv;
        else if (line.find("solver_type") != -1) sin >> config.solver_type;
        else if (line.find("CFL") != -1) sin >> config.CFL;
            
        else if (line.find("steps") != -1) sin >> steps;
    }

	delta = 8.0 / (5.0 * pow(PI, 0.5) * Kn);
	u_l = Mach * pow(gas_params.g * gas_params.Rg * T_l, 0.5);

	double n_s = n_l;
	double T_s = T_l;

	double p_s = gas_params.m * n_s * gas_params.Rg * T_s;

	double v_s = pow(2. * gas_params.Rg * T_s, 0.5);
	double mu_s = gas_params.mu(T_s);

	double l_s = delta * mu_s * v_s / p_s;

	Mesh mesh(mesh_path, l_s);
	if (cfg_path == "/home/egor/git/boltz-t-cpp/mesh-1d/1d.cfg") {
		n_r = (gas_params.g + 1.0) * Mach * Mach /
			((gas_params.g - 1.0) * Mach * Mach + 2.0) * n_l;
	    u_r = ((gas_params.g - 1.0) * Mach * Mach + 2.0) /
			((gas_params.g + 1.0) * Mach * Mach) * u_l;
	    T_r = (2.0 * gas_params.g * Mach * Mach - (gas_params.g - 1.0)) * ((gas_params.g - 1.0) * Mach * Mach + 2.0) /
			(pow(gas_params.g + 1.0, 2.0) * Mach * Mach) * T_l;
	}
	else if (cfg_path == "/home/egor/git/boltz-t-cpp/mesh-cyl/cyl.cfg") {
	    n_r = n_l;
	    u_r = u_l;
	    T_r = T_l;
	}

	double vmax = 22.0 * v_s;

	double hv = 2.0 * vmax / nv;

	double *vx_ = new double[nv];

	for (int i = 0; i < nv; ++i) {
		vx_[i] = - vmax + (hv / 2.0) + i * hv;
	}

	VelocityGrid<Tensor> v(nv, nv, nv, vx_, vx_, vx_);
	
	Tensor f_in = f_maxwell_t<Tensor>(v, n_l, u_l, 0.0, 0.0, T_l, gas_params.Rg);
	Tensor f_out = f_maxwell_t<Tensor>(v, n_r, u_r, 0.0, 0.0, T_r, gas_params.Rg);

	std::vector < double > params;

	params = comp_macro_params<Tensor>(f_in, v, gas_params);

	std::cout << "n " << (params[0] - n_l) / n_l << " = 0" << std::endl;
	std::cout << "ux " << (params[1] - u_l) / u_l << " = 0" << std::endl;
	std::cout << "uy " << (params[2]) << " = 0" << std::endl;
	std::cout << "uz " << (params[3]) << " = 0" << std::endl;
	std::cout << "T " << (params[4] - T_l) / T_l << " = 0" << std::endl;


	Problem<Tensor> problem;
	problem.init_tensor_list = {f_in, f_out};

	Tensor fmax = f_maxwell_t<Tensor>(v, 1.0, 0.0, 0.0, 0.0, T_w, gas_params.Rg);

	//                 SYMZ      INL   OUTL   WALL  SYMY      SYMX
	problem.bc_data = {Tensor(), f_in, f_out, fmax, Tensor()};
	problem.bc_types = {SYMMETRYZ, INLET, OUTLET, WALL, SYMMETRYY};

	Solution<Tensor> S(gas_params, mesh, v, problem, config);

	std::time_t t0 = std::time(nullptr);
	S.make_time_steps(config, steps);
	std::time_t t1 = std::time(nullptr);

	std::cout << "Time: " << t1 - t0 << " seconds." << std::endl;

	std::cout << "n = " << S.n[39] << std::endl;
	std::cout << "ux = " << S.ux[39] << std::endl;
	std::cout << "uy = " << S.uy[39] << std::endl;
	std::cout << "uz = " << S.uz[39] << std::endl;
	std::cout << "T = " << S.T[39] << std::endl;

/*
	int jf = 200;

	std::cout << "S.flux[jf]" << std::endl;
	std::cout << S.flux[jf] << std::endl;

	double *fluxjf = S.flux[jf].full();

	double s = 0.0;
	for (int i = 0; i < pow(44, 3); ++i) {
		s += fluxjf[i];
	}
	std::cout << s << std::endl;

	std::cout << "S.fm[jf]" << std::endl;
	std::cout << S.fm[jf] << std::endl;
	std::cout << "S.fp[jf]" << std::endl;
	std::cout << S.fp[jf] << std::endl;
	std::cout << "S.vn[jf]" << std::endl;
	std::cout << S.vn[jf] << std::endl;
	std::cout << "S.vn_abs[jf]" << std::endl;
	std::cout << S.vn_abs[jf] << std::endl;

	double *Sf = S.f[40].full();
	double *Sfm = S.fm[jf].full();
	double *Sfp = S.fp[jf].full();
	double *Svn = S.vn[jf].full();
	double *Svn_abs = S.vn_abs[jf].full();
	double *Sflux = S.flux[jf].full();

	std::cout << "Sf[i]" << std::endl;
	for (int i = 0; i < 10; ++i) {
		std::cout << Sf[i] << std::endl;
	}
	std::cout << "Sfm[i]" << std::endl;
	for (int i = 0; i < 10; ++i) {
		std::cout << Sfm[i] << std::endl;
	}
	std::cout << "Sfp[i]" << std::endl;
	for (int i = 0; i < 10; ++i) {
		std::cout << Sfp[i] << std::endl;
	}
	std::cout << "Svn[i]" << std::endl;
	for (int i = 0; i < 10; ++i) {
		std::cout << Svn[i] << std::endl;
	}
	std::cout << "Svn_abs[i]" << std::endl;
	for (int i = 0; i < 10; ++i) {
		std::cout << Svn_abs[i] << std::endl;
	}
	std::cout << "Sflux[i]" << std::endl;
	for (int i = 0; i < 10; ++i) {
		std::cout << Sflux[i] << std::endl;
	}

	std::ofstream file;
	file.open("/home/egor/git/boltz-t-cpp/run.txt");
	file << "Sflux[i]" << std::endl;
	for (int i = 0; i < S.v.nv; ++i) {
		file << Sflux[i] << std::endl;
	}

	std::cout << "S.vn_round[jf]" << std::endl;
	std::cout << round_t(S.vn[jf], 1e-3, 1e+6) << std::endl;

	file << "Svn[i]" << std::endl;
	for (int i = 0; i < S.v.nv; ++i) {
		file << Svn[i] << std::endl;
	}
/*
	for (int j = 0; j < S.mesh.nf; ++j) {
		std::cout << S.fm[j];
	}
*/
	return 0;
}
