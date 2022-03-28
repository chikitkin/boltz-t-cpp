#include "mesh.h"
#include "header.h"

#include "full.h"
#include "tucker.h"
#include "solver.h"

#include <ctime>

int main(int argc, char *argv[])
{
	std::shared_ptr < GasParams > gas_params = std::make_shared < GasParams > ();
    
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
	
	std::shared_ptr < Config > config = std::make_shared < Config > ();
	
	int steps;
		
	std::string mesh_path;
	
	
    std::string cfg_path = argv[1];
	std::ifstream cfg(cfg_path);
    std::string line;
    
    bool isTensorized;
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
        else if (line.find("isImplicit") != -1) sin >> config->isImplicit;
        else if (line.find("CFL") != -1) sin >> config->CFL;
            
        else if (line.find("steps") != -1) sin >> steps;
        else if (line.find("isTensorized") != -1) sin >> isTensorized;
    }
    //if (isTensorized) { typedef Tucker Tensor; }
    //else { typedef Full Tensor; }
    typedef Tucker Tensor;

	delta = 8.0 / (5.0 * pow(PI, 0.5) * Kn);
	u_l = Mach * pow(gas_params->g * gas_params->Rg * T_l, 0.5);

	double n_s = n_l;
	double T_s = T_l;

	double p_s = gas_params->m * n_s * gas_params->Rg * T_s;

	double v_s = pow(2. * gas_params->Rg * T_s, 0.5);
	double mu_s = gas_params->mu(T_s);

	double l_s = delta * mu_s * v_s / p_s;

	std::shared_ptr < Mesh > mesh = std::make_shared < Mesh > (mesh_path, l_s);

	if (mesh_path == "../mesh-1d/" || mesh_path == "../mesh-1d-tetra/") {
		n_r = (gas_params->g + 1.0) * Mach * Mach /
			((gas_params->g - 1.0) * Mach * Mach + 2.0) * n_l;
	    u_r = ((gas_params->g - 1.0) * Mach * Mach + 2.0) /
			((gas_params->g + 1.0) * Mach * Mach) * u_l;
	    T_r = (2.0 * gas_params->g * Mach * Mach - (gas_params->g - 1.0)) * ((gas_params->g - 1.0) * Mach * Mach + 2.0) /
			(pow(gas_params->g + 1.0, 2.0) * Mach * Mach) * T_l;
	}
	else if (cfg_path == "../mesh-cyl/") {
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

	std::shared_ptr < VelocityGrid<Tensor> > v = std::make_shared < VelocityGrid<Tensor> > (nv, nv, nv, vx_, vx_, vx_);
	
	Tensor f_in = f_maxwell_t<Tensor>(v, n_l, u_l, 0.0, 0.0, T_l, gas_params->Rg);
	Tensor f_out = f_maxwell_t<Tensor>(v, n_r, u_r, 0.0, 0.0, T_r, gas_params->Rg);

	std::vector < double > params;

	params = comp_macro_params<Tensor>(f_in, v, gas_params);

	std::cout << "n " << (params[0] - n_l) / n_l << " = 0" << std::endl;
	std::cout << "ux " << (params[1] - u_l) / u_l << " = 0" << std::endl;
	std::cout << "uy " << (params[2]) << " = 0" << std::endl;
	std::cout << "uz " << (params[3]) << " = 0" << std::endl;
	std::cout << "T " << (params[4] - T_l) / T_l << " = 0" << std::endl;


	std::shared_ptr < Problem<Tensor> > problem = std::make_shared < Problem<Tensor> > ();
	problem->initData = {f_in, f_out};

	Tensor fmaxwell = f_maxwell_t<Tensor>(v, 1.0, 0.0, 0.0, 0.0, T_w, gas_params->Rg);

	//                 SYMZ      INL   OUTL   WALL  SYMY      SYMX
//	problem->bcData = {Tensor(), f_in, f_out, fmaxwell, Tensor()};
//	problem->bcTypes = {SYMMETRYZ, INLET, OUTLET, WALL, SYMMETRYY};
	problem->bcData = {Tensor(), f_in, f_out, Tensor()}; // TODO fix
	problem->bcTypes = {SYMMETRYZ, INLET, OUTLET, SYMMETRYY};


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

	std::ofstream out;
	out.open("T.txt");
	for (int ic = 0; ic < S.mesh->nCells; ++ic) {
		out << S.mesh->cellCenters[ic][0] << " " << S.n[ic] << " " << S.ux[ic] << " " << S.T[ic] << "\n";
	}
	out.close();

	return 0;
}
