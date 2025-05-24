#include "mesh.h"
#include "header.h"

#include "full.h"
#include "tucker.h"
#include "solver.h"

#include <ctime>

template <class Tensor>
int check_velocity_grid(REAL n, REAL u, REAL T, 
        std::shared_ptr < VelocityGrid<Tensor> > v,
        std::shared_ptr < GasParams > gas_params)
{
    Tensor f = f_maxwell_t<Tensor>(v, n, u, 0.0, 0.0, T, gas_params->Rg);
    std::vector<REAL> params = comp_macro_params(f, v, gas_params);
    
    std::cout << "n:  " << (params[0] - n) / n << " = 0" << std::endl;
	std::cout << "ux: " << (params[1] - u) / u << " = 0" << std::endl;
	std::cout << "uy: " << params[2] << " = 0" << std::endl;
	std::cout << "uz: " << params[3] << " = 0" << std::endl;
	std::cout << "T:  " << (params[4] - T) / T << " = 0" << std::endl;
	
	return 0;
}

int main(int argc, char *argv[])
{
    typedef Tucker Tensor;
	std::shared_ptr < GasParams > gas_params = std::make_shared < GasParams > ();

	std::shared_ptr < Problem<Tensor> > problem = std::make_shared < Problem<Tensor> > ();
	std::shared_ptr < Config > config = std::make_shared < Config > ();
    
    REAL Mach;
	REAL Kn;
	REAL delta;
	REAL lambda;
	REAL Re;

	REAL l_s;
	
	REAL n_in;
	REAL u_in;
	REAL T_in;
	
	REAL n_out;
	REAL u_out;
	REAL T_out;

	int nvx;
	int nvy;
	int nvz;
	
	int steps;
		
	std::string mesh_path;
	
	
    std::string cfg_path = argv[1];
	std::ifstream cfg(cfg_path);
    std::string line;
    
    while (getline(cfg, line)) {
        std::istringstream line_stream(line.substr(line.find("=") + 1));
            
        if (line.find("Mach") != -1) { line_stream >> Mach; }
        else if (line.find("l_s") != -1) { line_stream >> gas_params->l_s; }
        else if (line.find("n_in") != -1) { line_stream >> n_in; }
        else if (line.find("u_in") != -1) { line_stream >> u_in; }
        else if (line.find("T_in") != -1) { line_stream >> T_in; }
        else if (line.find("n_out") != -1) { line_stream >> n_out; }
        else if (line.find("u_out") != -1) { line_stream >> u_out; }
        else if (line.find("T_out") != -1) { line_stream >> T_out; }
        
        else if (line.find("Mol") != -1) { line_stream >> gas_params->Mol; }
        else if (line.find("Pr") != -1) { line_stream >> gas_params->Pr; }
        else if (line.find("C") != -1) { line_stream >> gas_params->C; }
        else if (line.find("T_0") != -1) { line_stream >> gas_params->T_0; }
        else if (line.find("mu_0") != -1) { line_stream >> gas_params->mu_0; }
        else if (line.find("omega") != -1) { line_stream >> gas_params->omega; }
        else if (line.find("g") != -1) { line_stream >> gas_params->g; }
        
        else if (line.find("nvx") != -1) { line_stream >> nvx; }
        else if (line.find("nvy") != -1) { line_stream >> nvy; }
        else if (line.find("nvz") != -1) { line_stream >> nvz; }
        else if (line.find("CFL") != -1) { line_stream >> config->CFL; }
        else if (line.find("isImplicit") != -1) { line_stream >> config->isImplicit; }
        else if (line.find("isRusanov") != -1) { line_stream >> config->isRusanov; }
        else if (line.find("isImplicitIncrement") != -1) { line_stream >> config->isImplicitIncrement; }
        else if (line.find("tol") != -1) { line_stream >> config->tol; }
        else if (line.find("order") != -1) { line_stream >> config->order; }
        else if (line.find("steps") != -1) { line_stream >> steps; }
        
        else if (line.find("initType") != -1) { line_stream >> config->initType; }

		else if (line.find("saveTecStep") != -1) { line_stream >> config->saveTecStep; }
		else if (line.find("saveMacroStep") != -1) { line_stream >> config->saveMacroStep; }

        else if (line.find("boundary:") != -1) { break; }
	}
	
	gas_params->Rg = gas_params->Ru / gas_params->Mol; // = self.Ru / self.Mol  # J / (kg * K)
	gas_params->m = gas_params->Mol / gas_params->Na; // # kg
	
	mesh_path = argv[2];
	
	if (config->initType != 0) {
	    config->initFilename = argv[3];
	}

	gas_params->n_s = n_in;
	gas_params->T_s = T_in;

	gas_params->rho_s = gas_params->m * gas_params->n_s;
	gas_params->p_s = gas_params->rho_s * gas_params->Rg * gas_params->T_s;

	gas_params->v_s = pow(2. * gas_params->Rg * gas_params->T_s, 0.5);
	gas_params->mu_s = gas_params->mu_suth(gas_params->T_s);
	
	REAL c = pow(gas_params->g * gas_params->Rg * gas_params->T_s, 0.5);
	REAL S_inf = u_in / gas_params->v_s;
	
	delta = (gas_params->l_s * gas_params->p_s) / (gas_params->mu_s * gas_params->v_s);
	gas_params->Kn = 8.0 / (5.0 * pow(PI, 0.5)) / delta;
	Mach = u_in / c;
	Re = gas_params->rho_s * u_in * gas_params->l_s / gas_params->mu_s;
	
	// print parameters
	std::cout << "p^{star} = " << gas_params->p_s << std::endl;
	std::cout << "mu^{star} = " << gas_params->mu_s << std::endl;
	
	std::cout << "v^{star} = " << gas_params->v_s << std::endl;
	std::cout << "mu^{star} = " << gas_params->mu_s << std::endl;
	std::cout << "S^{inf} = " << S_inf << std::endl;
	
	std::cout << "Speed of sound,            c = " << c      << std::endl;
	std::cout << "Rarefaction parameter, delta = " << delta  << std::endl;
	std::cout << "Knudsen number,           Kn = " << gas_params->Kn << std::endl;
	std::cout << "Mean free path,       lambda = " << lambda << std::endl;
	std::cout << "Mach number,            Mach = " << Mach   << std::endl;
	std::cout << "Reynolds numbers,         Re = " << Re     << std::endl;

	std::shared_ptr < Mesh > mesh = std::make_shared < Mesh > (mesh_path, 1.0); // gas_params->l_s);

	std::cout << "START DIMENSIONLESS" << std::endl;

	REAL vmax = 20.0; // WAS 22.0
	
	REAL hvx = 2.0 * vmax / nvx;
	REAL *vx_ = new REAL[nvx];
	for (int i = 0; i < nvx; ++i) {
		vx_[i] = - vmax + (hvx / 2.0) + i * hvx;
	}
	
	REAL hvy = 2.0 * vmax / nvy;
	REAL *vy_ = new REAL[nvy];
	for (int i = 0; i < nvy; ++i) {
		vy_[i] = - vmax + (hvy / 2.0) + i * hvy;
	}
	
	REAL hvz = 2.0 * vmax / nvz;
	REAL *vz_ = new REAL[nvz];
	for (int i = 0; i < nvz; ++i) {
		vz_[i] = - vmax + (hvz / 2.0) + i * hvz;
	}
	
	std::cout << "v_min =  " << vx_[0]     << std::endl;
	std::cout << "v_max =  " << vx_[nvx-1] << std::endl;
	std::cout << "v step = " << hvx        << std::endl;
	
	std::shared_ptr < VelocityGrid<Tensor> > v = std::make_shared < VelocityGrid<Tensor> > (nvx, nvy, nvz, vx_, vy_, vz_);

	n_in /= gas_params->n_s;
	u_in /= gas_params->v_s;
	T_in /= gas_params->T_s;

	n_out /= gas_params->n_s;
	u_out /= gas_params->v_s;
	T_out /= gas_params->T_s;
	
	Tensor f_in  = f_maxwell_t<Tensor>(v, n_in,  u_in,  0.0, 0.0, T_in,  gas_params->Rg);
	Tensor f_out = f_maxwell_t<Tensor>(v, n_out, u_out, 0.0, 0.0, T_out, gas_params->Rg);

	problem->gas_params = gas_params;
	problem->v = v;
	problem->initData = {f_in, f_out};
	
	problem->params_in  = {n_in,  u_in,  0.0, 0.0, T_in};
	problem->params_out = {n_out, u_out, 0.0, 0.0, T_out};
	
    // Inlet
	// std::vector<REAL> params = comp_macro_params(problem->initData[0], v, gas_params, T_s);

    // std::cout << "Inlet" << std::endl;
	// std::cout << "n:  " << (params[0] - problem->params_in[0]) / (problem->params_in[0]) << " = 0" << std::endl;
	// std::cout << "ux: " << (params[1] - problem->params_in[1]) / (problem->params_in[1]) << " = 0" << std::endl;
	// std::cout << "uy: " << problem->params_in[2] << " = 0" << std::endl;
	// std::cout << "uz: " << problem->params_in[3] << " = 0" << std::endl;
	// std::cout << "T:  " << (params[4] - problem->params_in[4]) / (problem->params_in[4]) << " = 0" << std::endl;
	
	// Rankine-Hugoniot
	REAL n_rh = (gas_params->g + 1.) * Mach * Mach / ((gas_params->g - 1.) * Mach * Mach + 2.) * n_in;
	REAL u_rh = ((gas_params->g - 1.) * Mach * Mach + 2.) / ((gas_params->g + 1.) * Mach * Mach) * u_in;
	REAL T_rh = (2. * gas_params->g * Mach * Mach - (gas_params->g - 1.)) * ((gas_params->g - 1.) * Mach * Mach + 2.) / (pow(gas_params->g + 1, 2) * Mach * Mach) * T_in;
	std::cout << "Rankine-Hugoniot n, u, T" << std::endl;
	std::cout << n_rh << " " << u_rh << " " << T_rh << std::endl;

    REAL T_wall = 200.0/gas_params->T_s; // TODO magic number, should not fail if no wall
	{
		int tag;
		REAL n, ux, uy, uz, T;
		std::cout << "BC types are:" << std::endl;
		while (getline(cfg, line)) {
			std::string bc_type;
			std::istringstream bc_line_stream(line);
			std::istringstream bc_stream(line.substr(line.find(" ") + 1));
			bc_stream >> bc_type;
			std::cout << bc_type << std::endl;
			if (bc_type == "WALL") {
				bc_line_stream >> tag >> bc_type >> T_wall;
				T_wall /= gas_params->T_s;
				problem->bcTags.push_back(tag);
				problem->bcTypes.push_back(WALL);
				problem->bcData.push_back(f_maxwell_t<Tensor>(v, 1.0, 0.0, 0.0, 0.0, T_wall, gas_params->Rg));
			}
			else if (bc_type == "SYMMETRYX") {
				bc_line_stream >> tag >> bc_type;
				problem->bcTags.push_back(tag);
				problem->bcTypes.push_back(SYMMETRYX);
				problem->bcData.push_back(Tensor());
			}
			else if (bc_type == "SYMMETRYY") {
				bc_line_stream >> tag >> bc_type;
				problem->bcTags.push_back(tag);
				problem->bcTypes.push_back(SYMMETRYY);
				problem->bcData.push_back(Tensor());
			}
			else if (bc_type == "SYMMETRYZ") {
				bc_line_stream >> tag >> bc_type;
				problem->bcTags.push_back(tag);
				problem->bcTypes.push_back(SYMMETRYZ);
				problem->bcData.push_back(Tensor());
			}
			else if (bc_type == "INLET") {
				bc_line_stream >> tag >> bc_type >> n >> ux >> uy >> uz >> T;
				problem->bcTags.push_back(tag);
				problem->bcTypes.push_back(INLET);
				problem->bcData.push_back(
					f_maxwell_t<Tensor>(v, 
						n/gas_params->n_s,
						ux/gas_params->v_s, 
						uy/gas_params->v_s, 
						uz/gas_params->v_s,
						T/gas_params->T_s,
						gas_params->Rg
					)
				);
			}
			else if (bc_type == "OUTLET") {
				bc_line_stream >> tag >> bc_type >> n >> ux >> uy >> uz >> T;
				problem->bcTags.push_back(tag);
				problem->bcTypes.push_back(OUTLET);
				problem->bcData.push_back(
					f_maxwell_t<Tensor>(v, 
						n/gas_params->n_s,
						ux/gas_params->v_s, 
						uy/gas_params->v_s, 
						uz/gas_params->v_s,
						T/gas_params->T_s,
						gas_params->Rg
					)
				);
			}
			else if (bc_type == "") {
				break;
			}
			else {
				std::cout << "Wrong boundary condition in cfg" << std::endl;
				exit(-1);
			}
		}
    }
    
    std::cout << "Check v ranges" << std::endl;
	std::cout << "Inlet" << std::endl;
	check_velocity_grid(n_in, u_in, T_in, v, gas_params);
    std::cout << "Rankine-Hugoniot" << std::endl;
    check_velocity_grid(n_out, u_out, T_out, v, gas_params);
	std::cout << "Wall" << std::endl;
    check_velocity_grid(n_in, 0.0, T_wall, v, gas_params);

	Solution<Tensor> S(gas_params, mesh, v, problem, config);

	auto start = omp_get_wtime();
	S.make_time_steps(config, steps);
	auto end = omp_get_wtime();

	std::cout << "Time: " << end - start << " seconds." << std::endl;

	std::ofstream out;
	out.open("T.txt");
	for (int ic = 0; ic < S.mesh->nCells; ++ic) {
		out << S.mesh->cellCenters[ic][0] << " " << S.n[ic] << " " << S.ux[ic] << " " << S.T[ic] << "\n";
	}
	out.close();

	out.open("timings.txt");
	for (int i = 0; i < S.timings[RECONSTRUCTION].size(); ++i) {
		out << S.timings[RECONSTRUCTION][i] << " " <<
		S.timings[BOUNDARY_CONDITIONS][i] << " " << 
		S.timings[FLUXES][i] << " " << 
		S.timings[RHS][i] << " " << 
		S.timings[UPDATE][i] << "\n";
	}
	out.close();	

	return 0;
}
