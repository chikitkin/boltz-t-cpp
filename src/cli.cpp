#include "mesh.h"
#include "header.h"

#include "full.h"
#include "tucker.h"
#include "solver.h"

#include <ctime>

int main(int argc, char *argv[])
{
    typedef Tucker Tensor;
	std::shared_ptr < GasParams > gas_params = std::make_shared < GasParams > ();

	std::shared_ptr < Problem<Tensor> > problem = std::make_shared < Problem<Tensor> > ();
	std::shared_ptr < Config > config = std::make_shared < Config > ();
    
    REAL Mach;
	REAL Kn;
	REAL delta;
	REAL l_s;
	
	REAL n_in;
	REAL u_in;
	REAL T_in;
	
	REAL n_out;
	REAL u_out;
	REAL T_out;

	int nv;
	
	int steps;
		
	std::string mesh_path;
	
	
    std::string cfg_path = argv[1];
	std::ifstream cfg(cfg_path);
    std::string line;
    
    std::string init;
    while (getline(cfg, line)) {
        std::istringstream line_stream(line.substr(line.find("=") + 1));
            
        if (line.find("mesh_path") != -1) { line_stream >> mesh_path; }
        else if (line.find("Mach") != -1) { line_stream >> Mach; }
        else if (line.find("Kn") != -1) { line_stream >> Kn; }
        else if (line.find("delta") != -1) { line_stream >> delta; }
        else if (line.find("l_s") != -1) { line_stream >> l_s; }
        else if (line.find("n_in") != -1) { line_stream >> n_in; }
        else if (line.find("u_in") != -1) { line_stream >> u_in; }
        else if (line.find("T_in") != -1) { line_stream >> T_in; }
        else if (line.find("n_out") != -1) { line_stream >> n_out; }
        else if (line.find("u_out") != -1) { line_stream >> u_out; }
        else if (line.find("T_out") != -1) { line_stream >> T_out; }
        
        else if (line.find("nv") != -1) { line_stream >> nv; }
        else if (line.find("CFL") != -1) { line_stream >> config->CFL; }
        else if (line.find("isImplicit") != -1) { line_stream >> config->isImplicit; }
        else if (line.find("tol") != -1) { line_stream >> config->tol; }
        else if (line.find("steps") != -1) { line_stream >> steps; }

		else if (line.find("saveTecStep") != -1) { line_stream >> config->saveTecStep; }
		else if (line.find("saveMacroStep") != -1) { line_stream >> config->saveMacroStep; }

        else if (line.find("boundary:") != -1) { break; }
	}
	
	mesh_path = argv[2];

//	delta = 8.0 / (5.0 * pow(PI, 0.5) * Kn);

	REAL n_s = n_in;
	REAL T_s = T_in;

	// REAL p_s = gas_params->m * n_s * gas_params->Rg * T_s;

	REAL v_s = pow(2. * gas_params->Rg * T_s, 0.5);
	// REAL mu_s = gas_params->mu(T_s);

//	l_s = delta * mu_s * v_s / p_s;

	std::shared_ptr < Mesh > mesh = std::make_shared < Mesh > (mesh_path, l_s);

	REAL vmax = 22.0 * v_s;
	REAL hv = 2.0 * vmax / nv;
	REAL *vx_ = new REAL[nv];
	for (int i = 0; i < nv; ++i) {
		vx_[i] = - vmax + (hv / 2.0) + i * hv;
	}

	std::shared_ptr < VelocityGrid<Tensor> > v = std::make_shared < VelocityGrid<Tensor> > (nv, nv, nv, vx_, vx_, vx_);
	
	Tensor f_in = f_maxwell_t<Tensor>(v, n_in, u_in, 0.0, 0.0, T_in, gas_params->Rg);
	Tensor f_out = f_maxwell_t<Tensor>(v, n_out, u_out, 0.0, 0.0, T_out, gas_params->Rg);

	problem->gas_params = gas_params;
	problem->v = v;
	problem->initData = {f_in, f_out};

	{
		int tag;
		REAL n, ux, uy, uz, T;
		REAL T_wall;
		std::cout << "BC types are:" << std::endl;
		while (getline(cfg, line)) {
			std::string bc_type;
			std::istringstream bc_line_stream(line);
			std::istringstream bc_stream(line.substr(line.find(" ") + 1));
			bc_stream >> bc_type;
			std::cout << bc_type << std::endl;
			if (bc_type == "WALL") {
				bc_line_stream >> tag >> bc_type >> T_wall;
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
				problem->bcData.push_back(f_maxwell_t<Tensor>(v, n, ux, uy, uz, T, gas_params->Rg));
			}
			else if (bc_type == "OUTLET") {
				bc_line_stream >> tag >> bc_type >> n >> ux >> uy >> uz >> T;
				problem->bcTags.push_back(tag);
				problem->bcTypes.push_back(OUTLET);
				problem->bcData.push_back(f_maxwell_t<Tensor>(v, n, ux, uy, uz, T, gas_params->Rg));
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

//	std::vector < REAL > params;
//
//	params = comp_macro_params<Tensor>(f_in, v, gas_params);
//
//	std::cout << "n " << (params[0] - n_l) / n_l << " = 0" << std::endl;
//	std::cout << "ux " << (params[1] - u_l) / u_l << " = 0" << std::endl;
//	std::cout << "uy " << (params[2]) << " = 0" << std::endl;
//	std::cout << "uz " << (params[3]) << " = 0" << std::endl;
//	std::cout << "T " << (params[4] - T_l) / T_l << " = 0" << std::endl;


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
