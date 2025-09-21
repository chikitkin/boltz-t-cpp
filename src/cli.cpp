#include "mesh.h"
#include "header.h"

#include "full.h"
#include "tucker.h"
#include "solver.h"

#include <ctime>

template <class Tensor>
int check_velocity_grid(REAL n, REAL ux, REAL uy, REAL uz, REAL T, 
        std::shared_ptr < VelocityGrid<Tensor> > v,
        std::shared_ptr < GasParams > gas_params)
{
    Tensor f = f_maxwell_t<Tensor>(v, n, ux, uy, uz, T, gas_params->Rg);
    std::vector<REAL> params = comp_macro_params(f, v, gas_params);
    
    std::cout << "\tn:  " << abs(params[0] - n) / n << " = 0" << std::endl;
	std::cout << "\tux: " << abs(params[1] - ux) / ux << " = 0" << std::endl;
	std::cout << "\tuy: " << abs(params[2] - uy) / uy << " = 0" << std::endl;
	std::cout << "\tuz: " << abs(params[3] - uz) / uz << " = 0" << std::endl;
	std::cout << "\tT:  " << abs(params[4] - T) / T   << " = 0" << std::endl;
	
	return 0;
}

int main(int argc, char *argv[])
{
    typedef Tucker Tensor;
	std::shared_ptr < GasParams > gas_params = std::make_shared < GasParams > ();

	std::shared_ptr < Problem<Tensor> > problem = std::make_shared < Problem<Tensor> > ();
	std::shared_ptr < Config > config = std::make_shared < Config > ();
	
	REAL n_in;
	REAL ux_in = 0.0;
	REAL uy_in = 0.0;
	REAL uz_in = 0.0;
	REAL T_in;
	
	REAL n_out;
	REAL ux_out = 0.0;
	REAL uy_out = 0.0;
	REAL uz_out = 0.0;
	REAL T_out;

	REAL channel_length = -1.0; // IN Z DIRECTION

	REAL u_in;
	REAL u_out;

	REAL vmax = 22.0;
	int nvx;
	int nvy;
	int nvz;
	
	int steps;
		
	std::string mesh_path;
	
    std::string cfg_path = argv[1];
	std::ifstream cfg(cfg_path);
	cfg.precision(17); // TODO magic number
    std::string line;
    
    while (getline(cfg, line)) {
        std::istringstream line_stream(line.substr(line.find("=") + 1));
		line_stream.precision(17); // TODO magic number
            
        if      (line.find("l_s")   != std::string::npos) { line_stream >> gas_params->l_s; }
        else if (line.find("n_in")  != std::string::npos) { line_stream >> n_in; }
        else if (line.find("ux_in") != std::string::npos) { line_stream >> ux_in; }
        else if (line.find("uy_in") != std::string::npos) { line_stream >> uy_in; }
        else if (line.find("uz_in") != std::string::npos) { line_stream >> uz_in; }
        else if (line.find("T_in")  != std::string::npos) { line_stream >> T_in; }
        else if (line.find("n_out")  != std::string::npos) { line_stream >> n_out; }
        else if (line.find("ux_out") != std::string::npos) { line_stream >> ux_out; }
        else if (line.find("uy_out") != std::string::npos) { line_stream >> uy_out; }
        else if (line.find("uz_out") != std::string::npos) { line_stream >> uz_out; }
        else if (line.find("T_out")  != std::string::npos) { line_stream >> T_out; }
        
        else if (line.find("gas_params_Mol")   != std::string::npos) { line_stream >> gas_params->Mol; }
        else if (line.find("gas_params_Pr")    != std::string::npos) { line_stream >> gas_params->Pr; }
        else if (line.find("gas_params_C")     != std::string::npos) { line_stream >> gas_params->C; }
        else if (line.find("gas_params_T_0")   != std::string::npos) { line_stream >> gas_params->T_0; }
        else if (line.find("gas_params_mu_0")  != std::string::npos) { line_stream >> gas_params->mu_0; }
        else if (line.find("gas_params_omega") != std::string::npos) { line_stream >> gas_params->omega; }
        else if (line.find("gas_params_g")     != std::string::npos) { line_stream >> gas_params->g; }
        
        else if (line.find("nvx")   != std::string::npos) { line_stream >> nvx; }
        else if (line.find("nvy")   != std::string::npos) { line_stream >> nvy; }
        else if (line.find("nvz")   != std::string::npos) { line_stream >> nvz; }
        else if (line.find("vmax")  != std::string::npos) { line_stream >> vmax; }
        else if (line.find("CFL")   != std::string::npos) { line_stream >> config->CFL; }
        else if (line.find("tol")   != std::string::npos) { line_stream >> config->tol; }
        else if (line.find("order") != std::string::npos) { line_stream >> config->order; }
        else if (line.find("steps") != std::string::npos) { line_stream >> steps; }
        else if (line.find("isImplicit")  != std::string::npos) { line_stream >> config->isImplicit; }
        else if (line.find("isIncrement") != std::string::npos) { line_stream >> config->isIncrement; }

		else if (line.find("channel_length") != std::string::npos) { line_stream >> channel_length; }
        
        else if (line.find("initType")        != std::string::npos) { line_stream >> config->initType; }
		else if (line.find("saveTecStep")     != std::string::npos) { line_stream >> config->saveTecStep; }
		else if (line.find("saveMacroStep")   != std::string::npos) { line_stream >> config->saveMacroStep; }
		else if (line.find("saveRestartStep") != std::string::npos) { line_stream >> config->saveRestartStep; }
		else if (line.find("vnAbsRestart")    != std::string::npos) { line_stream >> config->vnAbsRestart; }

        else if (line.find("boundary:") != -1) { break; }
	}

	u_in = pow(ux_in*ux_in + uy_in*uy_in + uz_in*uz_in, 0.5);
	u_out = pow(ux_out*ux_out + uy_out*uy_out + uz_out*uz_out, 0.5);
	
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

	gas_params->v_s = pow(2.0 * gas_params->Rg * gas_params->T_s, 0.5);
	gas_params->mu_s = gas_params->mu_suth(gas_params->T_s);
	
	REAL c = pow(gas_params->g * gas_params->Rg * gas_params->T_s, 0.5);
	gas_params->S_inf = u_in / gas_params->v_s;
	
	REAL lambda       = (gas_params->mu_s / gas_params->p_s) * pow((PI * gas_params->Rg * gas_params->T_s) / 2.0, 0.5);
	gas_params->delta = (gas_params->l_s * gas_params->p_s) / (gas_params->mu_s * gas_params->v_s);
	gas_params->Kn    = 8.0 / (5.0 * pow(PI, 0.5)) / gas_params->delta;
	REAL Mach         = u_in / c;
	REAL Re           = gas_params->rho_s * u_in * gas_params->l_s / gas_params->mu_s;
	
	// print parameters
	std::cout << "n^{star}   = " << gas_params->n_s   << std::endl;
	std::cout << "T^{star}   = " << gas_params->T_s   << std::endl;
	std::cout << "rho^{star} = " << gas_params->rho_s << std::endl;
	std::cout << "p^{star}   = " << gas_params->p_s   << std::endl;
	std::cout << "v^{star}   = " << gas_params->v_s   << std::endl;
	std::cout << "mu^{star}  = " << gas_params->mu_s  << std::endl;
	std::cout << "S^{inf}    = " << gas_params->S_inf << std::endl;
	
	std::cout << "Speed of sound,            c = " << c                  << std::endl;
	std::cout << "Rarefaction parameter, delta = " << gas_params->delta  << std::endl;
	std::cout << "Knudsen number,           Kn = " << gas_params->Kn     << std::endl;
	std::cout << "Mean free path,       lambda = " << lambda             << std::endl;
	std::cout << "Mach number,            Mach = " << Mach               << std::endl;
	std::cout << "Reynolds numbers,         Re = " << Re                 << std::endl;

	std::shared_ptr < Mesh > mesh = std::make_shared < Mesh > (mesh_path, 1.0); // gas_params->l_s);

	std::cout << "START DIMENSIONLESS" << std::endl;
	
	REAL hvx = 2.0 * vmax / (nvx - 1);
	REAL *vx_ = new REAL[nvx];
	for (int i = 0; i < nvx; ++i) {
		vx_[i] = - vmax + i * hvx;
		std::cout << vx_[i] << " ";
	}
	std::cout << "\n";
	
	REAL hvy = 2.0 * vmax / (nvy - 1);
	REAL *vy_ = new REAL[nvy];
	for (int i = 0; i < nvy; ++i) {
		vy_[i] = - vmax + i * hvy;
	}
	
	REAL hvz = 2.0 * vmax / (nvz - 1);
	REAL *vz_ = new REAL[nvz];
	for (int i = 0; i < nvz; ++i) {
		vz_[i] = - vmax + i * hvz;
	}
	
	std::cout << "v_min =  " << vx_[0]     << std::endl;
	std::cout << "v_max =  " << vx_[nvx-1] << std::endl;
	std::cout << "v step = " << hvx        << std::endl;
	
	std::shared_ptr < VelocityGrid<Tensor> > v = std::make_shared < VelocityGrid<Tensor> > (nvx, nvy, nvz, vx_, vy_, vz_);

	n_in  /= gas_params->n_s;
	ux_in /= gas_params->v_s;
	uy_in /= gas_params->v_s;
	uz_in /= gas_params->v_s;
	T_in  /= gas_params->T_s;

	n_out  /= gas_params->n_s;
	ux_out /= gas_params->v_s;
	uy_out /= gas_params->v_s;
	uz_out /= gas_params->v_s;
	T_out  /= gas_params->T_s;

	u_in  /= gas_params->v_s;
	u_out /= gas_params->v_s;
	
	Tensor f_in  = f_maxwell_t<Tensor>(v, n_in, ux_in, uy_in, uz_in, T_in, gas_params->Rg);
	Tensor f_out = f_maxwell_t<Tensor>(v, n_out, ux_out, uy_out, uz_out, T_out, gas_params->Rg);

	// std::cout << "f_in string: " << f_in.to_string() << std::endl;
	std::cout << "f_in string error norm: " << (from_string(f_in.to_string(), f_in) - f_in).norm() / f_in.norm() << std::endl;

	// std::cout << "f_out string: " << f_out.to_string() << std::endl;
	std::cout << "f_out string error norm: " << (from_string(f_out.to_string(), f_in) - f_out).norm() / f_out.norm() << std::endl;

	std::cout << "Test minmod" << std::endl;
	std::cout << "f_in:  " << f_in  << std::endl;
	std::cout << "f_out: " << f_out << std::endl;
	std::cout << "f_in norm:  " << f_in.norm()  << std::endl;
	std::cout << "f_out norm: " << f_out.norm() << std::endl;
	std::cout << minmod(f_in, f_out, config->tol) << std::endl;

	problem->gas_params = gas_params;
	problem->v = v;
	problem->initData = {f_in, f_out};
	
	problem->params_in  = {n_in,  ux_in,  uy_in,  uz_in,  T_in};
	problem->params_out = {n_out, ux_out, uy_out, uz_out, T_out};
	
	// Rankine-Hugoniot
	REAL n_rh = (gas_params->g + 1.) * Mach * Mach / ((gas_params->g - 1.) * Mach * Mach + 2.) * n_in;
	REAL u_rh = ((gas_params->g - 1.) * Mach * Mach + 2.) / ((gas_params->g + 1.) * Mach * Mach) * u_in;
	REAL T_rh = (2. * gas_params->g * Mach * Mach - (gas_params->g - 1.)) * ((gas_params->g - 1.) * Mach * Mach + 2.) / (pow(gas_params->g + 1, 2) * Mach * Mach) * T_in;
	std::cout << "Rankine-Hugoniot n, u, T" << std::endl;
	std::cout << n_rh*gas_params->n_s << " " << u_rh*gas_params->v_s << " " << T_rh*gas_params->T_s << std::endl;

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
	cfg.close();
    
    std::cout << "Check v ranges" << std::endl;
	std::cout << "Inlet:" << std::endl;
	check_velocity_grid(n_in, ux_in, uy_in, uz_in, T_in, v, gas_params);
    std::cout << "Outlet:" << std::endl;
    check_velocity_grid(n_out, ux_out, uy_out, uz_out, T_out, v, gas_params);
    std::cout << "Rankine-Hugoniot:" << std::endl;
    check_velocity_grid(n_rh, u_rh, 0.0, 0.0, T_rh, v, gas_params);
	std::cout << "Wall:" << std::endl;
    check_velocity_grid(n_in, 0.0, 0.0, 0.0, T_wall, v, gas_params);
	std::cout << "Inlet 2:" << std::endl;
	std::vector<REAL> params = comp_macro_params(f_in, v, gas_params);
	std::cout << "\t" << abs(params[0] - n_in) / n_in << " = 0, " << abs(params[4] - T_in) / T_in << " = 0"  << std::endl;

	// WRITE MACRO START
	std::cout << "channel_length = " << channel_length << std::endl;
	if (channel_length > 0.0) {
		std::ofstream file;
		file.open("../macro_start.txt", std::ofstream::trunc);
		file.precision(17); // TODO magic number
		for (int ic = 0; ic < mesh->nCells; ++ic) {
			file << mesh->cellCenters[ic][0] << " " << mesh->cellCenters[ic][1] << " " << mesh->cellCenters[ic][2] << " ";
			file << 1.0 - (mesh->cellCenters[ic][2] / channel_length) << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << T_wall << "\n";
		}
		file.close();
	}

	Solution<Tensor> S(gas_params, mesh, v, problem, config);

	auto start = omp_get_wtime();
	S.make_time_steps(config, steps);
	auto end = omp_get_wtime();

	std::cout << "Time: " << end - start << " seconds." << std::endl;

	std::ofstream out;
	out.open("T.txt");
	out.precision(17); // TODO magic number
	for (int ic = 0; ic < S.mesh->nCells; ++ic) {
		out << 
		S.mesh->cellCenters[ic][0] << " " << 
		S.mesh->cellCenters[ic][1] << " " <<
		S.mesh->cellCenters[ic][2] << " " <<
		S.n[ic]  << " " << 
		S.ux[ic] << " " << 
		S.uy[ic] << " " << 
		S.uz[ic] << " " << 
		S.T[ic]  << "\n";
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
