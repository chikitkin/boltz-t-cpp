#include "read_mesh.h"
#include "header.h"

int main() {

	Mesh mesh("/home/egor/git/boltz-t-cpp/mesh-1d/");

	std::cout << mesh.nbf << std::endl;
	std::cout << mesh.nc << std::endl;
	std::cout << mesh.nv << std::endl;
	std::cout << mesh.nf << std::endl;

	std::cout << "bcface_vert_lists" << std::endl;
	for (auto bcface : mesh.bcface_vert_lists) {
		for (auto v : bcface) {
			std::cout << v << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "bcface_vert_set_lists" << std::endl;
	for (auto bcface : mesh.bcface_vert_set_lists) {
		for (auto v : bcface) {
			std::cout << v << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "bcface_bctype" << std::endl;
	for (auto type : mesh.bcface_bctype) {
		std::cout << type << " ";
	}
	std::cout << "\n";

	std::cout << "vert_list_for_cell" << std::endl;
	for (auto i : mesh.vert_list_for_cell) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "bf_for_each_bc" << std::endl;
	for (auto i : mesh.bf_for_each_bc) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "vert_coo" << std::endl;
	for (auto i : mesh.vert_coo) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cell_center_coo" << std::endl;
	for (auto i : mesh.cell_center_coo) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cell_volumes" << std::endl;
	for (auto a : mesh.cell_volumes) {
		std::cout << a << " ";
	}
	std::cout << "\n";

	std::cout << "cell_list_for_vertex" << std::endl;
	for (auto i : mesh.cell_list_for_vertex) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cell_num_for_vertex" << std::endl;
	for (auto a : mesh.cell_num_for_vertex) {
		std::cout << a << " ";
	}
	std::cout << "\n";

	std::cout << "cell_neighbors_list" << std::endl;
	for (auto i : mesh.cell_neighbors_list) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "face_vert_list" << std::endl;
	for (auto i : mesh.face_vert_list) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cell_face_list" << std::endl;
	for (auto i : mesh.cell_face_list) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "face_areas" << std::endl;
	for (auto a : mesh.face_areas) {
		std::cout << a << " ";
	}
	std::cout << "\n";

	std::cout << "face_normals" << std::endl;
	for (auto i : mesh.face_normals) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cell_face_normal_direction" << std::endl;
	for (auto i : mesh.cell_face_normal_direction) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "bound_face_info" << std::endl;
	for (auto i : mesh.bound_face_info) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cell_diam" << std::endl;
	for (auto a : mesh.cell_diam) {
		std::cout << a << " ";
	}
	std::cout << "\n";

	std::cout << "face_centers" << std::endl;
	for (auto i : mesh.face_centers) {
		for (auto j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

/*

	for (int jc = 0; jc < mesh.nc; ++jc) {
		std::cout << "mesh.cell_volume["<<jc<<"] = " << mesh.cell_volumes[jc] << std::endl;
	}

	for (int jf = 0; jf < mesh.nf; ++jf) {
		std::cout << "mesh.face_areas["<<jf<<"] = " << mesh.face_areas[jf] << std::endl;
		std::cout << "mesh.face_normals["<<jf<<"][0, 1, 2] = "
				<< mesh.face_normals[jf][0] << " "
				<< mesh.face_normals[jf][1] << " "
				<< mesh.face_normals[jf][2] << std::endl;
	}

//	std::cout << "Vertices" << std::endl;
//	for (int i = 0; i < (mesh.bcface_vert_lists.size() / 100); ++i) {
//		std::set <int> face = mesh.bcface_vert_set_lists[i];
//		for (auto it = mesh.bcface_vert_lists[i].begin(); it != mesh.bcface_vert_lists[i].end(); ++it)
//				std::cout << " " << *it;
//		std::cout << "\n";
//	}

	std::cout << "\nCell_face_normal_direction\n";
	for (auto cell : mesh.cell_face_normal_direction) {
		for (auto dir : cell) {
			std::cout << dir << " ";
		}
		std::cout << "\n";
	}

	std::vector < std::vector <double> > tetra;
	tetra.push_back(std::vector <double> {0, 0, 0});
	tetra.push_back(std::vector <double> {2, 5, -3});
	tetra.push_back(std::vector <double> {1, 4, -2});
	tetra.push_back(std::vector <double> {-7, 3, 0});
	std::cout << "\n" << mesh.compute_tetra_volume(tetra) << "\n";

*/

	return 0;
}
