#include "header.h"
#include "mesh.h"

int main() {

	Mesh mesh("/home/egor/git/boltz-t-cpp/mesh-1d/");

	std::ofstream file;
	file.open("/home/egor/git/boltz-t-cpp/mesh.txt");

	file << "Face normals" << "\n";
	for (auto face : mesh.face_normals) {
		for (auto info : face) {
			file << info << " ";
		}
		file << "\n";
	}

	file << "bcface_bctypes" << "\n";
	for (auto type : mesh.bcface_bctype) {
		file << type << "\n";
	}

	file << "cell_face_list" << "\n";
	for (auto cell : mesh.cell_face_list) {
		for (auto face : cell) {
			file << face << " ";
		}
		file << "\n";
	}

	file << "cell_face_normal_direction" << "\n";
	for (auto cell : mesh.cell_face_normal_direction) {
		for (auto dir : cell) {
			file << dir << " ";
		}
		file << "\n";
	}

	file << "face_centers" << "\n";
	for (auto face : mesh.face_centers) {
		for (auto dim : face) {
			file << dim << " ";
		}
		file << "\n";
	}

	file << "bound_face_info" << "\n";
	for (auto face : mesh.bound_face_info) {
		for (auto i : face) {
			file << i << " ";
		}
		file << "\n";
	}
/*
	for (int jf = 0; jf < mesh.nf; ++jf) {
		std::cout << "mesh.face_areas["<<jf<<"] = " << mesh.face_areas[jf] << std::endl;
		std::cout << "mesh.face_normals["<<jf<<"][0] = " << mesh.face_normals[jf][0] << std::endl;
		std::cout << "mesh.face_normals["<<jf<<"][1] = " << mesh.face_normals[jf][1] << std::endl;
		std::cout << "mesh.face_normals["<<jf<<"][2] = " << mesh.face_normals[jf][2] << std::endl;
	}
//	mesh.read_starcd("/home/egor/git/boltz-t-cpp/mesh-cyl/");

	std::cout << "Vertices" << std::endl;
	for (int i = 0; i < (mesh.bcface_vert_lists.size() / 100); ++i) {
		set <int> face = mesh.bcface_vert_set_lists[i];
		for (auto it = mesh.bcface_vert_lists[i].begin(); it != mesh.bcface_vert_lists[i].end(); ++it)
				std::cout << " " << *it;
		std::cout << "\n";
	}
	std::cout << "Boundary condition type" << std::endl;
	for (int n : mesh.bcface_bctype) {
		std::cout << n << " ";
	}
	std::cout << "\nBf for each bc\n";
	for (int i = 0; i < (mesh.bf_for_each_bc.size()); ++i) {
		std::cout << mesh.bf_for_each_bc[i][0] << std::endl;
	}

	std::cout << "\nCell center coo\n";
	for (int i = 0; i < (mesh.nc); ++i) {
		std::cout << mesh.cell_center_coo[i][1] << " ";
	}

	std::cout << "\n";
	for (int i = 0; i < (mesh.nc); ++i) {
		std::cout << mesh.cell_center_coo[i][1] << " ";
	}
	std::cout << "\n";
	for (int i = 0; i < (mesh.nc); ++i) {
		std::cout << mesh.cell_center_coo[i][2] << " ";
	}

	std::cout << "Cell volumes" << std::endl;
	for (double c : mesh.cell_volumes) {
		std::cout << c << " ";
	}

	std::cout << "\nFace areas" << std::endl;
	for (double c : mesh.face_areas) {
		std::cout << c << " ";
	}


	std::cout << "\n";
	for (auto c : mesh.face_areas) {
		std::cout << c << " ";
	}

	std::vector < std::vector < double > > data;
	data.resize(mesh.nv, {0.0, 0.0, 0.0});

	std::vector <string> var_names {"A", "B", "C"};

	mesh.write_tecplot(data, "file.dat", var_names);


	std::vector < std::vector <double> > tetra;
	tetra.push_back(std::vector <double> {0, 0, 0});
	tetra.push_back(std::vector <double> {2, 5, -3});
	tetra.push_back(std::vector <double> {1, 4, -2});
	tetra.push_back(std::vector <double> {-7, 3, 0});
	cout << "\n" << mesh.compute_tetra_volume(tetra) << "\n";
*/
	return 0;
}
