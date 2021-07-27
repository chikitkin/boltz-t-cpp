#include "mesh.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <numeric>
using namespace std;

int main() {

	Mesh mesh("/home/egor/git/boltz-t-cpp/mesh-1d/");

	ofstream file;
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
		cout << "mesh.face_areas["<<jf<<"] = " << mesh.face_areas[jf] << endl;
		cout << "mesh.face_normals["<<jf<<"][0] = " << mesh.face_normals[jf][0] << endl;
		cout << "mesh.face_normals["<<jf<<"][1] = " << mesh.face_normals[jf][1] << endl;
		cout << "mesh.face_normals["<<jf<<"][2] = " << mesh.face_normals[jf][2] << endl;
	}
//	mesh.read_starcd("/home/egor/git/boltz-t-cpp/mesh-cyl/");

	cout << "Vertices" << endl;
	for (int i = 0; i < (mesh.bcface_vert_lists.size() / 100); ++i) {
		set <int> face = mesh.bcface_vert_set_lists[i];
		for (auto it = mesh.bcface_vert_lists[i].begin(); it != mesh.bcface_vert_lists[i].end(); ++it)
				cout << " " << *it;
		cout << "\n";
	}
	cout << "Boundary condition type" << endl;
	for (int n : mesh.bcface_bctype) {
		cout << n << " ";
	}
	cout << "\nBf for each bc\n";
	for (int i = 0; i < (mesh.bf_for_each_bc.size()); ++i) {
		cout << mesh.bf_for_each_bc[i][0] << endl;
	}

	cout << "\nCell center coo\n";
	for (int i = 0; i < (mesh.nc); ++i) {
		cout << mesh.cell_center_coo[i][1] << " ";
	}

	cout << "\n";
	for (int i = 0; i < (mesh.nc); ++i) {
		cout << mesh.cell_center_coo[i][1] << " ";
	}
	cout << "\n";
	for (int i = 0; i < (mesh.nc); ++i) {
		cout << mesh.cell_center_coo[i][2] << " ";
	}

	cout << "Cell volumes" << endl;
	for (double c : mesh.cell_volumes) {
		cout << c << " ";
	}

	cout << "\nFace areas" << endl;
	for (double c : mesh.face_areas) {
		cout << c << " ";
	}


	cout << "\n";
	for (auto c : mesh.face_areas) {
		cout << c << " ";
	}

	vector < vector < double > > data;
	data.resize(mesh.nv, {0.0, 0.0, 0.0});

	vector <string> var_names {"A", "B", "C"};

	mesh.write_tecplot(data, "file.dat", var_names);


	vector < vector <double> > tetra;
	tetra.push_back(vector <double> {0, 0, 0});
	tetra.push_back(vector <double> {2, 5, -3});
	tetra.push_back(vector <double> {1, 4, -2});
	tetra.push_back(vector <double> {-7, 3, 0});
	cout << "\n" << mesh.compute_tetra_volume(tetra) << "\n";
*/
	return 0;
}
