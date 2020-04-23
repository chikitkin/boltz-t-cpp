#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>
using namespace std;

class Mesh{
public:
	int nbf;
	int nbc;
	int nc;
	int nv;

	vector < vector < int > > bcface_vert_lists; // TODO: Needs ()?
	vector < set < int > > bcface_vert_set_lists;
	vector < int > bcface_bctype;
	vector < vector < int > > vert_list_for_cell;
	vector < vector < int > > bf_for_each_bc;
	vector < vector < double > > vert_coo;
	vector < vector < double > > cell_center_coo;

	void read_starcd(const string& path, const double scale = 1.0) {
		int max_vert_in_face = 4;
		int max_vert_in_cell = 8;
/*
		Read vertex list and bc type for each boundary face
*/
		ifstream bnd_data;
		bnd_data.open(path + "star.bnd");
		if (bnd_data.fail())
		{
		    cout << "Could not open star.bnd" << endl;
		    exit(1);
		}
		string line;
		while (getline(bnd_data, line)) {
			istringstream iss(line); // TODO: is it good?
			int n, v1, v2, v3, v4, type;
			if (!(iss >> n >> v1 >> v2 >> v3 >> v4 >> type)) { break; }
			vector <int> bcface;
			set <int> bcface_set;
			bcface.push_back(v1 - 1);
			bcface_set.insert(v1 - 1);
			bcface.push_back(v2 - 1);
			bcface_set.insert(v2 - 1);
			bcface.push_back(v3 - 1);
			bcface_set.insert(v3 - 1);
			bcface.push_back(v4 - 1);
			bcface_set.insert(v4 - 1);
//			TODO: SHRINK_TO_FIT??
//			TODO: SHRINK_TO_FIT??
			bcface_vert_lists.push_back(bcface);
			bcface_vert_set_lists.push_back(bcface_set);
			bcface_bctype.push_back(type - 1);
		}
		bnd_data.close();
		bcface_vert_set_lists.shrink_to_fit();
		bcface_vert_lists.shrink_to_fit();
		bcface_bctype.shrink_to_fit();

		nbf = bcface_bctype.size();

		cout << "Number of boundary faces = " << nbf << endl;
/*
		Construct list of boundary faces indices for each bctype
*/
		set <int> bc(bcface_bctype.begin(), bcface_bctype.end());
		nbc = bc.size();
		cout << "Number of boundary conditions = " << nbc << endl;

		bf_for_each_bc.resize(nbc, vector <int>());
		for (int i = 0; i < nbf; ++i) {
			bf_for_each_bc[bcface_bctype[i]].push_back(i);
		}
		bf_for_each_bc.shrink_to_fit(); // TODO: Doesnt need?
/*
		Count number of cells
*/
		ifstream cel_data;
		cel_data.open(path + "star.cel");
		if (cel_data.fail())
		{
			cout << "Could not open star.cel" << endl;
			exit(1);
		}
		while (getline(cel_data, line)) {
			istringstream iss(line);
			int n, v0, v1, v2, v3, v4, v5, v6, v7;
			if (!(iss >> n >> v0 >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7)) { break; }
			vector <int> cell;
			cell.push_back(v0 - 1);
			cell.push_back(v1 - 1);
			cell.push_back(v2 - 1);
			cell.push_back(v3 - 1);
			cell.push_back(v4 - 1);
			cell.push_back(v5 - 1);
			cell.push_back(v6 - 1);
			cell.push_back(v7 - 1);
//			TODO: SHRINK_TO_FIT??
			vert_list_for_cell.push_back(cell);
		}
		cel_data.close();
		vert_list_for_cell.shrink_to_fit();
		nc = 0;
		for (int i = 0; i < vert_list_for_cell.size(); ++i) {
//			Check, is it "shell" cell?
			if ((vert_list_for_cell[i][4] != -1) && (vert_list_for_cell[i][5] != -1) &&
					(vert_list_for_cell[i][6] != -1) && (vert_list_for_cell[i][7] != -1)) {
				nc += 1; // else it is shell
				int v0 = vert_list_for_cell[i][0];
				int v1 = vert_list_for_cell[i][1];
				int v2 = vert_list_for_cell[i][2];
				int v3 = vert_list_for_cell[i][3];
				int v4 = vert_list_for_cell[i][4];
				int v5 = vert_list_for_cell[i][5];
				int v6 = vert_list_for_cell[i][6];
				int v7 = vert_list_for_cell[i][7];
				// Convert order to StarCD
				vert_list_for_cell[i][0] = v4;
				vert_list_for_cell[i][1] = v5;
				vert_list_for_cell[i][2] = v7;
				vert_list_for_cell[i][3] = v6;
				vert_list_for_cell[i][4] = v0;
				vert_list_for_cell[i][5] = v1;
				vert_list_for_cell[i][6] = v3;
				vert_list_for_cell[i][7] = v2;
				// Convert order to Gambit
				vert_list_for_cell[i][0] = v6;
				vert_list_for_cell[i][1] = v7;
				vert_list_for_cell[i][2] = v2;
				vert_list_for_cell[i][3] = v3;
				vert_list_for_cell[i][4] = v4;
				vert_list_for_cell[i][5] = v5;
				vert_list_for_cell[i][6] = v0;
				vert_list_for_cell[i][7] = v1;
			}
		}
		cout << "Number of cells = " << nc << endl;
/*
		Count number of vertices
*/
		ifstream vrt_data;
		vrt_data.open(path + "star.vrt");
		if (vrt_data.fail())
		{
			cout << "Could not open star.vrt" << endl;
			exit(1);
		}
		nv = 0;
		while (getline(vrt_data, line)) {
			nv += 1;
			istringstream iss(line);
			int n;
			double x1, x2, x3;
			if (!(iss >> n >> x1 >> x2 >> x3)) { break; }
			vector <double> vert;
			vert.push_back(scale * x1);
			vert.push_back(scale * x2);
			vert.push_back(scale * x3);
			vert_coo.push_back(vert);
		}
		vrt_data.close();
		vert_coo.shrink_to_fit();
		cout << "Number of vertices = " << nv << endl;
/*
		Calculate cell centers - arithmetic mean of vertises' coordinates
*/
		cell_center_coo.resize(nc, vector<double>());
		for (int i = 0; i < nc; ++i) {
			double x0 = 0;
			double x1 = 0;
			double x2 = 0;
			for (int j = 0; j < max_vert_in_cell; ++j) {
				x0 += vert_coo[vert_list_for_cell[i][j]][0];
				x1 += vert_coo[vert_list_for_cell[i][j]][1];
				x2 += vert_coo[vert_list_for_cell[i][j]][2];
			}
			cell_center_coo[i].push_back(x0 / max_vert_in_cell);
			cell_center_coo[i].push_back(x1 / max_vert_in_cell);
			cell_center_coo[i].push_back(x2 / max_vert_in_cell);
		}
		// TODO: Need shrink?
	}
};

int main() {

	Mesh mesh;

	mesh.read_starcd("/home/egor/eclipse-workspace/boltz-t-cpp/mesh-cyl/");
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
	cout << "\nbf_for_each_bc\n";
	for (int i = 0; i < (mesh.bf_for_each_bc.size()); ++i) {
		cout << mesh.bf_for_each_bc[i][0] << endl;
	}

	return 0;
}
