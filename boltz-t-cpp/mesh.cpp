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
//	TODO: set or vector?
	vector < set < int > > bcface_vert_lists;
	vector < int > bcface_bctype;
	vector < vector < int > > verts;

	void read_starcd(const string& path, const float scale = 1.0) {
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
			set <int> bcface;
			bcface.insert(v1 - 1);
			bcface.insert(v2 - 1);
			bcface.insert(v3 - 1);
			bcface.insert(v4 - 1);
			bcface_vert_lists.push_back(bcface);
			bcface_bctype.push_back(type - 1);
		}
		bnd_data.close();
		bcface_vert_lists.shrink_to_fit();
		bcface_bctype.shrink_to_fit();

		int nbf = bcface_bctype.size();

		cout << "Number of boundary faces = " << nbf << endl;
/*
		Construct list of boundary faces indices for each bctype
*/
		set <int> bc(bcface_bctype.begin(), bcface_bctype.end());
		int nbc = bc.size();
		cout << "Number of boundary conditions = " << nbc << endl;

//		TODO: bf_for_each_bc

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
			int n, v1, v2, v3, v4, v5, v6, v7, v8;
			if (!(iss >> n >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8)) { break; }
			vector <int> cell;
			cell.push_back(v1 - 1);
			cell.push_back(v2 - 1);
			cell.push_back(v3 - 1);
			cell.push_back(v4 - 1);
			cell.push_back(v5 - 1);
			cell.push_back(v6 - 1);
			cell.push_back(v7 - 1);
			cell.push_back(v8 - 1);
			verts.push_back(cell);
		}
		verts.shrink_to_fit();
		cel_data.close();
		int nc = 0;
		for (int i = 0; i < verts.size(); ++i) {
//			Check, is it "shell" cell?
			if ((verts[i][4] != -1) && (verts[i][5] != -1) &&
					(verts[i][6] != -1) && (verts[i][7] != -1)) {
				nc += 1; // else it is shell
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
		int nv = 0;
		while (getline(vrt_data, line)) {
			nv += 1;
		}
		vrt_data.close();
		cout << "Number of vertices = " << nv << endl;

	}
};

int main() {

	Mesh mesh;

	mesh.read_starcd("/home/egor/eclipse-workspace/boltz-t-cpp/mesh-cyl/");
	cout << "Vertices" << endl;
	for (int i = 0; i < (mesh.bcface_vert_lists.size() / 100); ++i) {
		set <int> face = mesh.bcface_vert_lists[i];
		for (auto it = mesh.bcface_vert_lists[i].begin(); it != mesh.bcface_vert_lists[i].end(); ++it)
		        cout << " " << *it;
		cout << "\n";
	}
	cout << "Boundary condition type" << endl;
	for (int i = 0; i < mesh.bcface_bctype.size(); ++i)
		cout << mesh.bcface_bctype[i] << " ";

	return 0;
}
