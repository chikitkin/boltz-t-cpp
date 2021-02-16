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
#include "mesh.h"
using namespace std;

double Mesh::compute_tetra_volume(vector < vector < double > > tetra) {
	// tetra - (4, 3)
	double a00 = tetra[1][0] - tetra[0][0];
	double a01 = tetra[1][1] - tetra[0][1];
	double a02 = tetra[1][2] - tetra[0][2];
	double a10 = tetra[2][0] - tetra[0][0];
	double a11 = tetra[2][1] - tetra[0][1];
	double a12 = tetra[2][2] - tetra[0][2];
	double a20 = tetra[3][0] - tetra[0][0];
	double a21 = tetra[3][1] - tetra[0][1];
	double a22 = tetra[3][2] - tetra[0][2];

	return ((a00 * a11 * a22 + a01 * a12 * a20 + a02 * a21 * a10) -
			(a02 * a11 * a20 + a00 * a12 * a21 + a01 * a10 * a22)) / 6.0;
}

vector < vector < int > > Mesh::get_faces_for_cell(int ic) {

	vector < vector < int > > faces;
	faces.reserve(6);
	vector <int> verts = vert_list_for_cell[ic];

	faces.push_back(vector <int> {verts[0], verts[1], verts[5], verts[4]});
	faces.push_back(vector <int> {verts[1], verts[3], verts[7], verts[5]});
	faces.push_back(vector <int> {verts[3], verts[2], verts[6], verts[7]});
	faces.push_back(vector <int> {verts[2], verts[0], verts[4], verts[6]});
	faces.push_back(vector <int> {verts[1], verts[0], verts[2], verts[3]});
	faces.push_back(vector <int> {verts[4], verts[5], verts[7], verts[6]});

	return faces;
}

void Mesh::read_starcd(const string& path, const double scale) {
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
		istringstream iss(line);
		int n, v0, v1, v2, v3, type;
		if (!(iss >> n >> v0 >> v1 >> v2 >> v3 >> type)) { break; }
		vector <int> bcface;
		set <int> bcface_set;
		bcface.push_back(v0 - 1);
		bcface_set.insert(v0 - 1);
		bcface.push_back(v1 - 1);
		bcface_set.insert(v1 - 1);
		bcface.push_back(v2 - 1);
		bcface_set.insert(v2 - 1);
		bcface.push_back(v3 - 1);
		bcface_set.insert(v3 - 1);
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
			// Convert order
			vert_list_for_cell[i][0] = v3;
			vert_list_for_cell[i][1] = v2;
			vert_list_for_cell[i][2] = v7;
			vert_list_for_cell[i][3] = v6;
			vert_list_for_cell[i][4] = v0;
			vert_list_for_cell[i][5] = v1;
			vert_list_for_cell[i][6] = v4;
			vert_list_for_cell[i][7] = v5;
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
		double x0, x1, x2;
		if (!(iss >> n >> x0 >> x1 >> x2)) { break; }
		vector <double> vert;
		vert.push_back(scale * x0);
		vert.push_back(scale * x1);
		vert.push_back(scale * x2);
		vert_coo.push_back(vert);
	}
	vrt_data.close();
	vert_coo.shrink_to_fit();
	cout << "Number of vertices = " << nv << endl;
/*
	Calculate cell centers - arithmetic mean of vertises' coordinates
*/
	cell_center_coo.resize(nc, vector <double>(3));
	for (int i = 0; i < nc; ++i) {
		double x0 = 0;
		double x1 = 0;
		double x2 = 0;
		for (int j = 0; j < max_vert_in_cell; ++j) {
			x0 += vert_coo[vert_list_for_cell[i][j]][0];
			x1 += vert_coo[vert_list_for_cell[i][j]][1];
			x2 += vert_coo[vert_list_for_cell[i][j]][2];
		}
		cell_center_coo[i][0] = (x0 / max_vert_in_cell);
		cell_center_coo[i][1] = (x1 / max_vert_in_cell);
		cell_center_coo[i][2] = (x2 / max_vert_in_cell);
	}
/*
	Calculate volume of each cell
*/
	cell_volumes.resize(nc, 0.0);
	for (int ic = 0; ic < nc; ++ic) {
		vector < vector < int > > faces = get_faces_for_cell(ic);
		// Loop over faces, for each face construct 4 tetras
		// and compute their volumes
		for (int jf = 0; jf < 6; ++jf) {
			vector <double> face_center = {0.0, 0.0, 0.0};
			for (int kf = 0; kf < 4; ++kf) {
				face_center[0] += (vert_coo[faces[jf][kf]][0] / 4.0);
				face_center[1] += (vert_coo[faces[jf][kf]][1] / 4.0);
				face_center[2] += (vert_coo[faces[jf][kf]][2] / 4.0);
			}
			vector <double> x1 = vert_coo[faces[jf][0]];
			vector <double> x2 = vert_coo[faces[jf][1]];
			vector <double> x3 = vert_coo[faces[jf][2]];
			vector <double> x4 = vert_coo[faces[jf][3]];
			vector <vector < double > > tetra1;
			tetra1.push_back(cell_center_coo[ic]);
			tetra1.push_back(x1);
			tetra1.push_back(x2);
			tetra1.push_back(face_center);
			cell_volumes[ic] += compute_tetra_volume(tetra1);
			vector <vector < double > > tetra2;
			tetra2.push_back(cell_center_coo[ic]);
			tetra2.push_back(x2);
			tetra2.push_back(x3);
			tetra2.push_back(face_center);
			cell_volumes[ic] += compute_tetra_volume(tetra2);
			vector <vector < double > > tetra3;
			tetra3.push_back(cell_center_coo[ic]);
			tetra3.push_back(x3);
			tetra3.push_back(x4);
			tetra3.push_back(face_center);
			cell_volumes[ic] += compute_tetra_volume(tetra3);
			vector <vector < double > > tetra4;
			tetra4.push_back(cell_center_coo[ic]);
			tetra4.push_back(x4);
			tetra4.push_back(x1);
			tetra4.push_back(face_center);
			cell_volumes[ic] += compute_tetra_volume(tetra4);
		}
	}
/*
	Construct for each vertex list of cells to which it belongs
*/
	cell_list_for_vertex.resize(nv, vector <int> {-1, -1, -1, -1, -1, -1, -1, -1}); // it may be > 8!
	cell_num_for_vertex.resize(nv, 0); // Number of cells, adjacent to each vertex
	for (int ic = 0; ic < nc; ++ic) {
		for (int jv = 0; jv < max_vert_in_cell; ++jv) {
			int vert = vert_list_for_cell[ic][jv]; // global vertex index
			cell_list_for_vertex[vert][cell_num_for_vertex[vert]] = ic;
			cell_num_for_vertex[vert] += 1;
		}
	}
/*
	Construct for each cell list of neighboring cells
*/
	cell_neighbors_list.resize(nc, vector <int> {-1, -1, -1, -1, -1, -1});
	face_vert_list.resize(6 * nc, vector <int> {0, 0, 0, 0});
	cell_face_list.resize(nc, vector <int> {-1, -1, -1, -1, -1, -1});
	nf = 0;
	for (int ic = 0; ic < nc; ic++) {
		vector < vector < int > > faces = get_faces_for_cell(ic);
		for (int jf = 0; jf < 6; ++jf) {
			if  (cell_neighbors_list[ic][jf] >= 0) { // if face is already assigned - skip
				continue;
			}
			face_vert_list[nf] = faces[jf]; // add face to global list
			cell_face_list[ic][jf] = nf;
			// loop over all vertices in face
			for (int iv = 0; iv < 4; ++iv) {
				// loop over all cells containing this vertex
				for (int kc = 0; kc < cell_num_for_vertex[faces[jf][iv]]; ++kc) {
					int icell = cell_list_for_vertex[faces[jf][iv]][kc];
					vector < vector < int > > faces_neigh = get_faces_for_cell(icell);
					// Now compare these faces with face
					vector <int> faces_sorted(4);
					vector <int> faces_neigh_sorted(4);
					for (int lf = 0; lf < 6; ++lf) {
						partial_sort_copy(faces[jf].begin(), faces[jf].end(), faces_sorted.begin(), faces_sorted.end());
						partial_sort_copy(faces_neigh[lf].begin(), faces_neigh[lf].end(), faces_neigh_sorted.begin(), faces_neigh_sorted.end());
						if (faces_sorted == faces_neigh_sorted) {
							cell_face_list[icell][lf] = nf;
							cell_neighbors_list[ic][jf] = icell;
							cell_neighbors_list[icell][lf] = ic;
						}
					}
				}
			}
			nf += 1;
		}
	}

	cout << "Number of faces = " << nf << endl;
	face_vert_list.resize(nf); // exclude extra rows
	cout << "Sum of volumes = " << accumulate(cell_volumes.begin(), cell_volumes.end(), 0.0) << endl;
/*
	Compute face areas and normals
*/
	face_areas.resize(nf);
	face_normals.resize(nf, vector <double>(3));
	for (int jf = 0; jf < nf; ++jf) {
		vector <int> verts = face_vert_list[jf];
		vector <double> verts_coo0 = vert_coo[verts[0]];
		vector <double> verts_coo1 = vert_coo[verts[1]];
		vector <double> verts_coo2 = vert_coo[verts[2]];
		vector <double> verts_coo3 = vert_coo[verts[3]];

		vector <double> vec1(3);
		vec1[0] = 0.5 * (verts_coo2[0] + verts_coo1[0]) - 0.5 * (verts_coo0[0] + verts_coo3[0]);
		vec1[1] = 0.5 * (verts_coo2[1] + verts_coo1[1]) - 0.5 * (verts_coo0[1] + verts_coo3[1]);
		vec1[2] = 0.5 * (verts_coo2[2] + verts_coo1[2]) - 0.5 * (verts_coo0[2] + verts_coo3[2]);

		vector <double> vec2(3);
		vec2[0] = 0.5 * (verts_coo3[0] + verts_coo2[0]) - 0.5 * (verts_coo1[0] + verts_coo0[0]);
		vec2[1] = 0.5 * (verts_coo3[1] + verts_coo2[1]) - 0.5 * (verts_coo1[1] + verts_coo0[1]);
		vec2[2] = 0.5 * (verts_coo3[2] + verts_coo2[2]) - 0.5 * (verts_coo1[2] + verts_coo0[2]);

		face_areas[jf] = sqrt(pow(vec1[1]*vec2[2] - vec1[2]*vec2[1], 2) +
				pow(vec1[2]*vec2[0] - vec1[0]*vec2[2], 2) +
				pow(vec1[0]*vec2[1] - vec1[1]*vec2[0], 2));

		// TODO: Complicated procedure to overcome problems when area is tiny
	}

	// face centers
	for (int jf = 0; jf < nf; ++jf) {
		vector <int> face_verts = face_vert_list[jf];
		vector <double> face_center = {0.0, 0.0, 0.0};
		for (int i = 0; i < 4; ++i) {
			face_center[0] += vert_coo[face_verts[i]][0] / 4.0;
			face_center[1] += vert_coo[face_verts[i]][1] / 4.0;
			face_center[2] += vert_coo[face_verts[i]][2] / 4.0;
		}
		face_centers.push_back(face_center);
	}
	face_centers.shrink_to_fit();

/*
	Compute orientation of face normals with respect to each cell

	+1 - outer normal, -1 - inner normal (directed in cell)
*/
	cell_face_normal_direction.resize(nc, vector <int>(6));
	for (int ic = 0; ic < nc; ++ic) {
		for (int jf = 0; jf < 6; ++jf) {
			int face = cell_face_list[ic][jf];
			vector <double> face_normal = face_normals[face];
			// Compute vector from cell center to center of face
			vector <double> vec = {face_centers[face][0] - cell_center_coo[ic][0],
					face_centers[face][1] - cell_center_coo[ic][1], face_centers[face][2] - cell_center_coo[ic][2]};
			double dot_prod = vec[0] * face_normal[0] + vec[1] * face_normal[1] + vec[2] * face_normal[2];
			if (dot_prod >= 0) {
				cell_face_normal_direction[ic][jf] = +1;
			}
			else {
				cell_face_normal_direction[ic][jf] = -1;
			}
		}
	}
	// TODO: comment
	bound_face_info.resize(nbf, vector <int>(3));
	for (int ibf = 0; ibf < nbf; ++ibf) {
		for (int jf = 0; jf < nf; ++jf) {
			vector <int> bcface_vert_sorted(4);
			vector <int> face_vert_sorted(4);
			partial_sort_copy(bcface_vert_lists[ibf].begin(), bcface_vert_lists[ibf].end(), bcface_vert_sorted.begin(), bcface_vert_sorted.end());
			partial_sort_copy(face_vert_list[jf].begin(), face_vert_list[jf].end(), face_vert_sorted.begin(), face_vert_sorted.end());
			if (bcface_vert_sorted == face_vert_sorted) {
				bound_face_info[ibf][0] = jf;
			}
		}
		for (int ic = 0; ic < nc; ++ic) {
			for (int jf = 0; jf < 6; ++jf) {
				if (cell_face_list[ic][jf] == bound_face_info[ibf][0]) {
					bound_face_info[ibf][2] = cell_face_normal_direction[ic][jf];
				}
			}
		}
		bound_face_info[ibf][1] = bcface_bctype[ibf];
	}
	// TODO: comment
	cell_diam.resize(nc, 0.0);
	for (int ic = 0; ic < nc; ++ic) {
		vector <double> face_diam(6);
		for (int jf = 0; jf < 6; ++jf) {
			int face = cell_face_list[ic][jf];
			vector <double> vec = {face_centers[face][0] - cell_center_coo[ic][0],
					face_centers[face][1] - cell_center_coo[ic][1], face_centers[face][2] - cell_center_coo[ic][2]};
			face_diam[jf] = 2.0 * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		}
		int dist = distance(face_diam.begin(), min_element(face_diam.begin(), face_diam.end()));
		cell_diam[ic] = face_diam[dist];
	}
	// end of function
}

void Mesh::write_tecplot(vector < vector <double> > data, string filename,
		vector <string> var_names, double time) {

	int nv = data[0].size();
	ofstream file;
	file.open(filename);
	if (file.fail())
	{
	    cout << "Could not open " << filename << endl;
	    exit(1);
	}
	file << "TITLE = \"VolumeData\"\n";
	file << "VARIABLES = \"x\" \"y\" \"z\" ";
	for (int iv = 0; iv < nv; ++iv) {
		file << " \"" << var_names[iv] << "\" ";
	}
	file << "\n";
	file << "ZONE T= my_zone, SolutionTime = " << time <<
            ", DATAPACKING=Block, ZONETYPE=FEBRICK Nodes= " << nv <<
           " Elements= " << nc;
	file << " VarLocation=([4-" << 3+nv << "]=CellCentered)";
	// Write vertices' coo;
	for (int i = 0; i < 3; ++i) {
		for (int iv = 0; iv < nv; ++iv) {
			file << vert_coo[iv][i] << "\n"; // TODO: format
		}
	}
	// Write values of variables
	for (int i = 0; i < nv; ++i) {
		for (int ic = 0; ic < nc; ++ic) {
			file << data[ic][i] << "\n"; // TODO: format
		}
	}
	// Write cell-to-vertices connectivity
	for (int ic = 0; ic < nc; ++ic) {
		vector <int> verts = vert_list_for_cell[ic]; // TODO: format
		file << verts[4] + 1 << " ";
		file << verts[5] + 1 << " ";
		file << verts[1] + 1 << " ";
		file << verts[0] + 1 << " ";
		file << verts[6] + 1 << " ";
		file << verts[7] + 1 << " ";
		file << verts[3] + 1 << " ";
		file << verts[2] + 1 << " ";
		file << "\n";
	}
}

int main() {

	Mesh mesh;

	mesh.read_starcd("/home/egor/git/boltz-t-cpp/mesh-cyl/");
/*
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
*/

	cout << "\n";
	for (auto c : mesh.face_areas) {
		cout << c << " ";
	}

	vector < vector < double > > data;
	data.resize(mesh.nv, {0.0, 0.0, 0.0});

	vector <string> var_names {"A", "B", "C"};

	mesh.write_tecplot(data, "file.dat", var_names);

/*
	vector < vector <double> > tetra;
	tetra.push_back(vector <double> {0, 0, 0});
	tetra.push_back(vector <double> {2, 5, -3});
	tetra.push_back(vector <double> {1, 4, -2});
	tetra.push_back(vector <double> {-7, 3, 0});
	cout << "\n" << mesh.compute_tetra_volume(tetra) << "\n";
*/
	return 0;
}
