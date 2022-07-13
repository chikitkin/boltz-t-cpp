#include "mesh.h"
#include "header.h"

// for unordered map
struct FaceHasher {
    int operator()(const std::vector<int> &face) const {
        int hash = 0;
        for(auto &i : face) {
            hash += i;
        }
        return hash;
    }
};

void Mesh::readVerts(std::ifstream &data, double scale) {
	std::string line;
	getline(data, line);
	std::istringstream(line) >> nVerts;
	std::cout << "Number of vertices = " << nVerts << std::endl;
	vertsCoo.reserve(nVerts);
	for (int iv = 0; iv < nVerts; ++iv) {
		getline(data, line);
		std::istringstream iss(line);
		double x0, x1, x2;
		int tag;
		if (!(iss >> x0 >> x1 >> x2 >> tag)) {
			std::cout << "Wrong vertice line!" << std::endl;
		}
		std::vector <double> vert;
		vert.push_back(scale * x0);
		vert.push_back(scale * x1);
		vert.push_back(scale * x2);
		vertsCoo.push_back(vert);
	}
}
//             0_____1
//            / |   / |
//           4__|__5  |
//           |  2__|__3
//           | /   | /
//           6_____7
void Mesh::readHexa(std::ifstream &data) {
	std::string line;
	getline(data, line);
	int nHexa;
	std::istringstream(line) >> nHexa;
	std::cout << "Number of hexa cells = " << nHexa << std::endl;
	if (cellVerts.empty()) {
		nCells = nHexa;
	}
	else {
		nCells += nHexa;
	}
	cellVerts.reserve(nCells);
	cellTypes.resize(nCells, HEXA);
	for (int i = 0; i < nHexa; ++i) {
		getline(data, line);
		std::istringstream iss(line);
		int v0, v1, v2, v3, v4, v5, v6, v7, tag;
		if (!(iss >> v0 >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> tag)) {
			std::cout << "Wrong hexa line!" << std::endl;
		}
		std::vector <int> cell;
		// reorder to gambit
		cell.push_back(v7 - 1);
		cell.push_back(v6 - 1);
		cell.push_back(v3 - 1);
		cell.push_back(v2 - 1);
		cell.push_back(v4 - 1);
		cell.push_back(v5 - 1);
		cell.push_back(v0 - 1);
		cell.push_back(v1 - 1);
		cellVerts.push_back(cell);
	}
}
//                  3
//                 / \
//                /   \
//               /     \
//              /       \
//             /    0    \
//            /           \
//           1_____________2
void Mesh::readTetra(std::ifstream &data) {
	std::string line;
	getline(data, line);
	int nTetra;
	std::istringstream(line) >> nTetra;
	std::cout << "Number of tetra cells = " << nTetra << std::endl;
	if (cellVerts.empty()) {
		nCells = nTetra;
	}
	else {
		nCells += nTetra;
	}
	cellVerts.reserve(nCells);
	cellTypes.resize(nCells, TETRA);
	for (int i = 0; i < nTetra; ++i) {
		getline(data, line);
		std::istringstream iss(line);
		int v0, v1, v2, v3, tag;
		if (!(iss >> v0 >> v1 >> v2 >> v3 >> tag)) {
			std::cout << "Wrong tetra line!" << std::endl;
		}
		std::vector <int> cell;
		// reorder so counter-clock TODO
		cell.push_back(v3 - 1);
		cell.push_back(v1 - 1);
		cell.push_back(v2 - 1);
		cell.push_back(v0 - 1);
		cellVerts.push_back(cell);
	}
}

void Mesh::readBoundaryQuads(std::ifstream &data) {
	std::string line;
	getline(data, line);
	int nBoundaryQuads;
	std::istringstream(line) >> nBoundaryQuads;
	std::cout << "Number of boundary quadrilaterals = " << nBoundaryQuads << std::endl;
	if (boundaryFaces.empty()) {
		nBoundaryFaces = nBoundaryQuads;
	}
	else {
		nBoundaryFaces += nBoundaryQuads;
	}
	boundaryFaces.reserve(nBoundaryFaces);
	boundaryFaceTypes.resize(nBoundaryFaces, QUAD);
	boundaryFaceTags.reserve(nBoundaryFaces);
	for (int i = 0; i < nBoundaryQuads; ++i) {
		getline(data, line);
		std::istringstream iss(line);
		int v0, v1, v2, v3, tag;
		if (!(iss >> v0 >> v1 >> v2 >> v3 >> tag)) {
			std::cout << "Wrong quad line!" << std::endl;
		}
		std::vector <int> bcface;
		bcface.push_back(v0 - 1);
		bcface.push_back(v1 - 1);
		bcface.push_back(v2 - 1);
		bcface.push_back(v3 - 1);
		boundaryFaces.push_back(bcface);
		boundaryFaceTags.push_back(tag);
	}
}

void Mesh::readBoundaryTriangles(std::ifstream &data) {
	std::string line;
	getline(data, line);
	int nBoundaryTriangles;
	std::istringstream(line) >> nBoundaryTriangles;
	std::cout << "Number of boundary triangles = " << nBoundaryTriangles << std::endl;
	if (boundaryFaces.empty()) {
		nBoundaryFaces = nBoundaryTriangles;
	}
	else {
		nBoundaryFaces += nBoundaryTriangles;
	}
	boundaryFaces.reserve(nBoundaryFaces);
	boundaryFaceTypes.resize(nBoundaryFaces, TRIANGLE);
	boundaryFaceTags.reserve(nBoundaryFaces);
	for (int i = 0; i < nBoundaryTriangles; ++i) {
		getline(data, line);
		std::istringstream iss(line);
		int v0, v1, v2, tag;
		if (!(iss >> v0 >> v1 >> v2 >> tag)) {
			std::cout << "Wrong triangle line!" << std::endl;
		}
		std::vector <int> bcface;
		bcface.push_back(v0 - 1);
		bcface.push_back(v1 - 1);
		bcface.push_back(v2 - 1);
		boundaryFaces.push_back(bcface);
		boundaryFaceTags.push_back(tag);
	}
}

double Mesh::computeTetraVolume(std::vector < std::vector < double > > tetra) {
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

std::vector<std::vector<int>> Mesh::computeFacesOfCell(int ic) {

	std::vector <int> verts = cellVerts[ic];
	std::vector<std::vector<int>> faces;

	if (cellTypes[ic] == HEXA) {
		faces.reserve(6);
		faces.push_back(std::vector <int> {verts[0], verts[1], verts[5], verts[4]});
		faces.push_back(std::vector <int> {verts[1], verts[3], verts[7], verts[5]});
		faces.push_back(std::vector <int> {verts[3], verts[2], verts[6], verts[7]});
		faces.push_back(std::vector <int> {verts[2], verts[0], verts[4], verts[6]});
		faces.push_back(std::vector <int> {verts[1], verts[0], verts[2], verts[3]});
		faces.push_back(std::vector <int> {verts[4], verts[5], verts[7], verts[6]});
	}
	else if (cellTypes[ic] == TETRA) {
		faces.reserve(4);
		faces.push_back(std::vector <int> {verts[2], verts[1], verts[3]});
		faces.push_back(std::vector <int> {verts[0], verts[1], verts[2]});
		faces.push_back(std::vector <int> {verts[0], verts[2], verts[3]});
		faces.push_back(std::vector <int> {verts[0], verts[3], verts[1]});
	}

	return faces;
}

Mesh::Mesh(const std::string& path, double scale) {
	read(path, scale);
}

void Mesh::read(const std::string& path, double scale) {

	std::ifstream data;
	data.open(path + "mesh.mesh");
	if (data.fail())
	{
		std::cout << "Could not open mesh data" << std::endl;
		exit(1);
	}

	std::string line;
	while (getline(data, line)) {
		if (line == "End") {
			break;
		}
		else if (line == "Vertices") {
			readVerts(data, scale);
		}
		else if (line == "Hexahedra") {
			readHexa(data);
		}
		else if (line == "Quadrilaterals") {
			readBoundaryQuads(data);
		}
		else if (line == "Triangles") {
			readBoundaryTriangles(data);
		}
		else if (line == "Tetrahedra") {
			readTetra(data);
		}
	}
	data.close();
	// Create a set of boundary tags
	std::set<int> boundaryFaceTagsSet_(boundaryFaceTags.begin(), boundaryFaceTags.end());
	std::copy(boundaryFaceTagsSet_.begin(), boundaryFaceTagsSet_.end(), std::back_inserter(boundaryFaceTagsSet));
	/*
		Calculate cell centers - arithmetic mean of vertises' coordinates
	*/
	cellCenters.resize(nCells, std::vector <double>(3));
	for (int i = 0; i < nCells; ++i) {
		int nVertsInCell = cellVerts[i].size();
		double x0 = 0;
		double x1 = 0;
		double x2 = 0;
		for (int j = 0; j < nVertsInCell; ++j) {
			x0 += vertsCoo[cellVerts[i][j]][0];
			x1 += vertsCoo[cellVerts[i][j]][1];
			x2 += vertsCoo[cellVerts[i][j]][2];
		}
		cellCenters[i][0] = (x0 / nVertsInCell);
		cellCenters[i][1] = (x1 / nVertsInCell);
		cellCenters[i][2] = (x2 / nVertsInCell);
	}
	/*
		Construct for each vertex list of cells to which it belongs
	*/
	vertsCells.resize(nVerts);
	for (int ic = 0; ic < nCells; ++ic) {
		for (int jv = 0; jv < cellVerts[ic].size(); ++jv) {
			int vert = cellVerts[ic][jv]; // global vertex index
			vertsCells[vert].push_back(ic);
		}
	}
	/*
		Construct for each cell list of neighboring cells
	*/
	faces.reserve(6 * nCells); // estimated
	faceTypes.reserve(6 * nCells);
	std::unordered_map<std::vector<int>, int, FaceHasher> facesMap;
	cellNeighbors.reserve(nCells);
	cellFaces.reserve(nCells);
	for (int i = 0; i < nCells; ++i) {
		if (cellTypes[i] == HEXA) {
			cellNeighbors.push_back(std::vector <int> {-1, -1, -1, -1, -1, -1});
			cellFaces.push_back(std::vector <int> {-1, -1, -1, -1, -1, -1});
		}
		else if (cellTypes[i] == TETRA) {
			cellNeighbors.push_back(std::vector <int> {-1, -1, -1, -1});
			cellFaces.push_back(std::vector <int> {-1, -1, -1, -1});
		}
	}
	/*
		Construct all faces
	*/
	nFaces = 0; // TODO follow jf
	for (int ic = 0; ic < nCells; ic++) {
		std::vector <std::vector<int>> localCellFaces = computeFacesOfCell(ic);
		for (int jf = 0; jf < localCellFaces.size(); ++jf) {
			if (cellNeighbors[ic][jf] >= 0) { // if face is already assigned - skip
				continue;
			}
			std::vector<int> face = localCellFaces[jf];
			faces.push_back(face); // add face to global list
			std::vector<int> faceSorted(face.size());
			std::partial_sort_copy(face.begin(), face.end(), faceSorted.begin(), faceSorted.end());
			cellFaces[ic][jf] = nFaces;
			facesMap[faceSorted] = nFaces;
			if (face.size() == 4) {
				faceTypes.push_back(QUAD);
			}
			else if (face.size() == 3) {
				faceTypes.push_back(TRIANGLE);
			}
			// loop over all vertices in face
			for (int iv = 0; iv < face.size(); ++iv) {
				// loop over all cells containing this vertex
				for (int kc = 0; kc < vertsCells[face[iv]].size(); ++kc) {
					// if looking on the same cell - skip
					int icell = vertsCells[face[iv]][kc];
					if (icell == ic) {
						continue;
					}
					std::vector <std::vector<int>> neighCellFaces = computeFacesOfCell(icell);
					// Now compare these faces with face
					for (int lf = 0; lf < neighCellFaces.size(); ++lf) {
						// if num of verts in faces not equal - skip
						if (neighCellFaces[lf].size() != face.size()) {
							continue;
						}
						std::vector<int> neighCellFaceSorted(neighCellFaces[lf].size());
						std::partial_sort_copy(neighCellFaces[lf].begin(), neighCellFaces[lf].end(), neighCellFaceSorted.begin(), neighCellFaceSorted.end());
						if (faceSorted == neighCellFaceSorted) {
							cellFaces[icell][lf] = nFaces;
							cellNeighbors[ic][jf] = icell;
							cellNeighbors[icell][lf] = ic;
						}
					}
				}
			}
			++nFaces;
		}
	}
	std::cout << "Number of faces = " << nFaces << std::endl;
	faces.resize(nFaces); // exclude extra rows
	faceTypes.resize(nFaces);
	/*
		Calculate volume of each cell
	*/
	cellVolumes.resize(nCells, 0.0);
	for (int ic = 0; ic < nCells; ++ic) {
		std::vector < std::vector < int > > localFaces = computeFacesOfCell(ic);
		if (cellTypes[ic] == HEXA) {
			// Loop over faces, for each face construct 4 tetras
			// and compute their volumes
			for (int jf = 0; jf < localFaces.size(); ++jf) {
				std::vector <double> localFaceCenter = {0.0, 0.0, 0.0};
				int vertsNum = localFaces[jf].size();
				for (int kv = 0; kv < vertsNum; ++kv) {
					localFaceCenter[0] += vertsCoo[localFaces[jf][kv]][0];
					localFaceCenter[1] += vertsCoo[localFaces[jf][kv]][1];
					localFaceCenter[2] += vertsCoo[localFaces[jf][kv]][2];
				}
				localFaceCenter[0] /= static_cast<double>(vertsNum);
				localFaceCenter[1] /= static_cast<double>(vertsNum);
				localFaceCenter[2] /= static_cast<double>(vertsNum);
				std::vector <double> x1 = vertsCoo[localFaces[jf][0]];
				std::vector <double> x2 = vertsCoo[localFaces[jf][1]];
				std::vector <double> x3 = vertsCoo[localFaces[jf][2]];
				std::vector <double> x4 = vertsCoo[localFaces[jf][3]];
				std::vector <std::vector<double>> tetra0;
				tetra0.push_back(cellCenters[ic]);
				tetra0.push_back(x1);
				tetra0.push_back(x2);
				tetra0.push_back(localFaceCenter);
				cellVolumes[ic] += computeTetraVolume(tetra0);
				std::vector <std::vector<double>> tetra1;
				tetra1.push_back(cellCenters[ic]);
				tetra1.push_back(x2);
				tetra1.push_back(x3);
				tetra1.push_back(localFaceCenter);
				cellVolumes[ic] += computeTetraVolume(tetra1);
				std::vector <std::vector<double>> tetra2;
				tetra2.push_back(cellCenters[ic]);
				tetra2.push_back(x3);
				tetra2.push_back(x4);
				tetra2.push_back(localFaceCenter);
				cellVolumes[ic] += computeTetraVolume(tetra2);
				std::vector <std::vector<double>> tetra3;
				tetra3.push_back(cellCenters[ic]);
				tetra3.push_back(x4);
				tetra3.push_back(x1);
				tetra3.push_back(localFaceCenter);
				cellVolumes[ic] += computeTetraVolume(tetra3);
			}
		}
		else if (cellTypes[ic] == TETRA) {
			std::vector <double> x1 = vertsCoo[cellVerts[ic][0]];
			std::vector <double> x2 = vertsCoo[cellVerts[ic][1]];
			std::vector <double> x3 = vertsCoo[cellVerts[ic][2]];
			std::vector <double> x4 = vertsCoo[cellVerts[ic][3]];
			std::vector <std::vector<double>> tetra;
			tetra.push_back(x1);
			tetra.push_back(x2);
			tetra.push_back(x3);
			tetra.push_back(x4);
			cellVolumes[ic] = computeTetraVolume(tetra); // TODO
		}
	}
	std::cout << "Sum of volumes = " << accumulate(cellVolumes.begin(), cellVolumes.end(), 0.0) << std::endl;
	/*
		Compute face areas and normals
	*/
	faceAreas.resize(nFaces);
	faceNormals.resize(nFaces, std::vector <double>(3));
	for (int jf = 0; jf < nFaces; ++jf) {
		std::vector <int> verts = faces[jf];
		if (faceTypes[jf] == QUAD) {
			std::vector <double> vert0 = vertsCoo[verts[0]];
			std::vector <double> vert1 = vertsCoo[verts[1]];
			std::vector <double> vert2 = vertsCoo[verts[2]];
			std::vector <double> vert3 = vertsCoo[verts[3]];

			std::vector <double> vec0(3);
			vec0[0] = 0.5 * (vert2[0] + vert1[0]) - 0.5 * (vert0[0] + vert3[0]);
			vec0[1] = 0.5 * (vert2[1] + vert1[1]) - 0.5 * (vert0[1] + vert3[1]);
			vec0[2] = 0.5 * (vert2[2] + vert1[2]) - 0.5 * (vert0[2] + vert3[2]);

			std::vector <double> vec1(3);
			vec1[0] = 0.5 * (vert3[0] + vert2[0]) - 0.5 * (vert1[0] + vert0[0]);
			vec1[1] = 0.5 * (vert3[1] + vert2[1]) - 0.5 * (vert1[1] + vert0[1]);
			vec1[2] = 0.5 * (vert3[2] + vert2[2]) - 0.5 * (vert1[2] + vert0[2]);
			// TODO: посмотреть готовую функцию для векторного произведения и для нормы
			faceAreas[jf] = sqrt(pow(vec0[1]*vec1[2] - vec0[2]*vec1[1], 2) +
					pow(vec0[2]*vec1[0] - vec0[0]*vec1[2], 2) +
					pow(vec0[0]*vec1[1] - vec0[1]*vec1[0], 2));

			faceNormals[jf][0] = (vec0[1]*vec1[2] - vec0[2]*vec1[1]) / faceAreas[jf];
			faceNormals[jf][1] = (vec0[2]*vec1[0] - vec0[0]*vec1[2]) / faceAreas[jf];
			faceNormals[jf][2] = (vec0[0]*vec1[1] - vec0[1]*vec1[0]) / faceAreas[jf];
			// Complicated procedure to overcome problems when area is tiny
			double len0 = sqrt(pow(vec0[0], 2) + pow(vec0[1], 2) + pow(vec0[2], 2));
			double len1 = sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
			std::vector<double> bec0(vec0);
			std::vector<double> bec1(vec1);
			for (int ie = 0; ie < 3; ++ie) {
				bec0[ie] /= std::max(1e-13, len0);
				bec1[ie] /= std::max(1e-13, len1);
			}
			std::vector<double> normal(3);
			normal[0] = (bec0[1]*bec1[2] - bec0[2]*bec1[1]);
			normal[1] = (bec0[2]*bec1[0] - bec0[0]*bec1[2]);
			normal[2] = (bec0[0]*bec1[1] - bec0[1]*bec1[0]);
			double length = sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2));
			if (length < 1e-10) {
				std::cout << "flag = False" << std::endl;
			}
			else {
				normal[0] /= length;
				normal[1] /= length;
				normal[2] /= length;
			}
			// end of procedure
		}
		else if (faceTypes[jf] == TRIANGLE) {
			std::vector <double> vert0 = vertsCoo[verts[0]];
			std::vector <double> vert1 = vertsCoo[verts[1]];
			std::vector <double> vert2 = vertsCoo[verts[2]];

			std::vector <double> vec0(3);
			vec0[0] = (vert2[0] - vert0[0]);
			vec0[1] = (vert2[1] - vert0[1]);
			vec0[2] = (vert2[2] - vert0[2]);

			std::vector <double> vec1(3);
			vec1[0] = (vert1[0] - vert0[0]);
			vec1[1] = (vert1[1] - vert0[1]);
			vec1[2] = (vert1[2] - vert0[2]);

			faceAreas[jf] = 0.5 * sqrt(pow(vec0[1]*vec1[2] - vec0[2]*vec1[1], 2) +
					pow(vec0[2]*vec1[0] - vec0[0]*vec1[2], 2) +
					pow(vec0[0]*vec1[1] - vec0[1]*vec1[0], 2));

			faceNormals[jf][0] = (vec0[1]*vec1[2] - vec0[2]*vec1[1]) / (2.0 * faceAreas[jf]);
			faceNormals[jf][1] = (vec0[2]*vec1[0] - vec0[0]*vec1[2]) / (2.0 * faceAreas[jf]);
			faceNormals[jf][2] = (vec0[0]*vec1[1] - vec0[1]*vec1[0]) / (2.0 * faceAreas[jf]);
			// TODO: Complicated procedure to overcome problems when area is tiny
		}
	}
	/*
		Compute face centers
	*/
	faceCenters.reserve(nFaces);
	for (int jf = 0; jf < nFaces; ++jf) {
		std::vector <int> face = faces[jf];
		std::vector <double> faceCenter = {0.0, 0.0, 0.0};
		int vertsNum = face.size();
		for (int i = 0; i < vertsNum; ++i) {
			faceCenter[0] += vertsCoo[face[i]][0];
			faceCenter[1] += vertsCoo[face[i]][1];
			faceCenter[2] += vertsCoo[face[i]][2];
		}
		faceCenter[0] /= static_cast<double>(vertsNum);
		faceCenter[1] /= static_cast<double>(vertsNum);
		faceCenter[2] /= static_cast<double>(vertsNum);
		faceCenters.push_back(faceCenter);
	}
	/*
		Compute orientation of face normals with respect to each cell and cell diameters
	*/
	isOuterNormal.resize(nCells);
	cellDiameters.resize(nCells);
	for (int ic = 0; ic < nCells; ++ic) {
		std::vector <double> faceDiameters(cellFaces[ic].size(), 0.0);
		for (int jf = 0; jf < cellFaces[ic].size(); ++jf) {
			int face = cellFaces[ic][jf];
			std::vector <double> faceNormal = faceNormals[face];
			// Compute vector from cell center to center of face
			std::vector <double> vec = {faceCenters[face][0] - cellCenters[ic][0],
					faceCenters[face][1] - cellCenters[ic][1], faceCenters[face][2] - cellCenters[ic][2]};

			faceDiameters[jf] = 2.0 * sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

			double dotProduct = vec[0] * faceNormal[0] + vec[1] * faceNormal[1] + vec[2] * faceNormal[2];
			if (dotProduct >= 0.0) {
				isOuterNormal[ic].push_back(true);
			}
			else {
				isOuterNormal[ic].push_back(false);
			}
		}
		int dist = distance(faceDiameters.begin(), min_element(faceDiameters.begin(), faceDiameters.end()));
		cellDiameters[ic] = faceDiameters[dist];
	}
	/*
		Construct list of boundary faces global indices
	*/
	std::unordered_map<std::vector<int>, int, FaceHasher> boundaryFaceGlobalIndexMap; // TODO mb set?
	for (int ibf = 0; ibf < nBoundaryFaces; ++ibf) {
		std::vector<int> boundaryFace = boundaryFaces[ibf];
		std::vector<int> sortedBoundaryFace(boundaryFace.size());
		std::partial_sort_copy(boundaryFace.begin(), boundaryFace.end(), sortedBoundaryFace.begin(), sortedBoundaryFace.end());
		boundaryFaceGlobalIndexMap[sortedBoundaryFace] = facesMap[sortedBoundaryFace];
	}
	boundaryFaceIndex.reserve(nBoundaryFaces);
	for (int ibf = 0; ibf < nBoundaryFaces; ++ibf) {
		std::vector<int> boundaryFace = boundaryFaces[ibf];
		std::vector<int> sortedBoundaryFace(boundaryFace.size());
		std::partial_sort_copy(boundaryFace.begin(), boundaryFace.end(), sortedBoundaryFace.begin(), sortedBoundaryFace.end());
		boundaryFaceIndex.push_back(boundaryFaceGlobalIndexMap.at(sortedBoundaryFace));
	}
	/*
		Construct a map of bfaces global indices for tag
	*/
	for (int ibf = 0; ibf < nBoundaryFaces; ++ibf) {
		boundaryFacesForEachTag[boundaryFaceTags[ibf]].push_back(boundaryFaceIndex[ibf]);
	}
	/*
		Left and right cell according to face, -1 if no cell
	*/
	leftRightCells.resize(nFaces, std::vector<int>{-1, -1});
	for (int ic = 0; ic < nCells; ++ic) {
		for (int j = 0; j < cellFaces[ic].size(); ++j) {
			if (isOuterNormal[ic][j]) {
				leftRightCells[cellFaces[ic][j]][0] = ic;
			}
			else {
				leftRightCells[cellFaces[ic][j]][1] = ic;
			}
		}
	}
	/*
		Boundary face directions (same as in the cell)
	*/
	isOuterNormalBoundary.reserve(nFaces);
	for (int ibf = 0; ibf < nBoundaryFaces; ++ibf) {
		int jf = boundaryFaceIndex[ibf];
		if (leftRightCells[jf][0] == -1) {
			isOuterNormalBoundary[jf] = false;
		}
		else if (leftRightCells[jf][1] == -1) {
			isOuterNormalBoundary[jf] = true;
		}
		else {
			std::cout << "Boundary face with no cell" << std::endl;
			exit(1);
		}
	}
	/*
		Coloring of the mesh
	*/
	cellColors.resize(nCells, -1);
	for (int ic = 0; ic < nCells; ++ic) {
		std::vector<bool> colors(100, true); // 100 is a number of available colors
		colors[0] = false;
		for (const int& inc : cellNeighbors[ic]) {
			if (inc != -1) {
			colors[cellColors[inc] + 1] = false;
			}
		}
		for (int icolor = 1; icolor < 100; ++icolor) {
			if (colors[icolor]) {
				cellColors[ic] = icolor - 1;
				break;
			}
		}
	}
	nColors = *std::max_element(cellColors.begin(), cellColors.end()) + 1;
	coloredCells.resize(nColors, std::vector<int>());
	rcoloredCells.resize(nColors, std::vector<int>());
	std::cout << "nColors = " << nColors << std::endl;
	for (int ic = 0; ic < nCells; ++ic) {
		coloredCells[cellColors[ic]].push_back(ic);
		rcoloredCells[cellColors[ic]].push_back(ic);
	}
	for (int color = 0; color < nColors; ++color) {
		std::reverse(rcoloredCells[color].begin(), rcoloredCells[color].end());
	}
	// end of function
}

int Mesh::getNumFacesOfCell(int ic) {
	return cellFaces[ic].size();
}

int Mesh::getFacesOfCell(int ic, int jf) {
	return cellFaces[ic][jf];
}

int Mesh::getNumTags() {
	return boundaryFaceTagsSet.size();
}

// TODO fix
std::vector<int> Mesh::getFacesForTag(int tag) {
	return boundaryFacesForEachTag[tag];
}
// 1 if out, 0 else
int Mesh::getOutIndex(int ic, int j) {
	return static_cast<int>(isOuterNormal[ic][j]);
}
// 1.0 if out, -1.0 else
double Mesh::getOutSign(int ic, int j) {
	return 2.0 * static_cast<double>(isOuterNormal[ic][j]) - 1.0;
}
// 1 if out, 0 else
int Mesh::getOutIndex(int jf) {
	return static_cast<int>(isOuterNormalBoundary[jf]);
}
// 1.0 if out, -1.0 else
double Mesh::getOutSign(int jf) {
	return 2.0 * static_cast<double>(isOuterNormalBoundary[jf]) - 1.0;
}

void Mesh::divideMesh(int nParts) {

	nPartitions = nParts;

	std::vector<int> partitionSizes(nPartitions, nCells / nPartitions);
	partitionSizes[nPartitions - 1] = nCells - (nPartitions - 1) * (nCells / nPartitions);

	cellPartitions.clear();
	cellPartitions.resize(nCells, -1);
	for (int ip = 0; ip < nPartitions; ++ip) {
		int ic;
		int partitionSize = 0;
		for (int i = 0; i < nCells; ++i) {
			if (cellPartitions[i] == -1) {
				ic = i;
				cellPartitions[ic] = ip;
				++partitionSize;
				break;
			}
		}
		while (partitionSize < partitionSizes[ip]) {
			int ic_next = -1;
			for (const int& inc : cellNeighbors[ic]) {
				if (inc != -1) {
					if (cellPartitions.at(inc) == -1) {
						cellPartitions[inc] = ip; ++partitionSize;
						ic_next = inc;
						if (partitionSize == partitionSizes[ip]) {
							break;
						}
					}
				}
			}
			ic = ic_next;
			if ((ic == -1) && (partitionSize < partitionSizes[ip])) {
				for (int i = 0; i < nCells; ++i) {
					if (cellPartitions[i] == -1) {
						ic = i;
						cellPartitions[ic] = ip;
						++partitionSize;
						break;
					}
				}
			}
		}
	}

//	C.resize(nPartitions, std::vector<std::vector<int>> (nColors));

	for (int partition = 0; partition < nPartitions; ++partition) {
		std::vector <std::vector <int>> v;
		for (int color = 0; color < nColors; ++color) {
			v.push_back(std::vector<int>());
		}
		C.push_back(v);
	}

	for (int ic = 0; ic < nCells; ++ic) {
		int partition = cellPartitions[ic];
		int color = cellColors[ic];
		C[partition][color].push_back(ic);
	}

}

void Mesh::write_tecplot(std::vector < std::vector <double> > data, std::string filename,
		std::vector <std::string> var_names, double time) {

	int nvar = data[0].size();
	std::ofstream file;
	file.open(filename);
	if (file.fail())
	{
		std::cout << "Could not open " << filename << std::endl;
		exit(1);
	}
	file << "TITLE = \"VolumeData\"\n";
	file << "VARIABLES = \"x\" \"y\" \"z\" ";
	for (int iv = 0; iv < nvar; ++iv) {
		file << " \"" << var_names[iv] << "\" ";
	}
	file << "\n";
	file << "ZONE T=\"my_zone\", SolutionTime=" << time <<
			", DATAPACKING=Block, ZONETYPE=FEBRICK, Nodes=" << nVerts <<
		   ", Elements=" << nCells;
	file << ", VarLocation=([4-" << 3+nvar << "]=CellCentered)";
	// Write vertices' coo;
	for (int i = 0; i < 3; ++i) {
		for (int iv = 0; iv < nVerts; ++iv) {
			file << vertsCoo[iv][i] << "\n"; // TODO: format
		}
	}
	// Write values of variables
	for (int i = 0; i < nvar; ++i) {
		for (int ic = 0; ic < nCells; ++ic) {
			file << data[ic][i] << "\n"; // TODO: format
		}
	}
	// Write cell-to-vertices connectivity
	for (int ic = 0; ic < nCells; ++ic) {
		std::vector <int> verts = cellVerts[ic]; // TODO: format
		if (cellTypes[ic] == HEXA) {
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
		else if (cellTypes[ic] == TETRA) {
			file << verts[1] + 1 << " ";
			file << verts[2] + 1 << " ";
			file << verts[3] + 1 << " ";
			file << verts[3] + 1 << " ";
			file << verts[0] + 1 << " ";
			file << verts[0] + 1 << " ";
			file << verts[0] + 1 << " ";
			file << verts[0] + 1 << " ";
			file << "\n";
		}
	}
}
