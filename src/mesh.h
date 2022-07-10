#ifndef MESH_H_
#define MESH_H_

#include "header.h"

enum bcType
{
    SYMMETRYZ = 0,
    INLET = 1,
    OUTLET = 2,
    WALL = 3,
    SYMMETRYY = 4,
    SYMMETRYX = 5 
};

enum cellType
{
	TETRA = 0,
	HEXA = 1
};

enum faceType
{
	TRIANGLE = 0,
	QUAD = 1
};

class Mesh {
public:
	// Constructor
	Mesh(const std::string& path, double scale = 1.0);
	void read(const std::string& path, double scale = 1.0);

	int nVerts;
	std::vector<std::vector<double>> vertsCoo;
	void readVerts(std::ifstream &data, double scale);
	std::vector<std::vector<int>> vertsCells;

	int nCells;
	std::vector<std::vector<int>> cellVerts;
	std::vector<cellType> cellTypes;
	void readHexa(std::ifstream &data);
	void readTetra(std::ifstream &data);
	std::vector<std::vector<double>> cellCenters;
	std::vector<double> cellVolumes;
	std::vector<std::vector<int>> computeFacesOfCell(int ic);
	std::vector<std::vector<int>> cellNeighbors;
	std::vector<std::vector<int>> cellFaces;
	std::vector<std::vector<bool>> isOuterNormal;
	int getOutIndex(int ic, int j); // 1 if out, 0 else
	double getOutSign(int ic, int j); // 1.0 if out, -1.0 else
	std::vector<double> cellDiameters;
	std::vector<int> cellColors;
	std::vector<std::vector<int>> coloredCells;
	std::vector<std::vector<int>> rcoloredCells;
	int nColors;

	int nBoundaryFaces;
	std::vector<std::vector<int>> boundaryFaces;
	std::vector<faceType> boundaryFaceTypes;
	std::vector<int> boundaryFaceTags;
	std::vector<int> boundaryFaceTagsSet;
	void readBoundaryQuads(std::ifstream &data);
	void readBoundaryTriangles(std::ifstream &data);
	std::vector<int> boundaryFaceIndex;
	std::map<int, std::vector<int>> boundaryFacesForEachTag; // TODO or unordered?
	std::vector<bool> isOuterNormalBoundary;
	int getOutIndex(int jf); // 1 if out, 0 else
	double getOutSign(int jf); // 1.0 if out, -1.0 else

	int nFaces;
	std::vector<std::vector<int>> faces;
	std::vector<faceType> faceTypes;
	std::vector<std::vector<int>> leftRightCells;
	std::vector<double> faceAreas;
	std::vector<std::vector<double>> faceNormals;
	std::vector<std::vector<double>> faceCenters;

	double computeTetraVolume(std::vector<std::vector<double>> tetra);


	void write_tecplot(std::vector < std::vector <double> > data, std::string filename,
			std::vector <std::string> var_names, double time = 0.0);

	int getNumFacesOfCell(int ic);
	int getFacesOfCell(int ic, int jf);
	int getNumTags();
	std::vector<int> getFacesForTag(int tag);
};

#endif /* MESH_H_ */
