#include "mesh2.h"
#include "header.h"

int main() {

	Mesh mesh("/home/egor/git/boltz-t-cpp/mesh-1d/");

	std::cout << mesh.nBoundaryFaces << std::endl;
	std::cout << mesh.nCells << std::endl;
	std::cout << mesh.nVerts << std::endl;
	std::cout << mesh.nFaces << std::endl;


	std::cout << "boundaryFaces" << std::endl;
	for (const auto &bcface : mesh.boundaryFaces) {
		for (const auto &v : bcface) {
			std::cout << v << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "boundaryFaceIndex" << std::endl;
	for (const auto &index : mesh.boundaryFaceIndex) {
			std::cout << index << "\n";
	}
	std::cout << "\n";

	std::cout << "boundaryFaceTypes" << std::endl;
	for (const auto &type : mesh.boundaryFaceTypes) {
			std::cout << type << " ";
	}
	std::cout << "\n";

	std::cout << "boundaryFaceTags" << std::endl;
	for (const auto &type : mesh.boundaryFaceTags) {
		std::cout << type << " ";
	}
	std::cout << "\n";

	std::cout << "cellVerts" << std::endl;
	for (const auto &i : mesh.cellVerts) {
		for (const auto &j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "boundaryFaceForEachTag" << std::endl;
	for (const auto &type : mesh.boundaryFaceTagsSet) {
		std::cout << "type: " << type << "\n";
		for (const auto &j : mesh.boundaryFacesForEachTag[type]) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "vertsCoo" << std::endl;
	for (const auto &i : mesh.vertsCoo) {
		for (const auto &j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cellCentersCoo" << std::endl;
	for (const auto &i : mesh.cellCenters) {
		for (const auto &j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cellVolumes" << std::endl;
	for (const auto &a : mesh.cellVolumes) {
		std::cout << a << " ";
	}
	std::cout << "\n";

	std::cout << "vertsCells" << std::endl;
	for (const auto &i : mesh.vertsCells) {
		for (const auto &j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cell_num_for_vertex" << std::endl;
	for (const auto &a : mesh.vertsCells) {
		std::cout << a.size() << " ";
	}
	std::cout << "\n";

	std::cout << "cellNeighbors" << std::endl;
	for (const auto &i : mesh.cellNeighbors) {
		for (const auto &j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "faces" << std::endl;
	for (const auto &i : mesh.faces) {
		for (const auto &j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cellFaces" << std::endl;
	for (const auto &i : mesh.cellFaces) {
		for (const auto &j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "faceAreas" << std::endl;
	for (const auto &a : mesh.faceAreas) {
		std::cout << a << " ";
	}
	std::cout << "\n";

	std::cout << "faceNormals" << std::endl;
	for (const auto &i : mesh.faceNormals) {
		for (const auto &j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "isOuterNormal" << std::endl;
	for (const auto &i : mesh.isOuterNormal) {
		for (const auto &j : i) {
			std::cout << (j) << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "bound_face_info" << std::endl;
	for (int ibf = 0; ibf < mesh.nBoundaryFaces; ++ibf) {
		std::cout << mesh.boundaryFaceIndex[ibf] << " ";
		std::cout << mesh.boundaryFaceTags[ibf] << " ";
		int jf = mesh.boundaryFaceIndex[ibf];
		std::cout << mesh.getOutIndex(jf) << " ";
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cellDiameters" << std::endl;
	for (const auto &a : mesh.cellDiameters) {
		std::cout << a << " ";
	}
	std::cout << "\n";

	std::cout << "faceCenters" << std::endl;
	for (const auto &i : mesh.faceCenters) {
		for (const auto &j : i) {
			std::cout << j << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "cellTypes" << std::endl;
	for (const auto &a : mesh.cellTypes) {
		std::cout << a << " ";
	}
	std::cout << "\n";

	return 0;
}
