#include "header.h"
#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"

//#include "mesh2.h"

using namespace MEDCoupling;

int main(int argc, char *argv[])
{
	MEDFileUMesh *mesh = MEDFileUMesh::New("/home/egor/git/boltz-t-cpp/mesh-cyl/mesh.med", "mesh-cyl");

	bool b = mesh->existsGroup("SYMMETRYZ");

	std::cout << b << std::endl;

	return 0;
}
