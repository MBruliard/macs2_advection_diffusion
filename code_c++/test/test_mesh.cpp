/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/

#include "../mesh.hpp"
#include <iostream>

int main ()
{
	MeshRegulier lectMesh("data/mesh/mesh1.dat");
	
	std::cout << lectMesh._nx << " " << lectMesh._ny << std::endl;
	std::cout << lectMesh._discretX[lectMesh._nx-2] << std::endl;
	std::cout << lectMesh._discretX[lectMesh._nx-1] << std::endl;
	lectMesh.genererGrapheGnuplot("plot/data_mesh1.dat", "plot/maillage_mesh1.gnu");
	lectMesh.lancerGraphe("plot/maillage_mesh1.gnu");
	return 0;
}
