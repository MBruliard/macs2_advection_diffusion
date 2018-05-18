/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/

#include "mesh.hpp"

/**--------------------------------------- CONSTRUCTEURS -------------------------------------**/
MeshCarre::MeshCarre (double a, double b, int nx, int ny)
{
	_xmax = a; 
	_ymax = b;
	_nx= nx;
	_ny= ny;
	
	/*on fait la discrétisation*/
	/*partie selon _nx*/
	for (int i=0; i<nx; i++)
	{
		_discretX.push_back(i*a/(nx-1));
	}
	/*partie selon ny*/
	for (int j=0; j<ny; j++)
	{
		_discretY.push_back(j*b/(ny -1));
	}
	
}


MeshCarre::MeshCarre (int nx, int ny)
{
	/*MeshCarre pour un domaine [0,1]*[0,1]*/
	_xmax = 1;
	_ymax = 1;
	_nx = nx;
	_ny = ny;
	
	/*partie selon _nx*/
	for (int i=0; i<nx; i++)
	{
		_discretX.push_back(i*_xmax/(nx-1));
	}
	/*partie selon ny*/
	for (int j=0; j<ny; j++)
	{
		_discretY.push_back(j*_ymax/(ny -1));
	}
	
}


MeshCarre::MeshCarre (const MeshCarre &u)
{
	_xmax = u._xmax;
	_ymax = u._ymax;
	_nx = u._nx;
	_ny = u._ny; 
	
	/*partie selon _nx*/
	for (int i=0; i<_nx; i++)
	{
		_discretX.push_back(u._discretX[i]);
	}
	/*partie selon ny*/
	for (int j=0; j<_ny; j++)
	{
		_discretY.push_back(u._discretY[j]);
	}
}

MeshCarre::MeshCarre(std::string chemin)
{
	std::ifstream fichier (chemin, std::ios::in);
	
	if (fichier)
	{
		double a, b;
		int nx, ny;
		
		fichier >> a;
		fichier >> b;
		fichier >> nx;
		fichier >> ny;
		
		MeshCarre mesh (a, b, nx, ny);
		
		fichier.close();
	}
	
	
}



/**---------------------------------------- ACTIONS SUR LES FICHIERS -------------------------------------**/
void MeshCarre::save(std::string chemin)
{
	/*
	on va enregistrer les infos dans l'ordre suivant: a, b \\ nx, ny 
	*/
	
	std::ofstream fichier(chemin, std::ios::out);
	
	if (fichier)
	{
		fichier << _xmax << "\n" << _ymax << "\n";
		fichier << _nx << "\n" << _ny << "\n";
		
		fichier.close();
	}
}



















