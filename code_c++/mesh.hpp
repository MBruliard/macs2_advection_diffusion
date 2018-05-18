/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/

#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <vector>

/*
Classe permrettant de définir un mesh regulier.
soit les parametres suivants:
@_nx est le nombre de discretisations en x
@_ny le nombre de discretisations en y
la distance @a est la longueur du cote selon x
la distance @b est la longueur du cote selon y
@
*/

#include <vector>
#include <string>
#include <fstream>

/*
@_xmax et @_ymax
@_nx et @_ny sont le nombre de discrétisations entre 0 et _xmax et 0 et _ymax
@_discretX sont les différents x discrétisées
@_discretY sont les différents y discrétisées
*/



class MeshCarre
{
	public:
		
		/*attributs*/
		double _xmax, _ymax;
		int _nx, _ny;
		std::vector<double> _discretX; 
		std::vector<double> _discretY;
		
		
		/*constructeurs*/
		MeshCarre () = default;
		
		MeshCarre (double a, double b, int nx, int ny);
		MeshCarre (std::string chemin);
		MeshCarre (int nx, int ny);
		MeshCarre (const MeshCarre &u);
		
		/*lecture et ecriture dans les fichiers .dat*/
		void save(std::string chemin);
		
		
		


};



#endif //__MESH_HPP__
