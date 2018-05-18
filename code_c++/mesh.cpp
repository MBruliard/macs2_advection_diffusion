/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/

#include "mesh.hpp"

/**--------------------------------------- CONSTRUCTEURS -------------------------------------**/
MeshRegulier::MeshRegulier (double a, double b, int nx, int ny)
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


MeshRegulier::MeshRegulier (int nx, int ny)
{
	/*MeshRegulier pour un domaine [0,1]*[0,1]*/
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


MeshRegulier::MeshRegulier (const MeshRegulier &u)
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

MeshRegulier::MeshRegulier (std::string chemin)
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
		
		_xmax = a;
		_ymax =b;
		_nx = nx;
		_ny = ny;
		
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
		
		
		fichier.close();
	}
	
	
}



/**---------------------------------------- ACTIONS SUR LES FICHIERS -------------------------------------**/
void MeshRegulier::save(std::string chemin)
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

void MeshRegulier::genererGrapheGnuplot (std::string chemindonnees, std::string chemingnuplot)
{
	std::ofstream fichier_donnees(chemindonnees, std::ios::out);
	std::ofstream fichier_gnuplot(chemingnuplot, std::ios::out);
	
	std::string nomdonnees=chemindonnees.substr(5); //il s'agit du nom de fichier des données
	
	
	std::string chemingraphe = chemingnuplot.substr(0, chemingnuplot.rfind(".gnu")) + ".svg";
	std::cout << chemingraphe << std::endl;
	
	
	if (fichier_donnees)
	{
		for (int i=0; i<_nx; i++)
		{
			for (int j=0; j<_ny; j++)
			{
				fichier_donnees << _discretX[i] << "\t" << _discretY[j] << "\n";
			}
		}
		
		std::cout << "donnees enregistrees dans : " << chemindonnees << std::endl;
		fichier_donnees.close();
	}
	else
	{
		std::cout << "** ERREUR: Le fichier " << chemindonnees << " ne s'est pas ouvert correctement...\nArret du programme **" << std::endl;
		std::exit(-1);
	}
	
	
	if (fichier_gnuplot)
	{
		
		fichier_gnuplot << "set terminal svg enhanced background rgb 'white' size 1000,800\n";
		fichier_gnuplot << "set output '" << chemingraphe << "'\n";
		fichier_gnuplot << "set title 'Maillage' \n";
		fichier_gnuplot << "set time\n" ;
		fichier_gnuplot << "set xlabel 'discretisation en x'\n";
		fichier_gnuplot << "set ylabel 'discretisation en y'\n";
		fichier_gnuplot << "plot '" << chemindonnees << "' with points\n";
		//fichier_gnuplot << "pause -1";
		
		
		
		std::cout << "graphe enregistre dans : " << chemingnuplot << std::endl;
		fichier_gnuplot.close();
	}
	else
	{
		std::cout << "** ERREUR: Le fichier " << chemindonnees << " ne s'est pas ouvert correctement...\nArret du programme **" << std::endl;
		std::exit(-1);
	}
}


void MeshRegulier::lancerGraphe (std::string chemingraphe)
{
	
	std::string str = "gnuplot " + chemingraphe;
	
	std::system(str.c_str());
	
}

















