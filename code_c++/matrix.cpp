/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/


#include "matrix.hpp"



/**-------------------------------------CONSTRUCTEURS DE LA CLASSE MATRICE --------------------------------**/
Matrice::Matrice (int lignes, int colonnes)
{
	/*
	constructeur de la classe Matrice. On initialise la matrice à 0 pour chaque élément
	*/
	
	_nbLignes = lignes;
	_nbColonnes = colonnes;
	
	for (int i=0; i<lignes; i++)
	{
		for (int j=0; j<colonnes; j++)
		{
			_values.push_back(0.0);
		}
	}
}


Matrice::Matrice (int lignes, int colonnes, std::vector<double> values)
{
	/*
	constructeur de la classe. On initialise les valeurs contenues dans values pour le tableau.
	@values : les valeurs sont bien rangé dans l'ordre A(i, j) = values[i + j*_nbColonnes]
	*/
	_nbLignes = lignes;
	_nbColonnes = colonnes;
	
	_values = values;
}

Matrice::Matrice (const Matrice &u)
{
	_nbLignes = u._nbLignes;
	_nbColonnes= u._nbColonnes;
	
	for (int i=0; i<_nbLignes; i++)
	{
		for (int j=0; j<_nbColonnes; j++)
		{
			_values.push_back(u._values[j+i*u._nbColonnes]);
		}
	}
}


Matrice::Matrice (const std::string chemin)
{
	 
	std::ifstream fichier(chemin);  // on ouvre le fichier en lecture

	int l, c;
	double tampon;
	

	if(fichier)  
	{
		fichier >> l;
		fichier >> c;
		
		_nbLignes = l;
		_nbColonnes = c;
		
		for (int i=0; i<l*c; i++)
		{
			fichier >> tampon;
			std::cout << tampon << " ";
			_values.push_back(tampon);
		}
		std::cout << std::endl;
		
		
		fichier.close();  // on ferme le fichier

	}
	else
	{
		std::cout << "ERREUR: Le fichier " << chemin << " ne s'est pas ouvert correctement" << std::endl;
	}

}



/**--------------------------------------- METHODES DE LA CLASSE MATRICE ------------------------------------------**/

std::ostream& Matrice::affichageMatrix (std::ostream& flux) const
{
	/*
	permet l'affichage d'une matrice dans un flux en écriture
	*/
	for (int i=0; i < _nbLignes; i++)
	{
		for (int j=0; j< _nbColonnes; j++)
		{
			flux << _values[i + j*_nbLignes]<<"\t";
			
		}
		flux << "\n";
	}
	return flux;
}





/**-------------------------------------------- OPERATEURS DE LA CLASSE MATRICE ----------------------------**/



Matrice Matrice::operator = (const Matrice &u)
{
	/*
	operateur egalité pour une matrice.
	/!\ Il faut adapter la taille de _values en fonction de u._values
	*/
	if (this != &u)
	{
				
		//il nous faut autant de place dans le _values que dans u._values donc il faut vérifier le nb de cases dans chacun
		int sizeUvalues = u._values.size();
		int sizeValues = this->_values.size();
		int i=0;
		
		if (sizeUvalues <= sizeValues)
		{
			/*
			*	il y a deja l'espace necessaire donc on recopie les valeurs tout simplement et on supprimera ensuite les cases en surplus de vector
			*/
			for (i=0; i<sizeUvalues; i++)
			{
				this->_values[i]  = u._values[i];
			}
			int j= sizeValues;
			while (j>i)
			{
				this->_values.pop_back();
				j = j-1;
			}
		}
		else
		{
			/*
			*	Dans ca cas on copie ce qu'on peut dans _values et on va rajouter les cases manquantes.
			*/ 
			for (i=0; i<sizeValues; i++)
			{
				this->_values[i] = u._values[i];
			}
			while (i < sizeUvalues)
			{
				this->_values.push_back(u._values[i]);
				i = i+1;
			}
			
		}
		
		this->_nbColonnes = u._nbColonnes;
		this->_nbLignes = u._nbLignes;
	}
	return *this;
}

Matrice Matrice::operator - (const Matrice &u)
{
	/*
	* operateur - pour rendre négatif une matrice
	*/
	Matrice res (u._nbLignes, u._nbColonnes);
	for (unsigned int i=0; i<u._values.size(); i++)
	{
		res._values[i] = - u._values[i];
	}
	
	return res;	
}

double Matrice::operator () (int i, int j)
{
	/*
	Operateur permettant d'appeler via ces coordonnées (i,j) la Matrice M 
	*/
	return _values[i + j*_nbLignes] ;
}

double Matrice::operator() (int i)
{
	/*
	Operateur permettant d'appeler via son indice i la valeur de la matrice M
	On a : i = i mod _nbLignes + (i - (i mod _nbLignes));
	*/
	
	return _values[i];
}

/**---------------------------------FONCTIONS ASSOCIEES A LA CLASSE MATRICE --------------------------------**/
Matrice sommeMatrix (const Matrice &u, const Matrice &v)
{
	/*
	Méthode permettant de faire la somme de deux matrices de meme tailles. Si ce n'est pas le cas, elle renvoie une matrice nulle de taille u
	*/
	
	Matrice res(u._nbLignes, u._nbColonnes); 
	
	if (!((v._nbColonnes == u._nbColonnes) && (v._nbLignes == u._nbLignes)))
	{
		std::cout << "***ERREUR: Les matrices ne sont pas de meme dimensions dans operator + (const Matrice &u)***" << std::endl;
		return res;
	}
	
	for (unsigned int i =0; i<res._values.size(); i++)
	{
		res._values[i] = v._values[i] + u._values[i];
	}
	
	return res;
}


Matrice multiplierMatrix (const Matrice &u, const Matrice &v)
{
	/*
	Méthode permettant de faire le produit de deux matrices. Si les dimensions ne correspondent pas, elle renvoie une matrice nulle de taille u. 
	Dim de u = (ul, uc) et dim de v= (vl, vc). On doit tester que uc et vl sont identiques
	*/
	
	Matrice res(u._nbLignes, v._nbColonnes); 
	
	if (!(u._nbColonnes == v._nbLignes ))
	{
		std::cout << "***ERREUR: Les matrices n'ont pas les bonnes dimensions dans operator * (const Matrice &u, const Matrice &v)***" << std::endl;
		return res;
	}
	
	for (int i =0; i<u._nbLignes; i++)
	{
		for (int j=0; j<v._nbColonnes; j++)
		{
			double somme = 0.0;
			for (int k=0; k<u._nbColonnes; k++)
			{
				somme = somme + u._values[i+k*u._nbColonnes]*v._values[k+j*v._nbLignes];
			}
			res._values[j+ i*u._nbLignes]= somme;
		}
	}
	
	return res;
}








