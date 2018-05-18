/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/


#include "sparse_matrix.hpp"



/**-------------------------------------CONSTRUCTEURS DE LA CLASSE SPARSE MATRICE --------------------------------**/
SparseMatrice::SparseMatrice (int lignes, int colonnes)
{
	/*
	constructeur de la classe  SparseMatrice. On initialise la matrice à 0 pour chaque élément. Donc dans le cas d'une sparseMatrix on ne remplit pas les std::vector _values et _ind_values
	*/
	
	_nbLignes = lignes;
	_nbColonnes = colonnes;
	
}


SparseMatrice::SparseMatrice (int lignes, int colonnes, std::vector<double> values, std::vector<int> indices)
{
	/*
	constructeur de la classe. On initialise les valeurs contenues dans values.
	*/
	_nbLignes = lignes;
	_nbColonnes = colonnes;
	
	_values = values;
	_ind_values = indices;
}

SparseMatrice::SparseMatrice (const SparseMatrice &u)
{
	_nbLignes = u._nbLignes;
	_nbColonnes= u._nbColonnes;
	
	for (unsigned int i=0; i<u._values.size(); i++)
	{
		_values.push_back(u._values[i]);
		_ind_values.push_back(u._ind_values[i]);
	}
}


/**--------------------------------------- METHODES DE LA CLASSE SPARSE MATRICE ------------------------------------------**/

int SparseMatrice::indInIndValues (int ind)
{
	/*
	on vérifie si l'indice @ind donné est dans la liste des indices non nuls de la matrice. 
	Si vrai, on renvoie la place de ind dans _ind_values
	sinon on renvoie -1
	*/
	for (unsigned int i=0; i< _ind_values.size(); i++)
	{
		if (ind == _ind_values[i])
		{
			return i;
		}
	}
	return -1;
}


std::ostream& SparseMatrice::affichageMatrix (std::ostream& flux)
{
	/*
	permet l'affichage d'une sparse matrice dans un flux en écriture
	@k est l'indice permettant de parcourir _ind_values et _values
	@k_max est l'indice maximal de k
	*/
	
	for (int i=0; i < _nbLignes; i++)
	{
		for (int j=0; j< _nbColonnes; j++)
		{
			int ind_in_ind_values = indInIndValues(i+j*_nbLignes);
			if (ind_in_ind_values >= 0 )
			{
				flux << _values[ind_in_ind_values] << "\t";
			}
			else
			{
				flux << 0.0 << "\t";
			}
		}
		flux << "\n";
	}
	return flux;
}


/**-------------------------------------------- OPERATEURS DE LA CLASSE SPARSE MATRICE ----------------------------**/



SparseMatrice SparseMatrice::operator = (const SparseMatrice &u)
{
	/*
	operateur egalité pour une matrice.
	/!\ Il faut adapter la taille de _values en fonction de u._values
	*/
	
	if (this != &u)
	{
				
		//il nous faut autant de place dans le _values que dans u._values donc il faut vérifier le nb de cases dans chacun
		unsigned int sizeUvalues = u._values.size();
		unsigned int sizeValues = this->_values.size();
		unsigned int i=0;
		
		if (sizeUvalues <= sizeValues)
		{
			/*
			*	il y a deja l'espace necessaire donc on recopie les valeurs tout simplement et on supprimera ensuite les cases en surplus de vector
			*/
			
			for (i=0; i<sizeUvalues; i++)
			{
				this->_values[i]  = u._values[i];
				this->_ind_values[i] = u._ind_values[i];
			}
			
			unsigned int j= sizeValues;
			while (j>i)
			{
				this->_values.pop_back();
				this->_ind_values.pop_back();
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
				this->_ind_values[i] = u._ind_values[i];
			}
			while (i < sizeUvalues)
			{
				this->_values.push_back(u._values[i]);
				this->_ind_values.push_back(u._ind_values[i]);
				i = i+1;
			}
			
		}
		
		this->_nbColonnes = u._nbColonnes;
		this->_nbLignes = u._nbLignes;
	}
	return *this;
}


SparseMatrice SparseMatrice::operator - (const SparseMatrice &u)
{
	/*
	* operateur - pour rendre négatif une matrice
	*/
	SparseMatrice res(u);
	for (unsigned int i=0; i<u._values.size(); i++)
	{
		res._values[i] = - u._values[i];
	}
	
	return res;	
}

double SparseMatrice::operator () (int i, int j)
{
	/*
	Operateur permettant d'appeler via ces coordonnées (i,j) la sparse Matrice M 
	*/
	
	//on cherche si l'indice demandé appartient à _ind_values
	int k = this->indInIndValues(i+j*_nbLignes);
	if (k >=0)
	{
		return _values[k];
	}
	return 0;
}

double SparseMatrice::operator() (int i)
{
	/*
	Operateur permettant d'appeler via son indice i la valeur de la sparse matrice M
	On a : i = i mod _nbLignes + (i - (i mod _nbLignes));
	*/
	unsigned int k=0;
	while (k < _ind_values.size())
	{
		if (i == _ind_values[k])
		{
			return _values[k];
		}
		k=k+1;
	}
	return 0. ;
	
}


/**---------------------------------FONCTIONS ASSOCIEES A LA CLASSE MATRICE --------------------------------**/
SparseMatrice sommeSparseMatrix (const SparseMatrice &u, const SparseMatrice &v)
{
	/*
	Méthode permettant de faire la somme de deux matrices de meme tailles. Si ce n'est pas le cas, elle renvoie une matrice nulle de taille u
	*/
	
	SparseMatrice res(u._nbLignes, u._nbColonnes); 
	
	if (!((v._nbColonnes == u._nbColonnes) && (v._nbLignes == u._nbLignes)))
	{
		std::cout << "***ERREUR: Les matrices ne sont pas de meme dimensions dans operator + (const SparseMatrice &u)***" << std::endl;
		std::exit(-1);
	}
	
	//on va inserer les elements de B manquants
	unsigned int l=0; //parcours _ind_values de v
	unsigned int l_max = v._ind_values.size();
	
	for (unsigned int k=0; k< u._ind_values.size(); k++)
	{
		while((u._ind_values[k] < v._ind_values[l]) && (l < l_max))
		{
			//on est dans le cas ou A(_ind_values[k]) = 0 et B(_ind_values[l] !=0) -> res(_ind_values[l]) = B(_ind_values[l] !=0)
			res._ind_values.push_back(v._ind_values[l]);
			res._values.push_back(v._values[l]);
			l = l+1;
		}
		
		if ((u._ind_values[k] == v._ind_values[l]) && (l<l_max))
		{
			res._ind_values.push_back(v._ind_values[l]);
			res._values.push_back(u._values[k] + v._values[l]);
			l=l+1;
		}
	}
	
	
	return res;
}

/*
SparseMatrice multiplierSparseMatrix (const SparseMatrice &u, const SparseMatrice &v)
{
	*//*
	/!\ NE FONCTIONNE PAS !!!!!	
	*/
	/*
	Méthode permettant de faire le produit de deux matrices. Si les dimensions ne correspondent pas, elle renvoie une matrice nulle de taille u. 
	Dim de u = (ul, uc) et dim de v= (vl, vc). On doit tester que uc et vl sont identiques
	*/
	/*SparseMatrice res(u._nbLignes, v._nbColonnes); 
	
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
			if (somme != 0.0)
			{
				res._ind_values.push_back(j+ i*u._nbLignes);
				res._values.push_back(somme);
			}
		}
	}
	
	return res;
}
*/

























