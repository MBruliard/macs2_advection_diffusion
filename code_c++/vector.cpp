/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/


#include "vector.hpp"


/**------------------------CONSTRUCTEURS VECTOR----------------------------**/
/*
Vecteur::Vecteur (const Vecteur &u)
{

	_nbLignes = u._nbLignes;
	_nbColonnes = u._nbColonnes;

	for (int i=0; i<_nbLignes; i++)
	{
		_values[i] = u._values[i];
	}
}*/


/**------------------ OPERATEURS DE VECTEUR ------------------------------**/

double Vecteur::operator () (int i)
{
	return _values[i];
}

Vecteur Vecteur::operator - ()
{
	Vecteur res(_nbLignes);
	
	for (int i=0; i<_nbLignes; i++)
	{
		res._values[i] = - _values[i];
	}
	return res;
}



/**---------------------------- FONCTIONS ASSOCIEES A VECTEUR-------------------------**/
Vecteur sommeVector (const Vecteur &u, const Vecteur &v)
{
	/*
	Somme de deux vecteurs de même taille. Si ce n'est pas possible alors on renvoie un vecteur nul de la taille de @u et on affiche une erreur
	*/
	
	Vecteur res(u._nbLignes);
	
	if (u._nbLignes != v._nbLignes)
	{
		std::cout << "***ERREUR: Les vecteurs ne sont pas de meme dimensions dans operator + (const Vecteur &u, const Vecteur &v)***" << std::endl;
		return res;
	}
	
	for (int i=0; i<u._nbLignes; i++)
	{
		res._values[i] = u._values[i] + v._values[i];
	}
	
	return res;
}



Vecteur produitMatVect (const Matrice &u, const Vecteur &v)
{
	/*
	* Produit entre une matrice et un vecteur. 
	* /!\ il faut prendre en compte les dimensions. -> si oui on fait le calcul et ->sinon alors on renvoie un vecteur nul de la taille de la matrice et on affiche un message d'erreur
	*/
	
	Vecteur res(v._nbLignes);
	if (u._nbColonnes != v._nbLignes)
	{
		std::cout << "***ERREUR: Les dimensions sont conflictuelles dans operator * (const Matrice &u, const Vecteur &v)***" << std::endl;
		return res;
	}

	for (int i=0; i<u._nbLignes; i++)
	{
		double somme = 0;
		for (int k=0; k<u._nbColonnes; k++)
		{
			somme = somme + u._values[i + k*u._nbLignes]*v._values[k];
		}
		res._values[i] = somme;
	}

	return res;
}



double produitScalaire (const Vecteur &u, const Vecteur &v)
{
	if (u._nbLignes == v._nbLignes)
	{
		double somme = 0; 
	
		for (int i=0; i<u._nbLignes; i++ )
		{
			somme = somme + u._values[i]*v._values[i];
		}
		return somme;
	}
	else
	{
		std::cout << "***ERREUR: Les vecteurs n'ont pas les bonnes dimensions dans operator * (const Vecteur &u, const Vecteur &v)***" << std::endl;
		return 0.0;
	}
		
}







