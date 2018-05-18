/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/


#ifndef __VECTOR_HPP__
#define __VECTOR_HPP__


#include "matrix.hpp"

/*----------------------- CLASSE VECTEUR --------------------------*/
class Vecteur: public Matrice
{
	/*
	Par définition, on pose les vecteurs comme des matrices de 1 colonne.
	*/
	
	
	
	public:
		Vecteur();
		Vecteur (int taille): Matrice(taille, 1) {}
		Vecteur (int taille, std::vector<double> initialValues): Matrice(taille, 1, initialValues) {}
		Vecteur (const Vecteur &u): Matrice(u) {}
		
		
		/*methodes*/
		
		
		/*operateurs*/
		//Vecteur operator = (const Vecteur &u);
		double operator () (int i);
		double operator () (int i, int j);
		Vecteur operator - ();
		
};

Vecteur sommeVector (const Vecteur &u, const Vecteur &v);
Vecteur produitMatVect (const Matrice &u, const Vecteur &v);
double produitScalaire (const Vecteur &u, const Vecteur &v);

inline double operator * (const Vecteur &u, const Vecteur &v) { return produitScalaire(u, v);}
inline Vecteur operator + (const Vecteur &u, const Vecteur &v) { return sommeVector(u, v) ;}
inline Vecteur operator * (const Matrice &u, const Vecteur &v) { return produitMatVect (u, v) ; }
inline Vecteur operator * (const double alpha, const Vecteur &u) 
{
	Vecteur res(u._nbLignes);
	for (int i=0; i<u._nbLignes; i++)
	{
		res._values[i] = alpha* u._values[i];
	}
	return res;
}


#endif //__VECTOR_HPP__
