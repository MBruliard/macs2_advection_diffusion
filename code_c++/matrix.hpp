/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/

#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__


#include <vector>
#include <iostream>
#include <fstream>
#include <string>

//#include "sparse_matrix.hpp"

/*-------------------- CLASSE MATRICE ------------------------*/
class Matrice
{
	/*
	Cette classe définie les matrices pleines pour les petits tests.
	On définie A(i, j) comme suit A(i, j) = _values[j + i*_nbLignes] -> donc on stocke la matrice colonne par colonne
	*/
	private:
	
	public:
		/*ATTRIBUTS*/
		int _nbLignes;
		int _nbColonnes;
		std::vector<double> _values;
		
		
		/*constructeurs*/
		Matrice();
		Matrice(int lignes, int colonnes);
		Matrice(int lignes, int colonnes, std::vector<double> values);
		Matrice (const std::string chemin);
		
		Matrice (const Matrice &u);
		
		
		/*methodes*/
		std::ostream& affichageMatrix (std::ostream& flux) const;
		
		//transformer une matrice en sparse
		//SparseMatrice toSparse();
		
		
		/*operateurs*/
		Matrice operator = (const Matrice &u);
		Matrice operator - (const Matrice &u);
		double operator () (int i, int j);
		double operator () (int i);
		

};

Matrice sommeMatrix (const Matrice &u, const Matrice &v);
Matrice multiplierMatrix (const Matrice &u, const Matrice &v);

inline Matrice operator + (const Matrice &u, const Matrice &v) { return sommeMatrix(u, v); }
inline Matrice operator * (const Matrice &u, const Matrice &v) { return multiplierMatrix(u, v); }
inline std::ostream& operator << (std::ostream& out, const Matrice& mat) { return mat.affichageMatrix(out); }


#endif //__MATRIX_HPP__
