/*
*	BRULIARD - RIGAL
*	Projet Numérique - MACS 2
*/

#ifndef __SPARSE_MATRIX_HPP__
#define __SPARSE_MATRIX_HPP__


#include <vector>
#include <iostream>

/*-------------------- CLASSE MATRICE ------------------------*/
class SparseMatrice
{
	/*
	Cette classe définie les matrices "SPARSE" pour les calculs.
	@_nbLignes et @_nbColonnes définissent la taille de la matrice
	@_values est un tableau contennant les valeurs non nulles de la matrice
	@_ind_values est le tableau contennant les indices des valeurs non nulles de la matrice: on rappelle que les elements sont stockés colonne par colonne: i = i mod _nbLignes + (i - (i mod _nbLignes));
	
	Ainsi supposons la matrice suivante:
	1	0
	2	1
	
	on obtiendra _values = {1,2,1} et _ind_values = {0, 1, 3}	
	*/
	private:
	
	public:
		/*ATTRIBUTS*/
		int _nbLignes;
		int _nbColonnes;
		std::vector<double> _values;
		std::vector<int> _ind_values;
		
		
		/*constructeurs*/
		SparseMatrice();
		SparseMatrice(int lignes, int colonnes);
		SparseMatrice(int lignes, int colonnes, std::vector<double> values, std::vector<int> indices);
		
		SparseMatrice (const SparseMatrice &u);
		
		
		/*methodes*/
		int indInIndValues (int ind);
		std::ostream& affichageMatrix (std::ostream& flux);
		
		
		/*operateurs*/
		SparseMatrice operator = (const SparseMatrice &u);
		SparseMatrice operator - (const SparseMatrice &u);
		double operator () (int i, int j);
		double operator () (int i);
		

};

SparseMatrice sommeSparseMatrix (const SparseMatrice &u, const SparseMatrice &v);
//SparseMatrice multiplierSparseMatrix (const SparseMatrice &u, const SparseMatrice &v);

inline SparseMatrice operator + (const SparseMatrice &u, const SparseMatrice &v) { return sommeSparseMatrix(u, v); }
//inline SparseMatrice operator * (const SparseMatrice &u, const SparseMatrice &v) { return multiplierSparseMatrix(u, v); }
inline std::ostream& operator << (std::ostream& out, SparseMatrice& mat) { return mat.affichageMatrix(out); }











#endif //__SPARSE_MATRIX_HPP__
