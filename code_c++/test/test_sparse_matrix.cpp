
#include <vector>
#include<iostream>

#include "../sparse_matrix.hpp"


int main ()
{
	std::vector<double> initialValues = {1.1, 2.0, 1.0, 4.0};
	std::vector<int> indvalues = {1,3,5,8};
	SparseMatrice m(3,3, initialValues, indvalues);
	SparseMatrice nullmatrice(2,2);
	
	std::cout<< "test sparse matrice initialisee: m=\n" << m << std::endl;
	std::cout << "test sparse matrice nulle\n" << nullmatrice << std::endl;
	
	nullmatrice = m;
	std::cout << "test operator =: nullmatrice = m:\n" << nullmatrice << std::endl;
	
	//mat = m;
	//std::cout << "verification operateur =:\n" << mat << std::endl;
	
	//std::cout << "test operateur somme: on somme 2x m:\n" << m+m << std::endl;
	
	std::cout << "test operateur () - lecture de m(0,1) = 2.0 ?:\n" << m(0,1) << std::endl;
	std::cout << "test operateur () - lecture de m(1) = m(1, 0) = 1.1 ?:\n" << m(1) << std::endl;
	std::cout << "test operateur () - lecture de m(2,1) = 0?:\n" << m(2,1) << std::endl;
	std::cout << "test operateur () - lecture de m(5) = 1.0:\n" << m(5) << std::endl;
	
	return 0;
}
