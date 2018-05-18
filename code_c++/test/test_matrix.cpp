
#include <vector>
#include<iostream>

#include "../matrix.hpp"


int main ()
{
	Matrice mat(3,3);
	Matrice matLecture("test/essai_lecture_matrice.dat");
	std::vector<double> initialValues = {1.0, 2.0, 3.0, 4.0};
	Matrice m(2,2, initialValues);
	
	std::cout << "test matrice nulle\n"<< mat << std::endl;
	std::cout<< "test matrice initialisee: m=\n" << m << std::endl;
	
	//mat = m;
	//std::cout << "verification operateur =:\n" << mat << std::endl;
	
	//std::cout << "test operateur somme: on somme 2x m:\n" << m+m << std::endl;
	
	std::cout << "test operateur () - lecture de m(0,1):\n" << matLecture(2) << " " << matLecture(0, 1) << std::endl;
	
	std::cout << "test lecture matrice stockee dans un fichier:\n" << matLecture << std::endl;
	
	return 0;
}
