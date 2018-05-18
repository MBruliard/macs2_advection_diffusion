
#include <vector>
#include<iostream>

#include "../vector.hpp"


int main ()
{
	Vecteur vect(3);
	std::vector<double> val = {1.0, 1.1, 2.1};
	Vecteur vect2(3, val);
	std::vector<double> valMat = {2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0., 0.0, 2.0};
	Matrice m(3,3,valMat);

	
	std::cout << "test vect nulle\n"<< vect << std::endl;
	std::cout<< "test vect initialisee: vect2=\n" << vect2 << std::endl;
	
	vect = vect2;
	
	std::cout << "test operateur somme: on somme 2x m:\n" << vect2+vect2 << std::endl;
	
	std::cout << "test operateur () - lecture de vect2(1):\n" << vect2(1) << std::endl;
	
	std::cout << "matrice m:\n" << m << "test produit matrice vecteur:\n" << m*vect2 << std::endl;
	
	
	std::cout << "test produit vecteur*alpha:\n" << 5.2*vect2 << std::endl;
	
	
	
	
	return 0;
}
