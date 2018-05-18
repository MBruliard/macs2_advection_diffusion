##
# BRULIARD Margaux/ RIGAL Mathieu
# Projet Numérique 
##


# variables
header = matrix.hpp sparse_matrix.hpp vector.hpp mesh.hpp

# commandes du makefile
default: main.exe


test: test_matrix.exe test_vector.exe test_mesh.exe test_sparse_matrix.exe


clean: 
	rm *.o test_*.exe
	
mrproper:
	rm *.o *.exe
	

# compilation des fichiers .cpp associes aux .hpp
matrix.o: matrix.cpp matrix.hpp
	g++ -Wall -c -std=c++14 matrix.cpp
	
sparse_matrix.o: sparse_matrix.cpp sparse_matrix.hpp
	g++ -Wall -c -std=c++14 sparse_matrix.cpp
	
vector.o: vector.cpp vector.hpp matrix.hpp
	g++ -Wall -c -std=c++14 vector.cpp

mesh.o: mesh.cpp mesh.hpp
	g++ -Wall -c -std=c++14 mesh.cpp



# compilation du fichier principal
main.exe:
	g++ -Wall -o main.exe -std=c++14 main.cpp



# compilation des fichiers test
test_mesh.o: test/test_mesh.cpp mesh.hpp
	g++ -c -std=c++14 test/test_mesh.cpp
	
test_mesh.exe: test_mesh.o mesh.o mesh.hpp
	g++ -Wall -o test_mesh.exe -std=c++14 test_mesh.o mesh.o

test_matrix.o: test/test_matrix.cpp matrix.hpp
	g++ -c -std=c++14 test/test_matrix.cpp

test_matrix.exe: test_matrix.o matrix.o matrix.hpp
	g++ -Wall -o test_matrix.exe -std=c++14 test_matrix.o matrix.o
	
	
test_sparse_matrix.o: test/test_sparse_matrix.cpp sparse_matrix.hpp
	g++ -c -std=c++14 test/test_sparse_matrix.cpp

test_sparse_matrix.exe: test_sparse_matrix.o sparse_matrix.o sparse_matrix.hpp
	g++ -Wall -o test_sparse_matrix.exe -std=c++14 test_sparse_matrix.o sparse_matrix.o


	
test_vector.o: test/test_vector.cpp vector.hpp matrix.hpp
	g++ -c -std=c++14 test/test_vector.cpp

test_vector.exe: test_vector.o vector.o matrix.o
	g++ -Wall -o test_vector.exe -std=c++14 test_vector.o vector.o matrix.o
	

