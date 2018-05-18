 #### README pour code en c++ ####


### Utilisation du makefile ###
>> make : lance l'éxécution du main

>> make test: lance toutes les codes de test pour vérifier les fichiers hpp et cpp

>> make clean: nettoie les fichiers de compilation et les executables test

>> make mrproper: nettoie tous les fichiers objets et les executables


### conventions de programmation

Pour plus de simplicité, nous utiliserons les conventions pour les noms de variables, classes, etc...

* les variables simples ont leurs noms commencant en minuscules: (ex) toto, totoNext,...
* les classes ont leurs noms dont la première lettre est une majuscule : (ex) Matrix, Vecteur, 
* les attributs des classes ont leurs noms en minuscules comme les variables mais précédés de _: (ex) _prenom, _livreBibio,
* les noms de namespace sont definis en minuscules



### Programmation ###

Nous allons créer dans un premier temps un mesh -> fonctions de création de mesh 

Nous allons créer des matrices pleines et réutiliser les matrices sparses créées par Victor

Générer une classe fille des matrices pleines = les vecteurs. 

- classe PDE qui va contenir toutes les informations:
	* longueur de coté selon x
	* longueur de coté selon y
	* Temps final T
	* la matrice carré de R^2 de diffusion D 
	* K est la matrice de permeabilité
	* element source f1
	* 
	
	
### Fichiers data sur les matrices ###

** Rmq: Tous les fichiers data sont au format .dat ** 

le fichier est construit comme suit: 
_nbLignes
_nbColonnes
_values enregistrées sur le principe des indices i (ie: par colonne)

### Fichiers data sur les mesh

_xmax
_ymax
_nx
_ny
	
	
