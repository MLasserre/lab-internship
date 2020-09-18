/**********************************************************************
 * Nom ............ : io.c
 * Role ........... : La fonction entete ecrit le nombre de points
 *                    en x et d'instants au
 *                    debut du fichier "filename".
 *                    La fonction ecrit ajoute un tableau de taille
 *                    size_x dans un fichier filename en ecrivant en
 *                    debut et en fin respectivement la condition au
 *                    bord gauche et la condition au bord droit
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "real_precision.h"
#include "io.h"


void entete(const char* filename, const int nb_x, const int nb_t){

	FILE* file = NULL;
	file = fopen(filename, "w");

	fprintf(file, "%d  %d\n", nb_x, nb_t);

	fclose(file);
	return;
}
	
void ecrit(const char* filename, const int size_x, real* mat){
	
	int i = 0;
	FILE* file = NULL;

	file = fopen(filename, "a");

	for(i=0;i<size_x;i++){
		fprintf(file, "%g  ", mat[i]);
	}
	fprintf(file, "\n");

	fclose(file);
}
