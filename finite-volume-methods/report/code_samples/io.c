/**********************************************************************
 * Nom ............ : io.c
 * Role ........... : La fonction entete ecrit en debut de fichier
 *                    le nombre de points en x et en temps, les
 *                    bornes du domaine en espace, le temps final
 *                    ainsi que la valeur des pas.
 *                    La fonction ecrit ecrit un tableau de taille
 *                    size_x dans un fichier en ecrivant en
 *                    debut et en fin respectivement la condition au
 *                    bord gauche et la condition au bord droit
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "real_precision.h"
#include "io.h"


void entete(const char* filename, const int nb_x, const int nb_t,
            const real x_left, const real x_right, const real t_f,
            const real x_step, const real t_step){

	FILE* file = NULL;
	file = fopen(filename, "w");

	fprintf(file, "%d  %d %g %g %g %g %g\n", nb_x, nb_t, x_left,
            x_right, t_f, x_step, t_step);

	fclose(file);
	return;
}
	
void ecrit(const char* filename, const int size_x, real* mat){
	
	int i = 0;
	FILE* file = NULL;

	file = fopen(filename, "a");

	for(i=1;i<size_x-1;i++){
		fprintf(file, "%g  ", mat[i]);
	}
	fprintf(file, "\n");

	fclose(file);
}
