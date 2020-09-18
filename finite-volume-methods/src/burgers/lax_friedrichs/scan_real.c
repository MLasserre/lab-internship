/**********************************************************************
 * Nom ............ : scan_real.c
 * Role ........... : Fonction de lecture pour le type real definit
 *                    dans mncs_precision
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "real_precision.h"
#include "scan_real.h"

void scan_real(real* in){

	if(sizeof(real)==sizeof(float))
		scanf("%f", (float*)in);
	else if(sizeof(real)==sizeof(double))
		scanf("%lf", (double*)in);
	else
		fprintf(stderr, "Erreur sur la taille de real\n");
}
