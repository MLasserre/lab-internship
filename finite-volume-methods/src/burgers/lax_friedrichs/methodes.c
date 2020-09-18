/**********************************************************************
 * Nom ............ : methodes.c
 * Role ........... : Implementation de la methode de Lax-Friedrichs
 *                    pour la resolution de l'equation de Burgers.
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "methodes.h"


void lax_friedrichs(const int n, real alpha, real l_border,
                    real r_border, real* anc, real* nouv){

	int j = 0;
	for(j=1;j<n-1;j++){
		nouv[j] = (anc[j-1] + anc[j+1])/2
                  - alpha/4 * (anc[j+1]*anc[j+1] - anc[j-1]*anc[j-1]);
	}
	nouv[0] = l_border;
	nouv[n-1] = r_border;
}

