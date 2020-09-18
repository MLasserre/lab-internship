/**********************************************************************
 * Nom ............ : methodes.c
 * Role ........... : Definition de differents schemas numeriques aux
 *                    differences finies : Backward-Euler, One-sided,
 *                    Lax-Friedrichs, Lax-Wendroff et Beam-Warming      
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "methodes.h"


void backward_euler(const int n, real u_gauche, real u_droite, real M,
                    real k, real h, real* vect_r, real* anc, real* nouv){
	int j = 0;
	for(j=1;j<n-1;j++){
		nouv[j] = anc[j] - k*anc[j]*(1-2*M/vect_r[j])
                            *(anc[j+1] - anc[j-1])/2/h
                         + M*k*(anc[j]*anc[j] -1)/vect_r[j]/vect_r[j];
	}
	nouv[0] = u_gauche;
	nouv[n-1] = u_droite;
}
