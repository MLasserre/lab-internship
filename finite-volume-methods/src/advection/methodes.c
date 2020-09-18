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


void backward_euler(const int n, real u_gauche, real u_droite,
                    real alpha, real* anc, real* nouv){
	int j = 0;
	for(j=1;j<n-1;j++){
		nouv[j] = anc[j] - alpha*(anc[j+1] - anc[j-1])/2;
	}
	nouv[0] = u_gauche;
	nouv[n-1] = u_droite;
}

void one_sided_left(const int n, real u_gauche, real u_droite,
                    real alpha, real* anc, real* nouv){
	int j = 0;
	for(j=1;j<n-1;j++){
		nouv[j] = anc[j] - alpha*(anc[j] - anc[j-1]);
	}
	nouv[0] = u_gauche;
	nouv[n-1] = u_droite;
}

void one_sided_right(const int n, real u_gauche, real u_droite,
                     real alpha, real* anc, real* nouv){
	int j = 0;
	for(j=1;j<n-1;j++){
		nouv[j] = anc[j] - alpha*(anc[j+1] - anc[j]);
	}
	nouv[0] = u_gauche;
	nouv[n-1] = u_droite;
}

void lax_friedrichs(const int n, real u_gauche, real u_droite,
                    real alpha, real* anc, real* nouv){
	int j = 0;
	for(j=1;j<n-1;j++){
		nouv[j] = (anc[j-1] + anc[j+1])/2 -
                          alpha*(anc[j+1] - anc[j-1])/2;
	}
	nouv[0] = u_gauche;
	nouv[n-1] = u_droite;
}

void lax_wendroff(const int n, real u_gauche, real u_droite,
                  real alpha, real* anc, real* nouv){
	int j = 0;
	for(j=1;j<n-1;j++){
		nouv[j] = anc[j] - alpha*(anc[j+1] - anc[j-1])/2
                        + alpha*alpha*(anc[j+1] -2*anc[j] + anc[j-1])/2;
	}
	nouv[0] = u_gauche;
	nouv[n-1] = u_droite;
}

void beam_warming(const int n, real u_gauche_1, real u_gauche_2,
                  real alpha, real* anc, real* nouv){
	int j = 0;
	for(j=2;j<n;j++){
		nouv[j] = anc[j]
                        - alpha*(3*anc[j] - 4*anc[j-1] + anc[j-2])/2
                        + alpha*alpha*(anc[j] -2*anc[j-1] + anc[j-2])/2;
	}
	nouv[0] = u_gauche_1;
	nouv[1] = u_gauche_2;
}


