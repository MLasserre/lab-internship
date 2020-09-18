/**********************************************************************
 * Nom ............ : methodes.c
 * Role ........... : 
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "methodes.h"

void godunov(const int N, real alpha, real M, real l_border,
             real r_border, real r_0, real* anc, real* nouv){

	int i = 0;   /* Compteur pour le vecteur flux de taille N-1 */
	real s = 0.; /* Vitesse du choc */

	real f[N-1]; /* Vecteur contenant le flux aux interfaces */
	for(i=0;i<N-1;i++){
		f[i] = 0.;
	}
	

	for(i=0;i<N-1;i++){
		if(anc[i]>anc[i+1]){  /* Cas d'un choc */
			s = (1-2*M/r_0)
               *(1-2*M/r_0)
               *(1-2*M/r_0)
               *(anc[i] + anc[i+1]) / 2;
			if(s<0){
				f[i] = anc[i+1]*anc[i+1]/2;
			}
			else{
				f[i] = anc[i]*anc[i]/2;
			}
		}
		else if(anc[i]<=anc[i+1]){ /* Cas d'une rarefaction */
			if((1-2*M/r_0)*(1-2*M/r_0)*(1-2*M/r_0)*anc[i] > 0){
				f[i] = anc[i]*anc[i]/2;
			}
			else if((1-2*M/r_0)*(1-2*M/r_0)*(1-2*M/r_0)*anc[i+1] < 0){
				f[i] = anc[i+1]*anc[i+1]/2;
			}
			else{
				f[i] = 0.;
			}
		}
		else{
			fprintf(stderr, "Erreur dans la methode de Godunov\n");
			exit(EXIT_FAILURE);
		}
	}

	for(i=1;i<N-1;i++){  /* Resolution */
		nouv[i] = anc[i] - alpha*(f[i] - f[i-1]);
	}

	nouv[0] = l_border;
	nouv[N-1] = r_border;
}
