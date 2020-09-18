/**********************************************************************
 * Nom ............ : methodes.c
 * Role ........... : Implementation de la methode de Glimm.
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "methodes.h"


static float compteur = 0;

void glimm(const int N, real l_border, real r_border, real dt,
           real dx, real* anc, real* nouv){

	int j = 0;         /* Compteur pour l'espace */

	real x = 0.;       /* Variable aleatoire a valeurs dans */
	real theta = 0.;   /* Variable aleatoire a valeurs dans [0,1] */

	real s = 0.;       /* Vitesse du choc */

	/* Tirage d'un nombre entre -0.5 et 0.5 */
	theta = (real)(rand()) / (real)(RAND_MAX) - (real)0.5; 

	/* Passage de t a t+1/2 */
	printf("Premiere Phase\n");
	for(j=0;j<N-1;j++){
		x = theta * dx;
		if(anc[j]>anc[j+1]){         /* Cas d'un choc */
			s = (anc[j] + anc[j+1]) / 2;
			printf("sigma[%d] = %g\n", j, s*dt/2);
			if(x<=(s*dt/2)){
				nouv[j] = anc[j];
			}
			else{
				nouv[j] = anc[j+1];
			}
		}
		else{   /* Cas d'une rarefaction */
			printf("Rarefaction !\n");
			if(x<anc[j]*dt/2){
				nouv[j] = anc[j];
			}
			else if(x>anc[j+1]*dt/2){
				nouv[j] = anc[j+1];
			}
			else{
				nouv[j] = 2*x/dt;
			}
		}
	}

	for(j=0;j<N;j++){
		anc[j] = nouv[j];
	}

    /* Tirage d'un autre nombre entre -0.5 et 0.5 */
	theta = (real)(rand())/(real)(RAND_MAX) - (real) 0.5; 


	/* Passage de t a t+1/2 */
	printf("Deuxieme Phase\n");
	nouv[0] = l_border;
	for(j=0;j<N-2;j++){
		x = theta * dx;
		if(anc[j]>anc[j+1]){         /* Cas d'un choc */
			s = (anc[j] + anc[j+1]) / 2;
			printf("sigma[%d] = %g\n", j, s*dt/2);
			if(x<=(s*dt/2)){
				nouv[j+1] = anc[j];
			}
			else{
				nouv[j+1] = anc[j+1];
			}
		}
		else{   /* Cas d'une rarefaction */
			if(x<anc[j]*dt/2){
				nouv[j+1] = anc[j];
			}
			else if(x>anc[j+1]*dt/2){
				nouv[j+1] = anc[j+1];
			}
			else{
				nouv[j+1] = 2*x/dt;
			}
		}
	}
	nouv[N-1] = r_border;
	
	for(j=0;j<N;j++){
		printf("nouv[%d] = %g\n", j, nouv[j]);
	}
	
	compteur++;
}
