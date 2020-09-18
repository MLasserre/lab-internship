/**********************************************************************
 * Nom ............ : methodes.c
 * Role ........... : Implementation de la methode de Glimm.
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "methodes.h"
#include "functions.h"
#include "parametres.h"

/* Second membre de l'EDO pour trouver la forme du choc */
real scd_mbr(real t, real y, real v_l, real v_r, real r_dis){
	return (1-2*M/y)*(V(v_l, r_dis, y)+V(v_r, r_dis, y))/2;
}

/* Methode de Runge-Kutta d'ordre 4 */
real runge_kutta4(real u_i, real t_i, real dt, real v_l, real v_r,
                  real r_dis){
	real k1 = 0.;
	real k2 = 0.;
	real k3 = 0.;
	real k4 = 0.;

	k1 = scd_mbr(t_i, u_i, v_l, v_r, r_dis);
	k2 = scd_mbr(t_i + dt/2, u_i + k1*dt/2, v_l, v_r, r_dis);
	k3 = scd_mbr(t_i + dt/2, u_i + k2*dt/2, v_l, v_r, r_dis);
	k4 = scd_mbr(t_i + dt, u_i + k3*dt, v_l, v_r, r_dis);

	return u_i + (k1 + 2*(k2 + k3) + k4)*dt/6;
}


void glimm(const int nr, real dt, real dr, real r_dis, real t,
           real* vect_R, real* vect_r, real* anc, real* nouv){


	int i = 0;      /* Compteur pour l'espace */

	real K = 0;
	real l_border = anc[0];
	real r_border = anc[nr-1];

	real left = 0;
	real right = 0;

	real theta = 0.;   /* Variable aleatoire a valeurs
                          dans [-0.5, 0.5] */
	real sigma_l = 0;
	real sigma_r = 0;
	real sigma = 0;   /* solution des trajectoires pour chaque choc
                             entre deux cellules */

	real* Rnew = NULL;
	real* R = NULL;


	/* Allocation memoire des tableaux */
	R = (real*)calloc((size_t)nr-1, sizeof(real));
	Rnew = (real*)calloc((size_t)nr-2, sizeof(real));

	/* Tirage d'un nombre entre -0.5 et 0.5 */
	theta = (real)(rand()) / (real)(RAND_MAX) - (real)0.5; 

	/* Remplissage du vecteur R */
	for(i=0;i<nr-1;i++){
		R[i] = vect_r[i] + dr/2 + theta*dr;
	}

    /* ------------------- METHODE DE GLIMM ----------------------- */

	/* Passage de t a t+1/2 */
	for(i=0;i<nr-1;i++){
		left = V(anc[i], vect_R[i], vect_r[i] + dr/2);
		right = V(anc[i+1], vect_R[i+1], vect_r[i+1] - dr/2);

		/* *************** CAS D'UN CHOC **************** */

		if(left >= right || fabs(left-right) < 0.0000000001){
			sigma = runge_kutta4(vect_r[i] + dr/2, 0, dt/2,
                                    left, right,
                                    vect_r[i] + dr/2);
			if(R[i]<=sigma){
				nouv[i] = V(anc[i], vect_R[i], R[i]);
			}
			else{
				nouv[i] = V(anc[i+1],vect_R[i+1], R[i]);
			}
		}


		/* ************** CAS D'UNE RAREFACTION *********** */

		else{}
	}

	for(i=0;i<nr-1;i++){
		anc[i] = nouv[i];
	}

    /* Tirage d'un autre nombre entre -0.5 et 0.5 */
	theta = (real)(rand())/(real)(RAND_MAX) - (real) 0.5; 
	for(i=0;i<nr-2;i++){
		Rnew[i] = vect_r[i+1] + theta*dr;
	}

	/* Passage de t+1/2 a t+1 */
	nouv[0] = l_border;
	for(i=1;i<nr-1;i++){
		left = V(anc[i-1], R[i-1], vect_r[i]);
		right = V(anc[i], R[i], vect_r[i]);


		/* *************** CAS D'UN CHOC **************** */

		if(left >= right || fabs(left-right) < 0.000000000001){
			sigma = runge_kutta4(vect_r[i], dt/2, dt/2,
                                    left, right, vect_r[i]);
			if(Rnew[i-1]<=sigma){
				nouv[i] = V(anc[i-1], R[i-1], Rnew[i-1]);
			}
			else{
				nouv[i] = V(anc[i], R[i], Rnew[i-1]);
			}
		}

		/* ************** CAS D'UNE RAREFACTION *********** */

		else{}
	}
	nouv[nr-1] = r_border;

	vect_R[0] = vect_r[0];
	for(i=1;i<nr-1;i++){
		vect_R[i] = Rnew[i-1];
	}
	vect_R[nr-1] = vect_r[nr-1];

	free(R);
	free(Rnew);
}
