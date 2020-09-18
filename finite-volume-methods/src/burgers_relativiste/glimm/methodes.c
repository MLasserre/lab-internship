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
#include "functions.h"
#include "parametres.h"

/* Methode de dichotomie pour resoudre l'equation pour K dans le cas
   de la rarefaction */
real dichotomie(real epsilon, real inf, real sup, real r,
                real t, real r_dis){
	real mid = 0.;
	real sgn = (r-r_dis >= 0) ? 1. : -1.;

	if(sgn*(R(r,inf) - R(r_dis,inf))/t -1 == 0){
		return inf;
	}
	else if(sgn*(R(r,sup) - R(r_dis,sup))/t -1 == 0){
		return sup;
	}
	else{
		while(fabs(sup-inf) > epsilon){
			mid = (inf + sup) / 2;
			if(sgn*(R(r,mid) - R(r_dis,mid))/t -1 > 0){
				sup = mid;
			}
			else if(sgn*(R(r,mid) - R(r_dis,mid))/t -1 < 0){
				inf = mid;
			}
			else{
				return mid;
			}
		}
	}
	return mid;
}

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

//static int stop = 0;

void glimm(const int nr, real dt, real dr, real r_dis, real t,
           real* vect_R, real* vect_r, real* anc, real* nouv){


	int i = 0;         /* Compteur pour l'espace */
	//int condition = 2;

	real K = 0;
	real l_border = anc[0];
	real r_border = anc[nr-1];

	//printf("l_border = %g\n", l_border);
	//printf("r_border = %g\n", r_border);

	real left = 0;
	real right = 0;

	real theta = 0.;   /* Variable aleatoire a valeurs
                          dans [-0.5,0.5] */
	real sigma_l = 0;
	real sigma_r = 0;
	real sigma = 0;   /* solution des trajectoires pour chaque choc
                             entre deux cellules */

	real* Rnew = NULL;
	real* R = NULL;

	//FILE* file = fopen("sigma.dat", "w");


	/* Allocation memoire des tableaux */
	R = (real*)calloc((size_t)nr-1, sizeof(real));
	Rnew = (real*)calloc((size_t)nr-2, sizeof(real));

	/* Tirage d'un nombre entre -0.5 et 0.5 */
	theta = (real)(rand()) / (real)(RAND_MAX) - (real)0.5; 
	//printf("theta = %g\n", theta);

	/* Remplissage du vecteur R */
	for(i=0;i<nr-1;i++){
		R[i] = vect_r[i] + dr/2 + theta*dr;
	}
	/*
	for(i=0;i<nr;i++){
		printf("vect_R[%d] = %g\n", i, vect_R[i]);
	}
	printf("\n");
	*/


    /* ------------------- METHODE DE GLIMM ----------------------- */

	/* Passage de t a t+1/2 */
	//printf("Premiere phase\n");
	for(i=0;i<nr-1;i++){
		left = V(anc[i], vect_R[i], vect_r[i] + dr/2);
		right = V(anc[i+1], vect_R[i+1], vect_r[i+1] - dr/2);

		//printf("left = %g\n", left);
		//printf("right = %g\n", right);

		/*
		condition = left >= right;
		printf("condition = %d\n", condition);
		if(!condition){printf("left - right = %g\n", left - right);}
		*/

		/* *************** CAS D'UN CHOC **************** */

		//if(left >= right){
		if(left >= right || fabs(left-right) < 0.0000000001){
			sigma = runge_kutta4(vect_r[i] + dr/2, 0, dt/2,
                                    left, right,
                                    vect_r[i] + dr/2);
			//fprintf(file, "%g %g\n", vect_r[i] + dr/2, sigma);
			//printf("R[%d] = %g\n", i, R[i]);
			//printf("sigma[%d] = %g\n", i, sigma);
			
			//printf("anc[%d] = %g\n", i, anc[i]);
			if(R[i]<=sigma){
				nouv[i] = V(anc[i], vect_R[i], R[i]);
				//printf("nouv[%d] = %g\n\n", i, nouv[i]);
			}
			else{
				nouv[i] = V(anc[i+1],vect_R[i+1], R[i]);
				//printf("nouv[%d] = %g\n\n", i, nouv[i]);
			}
		}


		/* ************** CAS D'UNE RAREFACTION *********** */

		else{
			//printf("Rarefaction !\n");
			//printf("Pour i = %d et premiere phase\n\n", i);
			sigma_l = runge_kutta4(vect_r[i] + dr/2, 0, dt/2,
                                    left, left,
                                    vect_r[i] + dr/2);
			//printf("sigma_l[%d] = %g\n", i, sigma_l);
			sigma_r = runge_kutta4(vect_r[i] + dr/2, 0, dt/2,
                                    right, right,
                                    vect_r[i] + dr/2);
			//printf("sigma_r[%d] = %g\n", i, sigma_r);
			if(R[i] < sigma_l){
				nouv[i] = V(anc[i], vect_R[i], R[i]);
				//printf("nouv[%d] = %g\n\n", i, nouv[i]);
			}
			else if(R[i] > sigma_r){
				nouv[i] = V(anc[i+1],vect_R[i+1], R[i]);
				//printf("nouv[%d] = %g\n\n", i, nouv[i]);
			}
			else{
				real sgn = (R[i] >= r_dis) ? 1. : -1;
				K = dichotomie(0.0000001, 0., 1., R[i], t, r_dis);
				nouv[i] = sgn*sqrt(1 - K*K*(1-2*M/R[i]));
				//nouv[i] = 2*R[i]/dt;
			}
		}
	}

	for(i=0;i<nr-1;i++){
	//for(i=0;i<nr-1;i++){
		anc[i] = nouv[i];
	}

    /* Tirage d'un autre nombre entre -0.5 et 0.5 */
	theta = (real)(rand())/(real)(RAND_MAX) - (real) 0.5; 
	//printf("theta = %g\n", theta);
	for(i=0;i<nr-2;i++){
		Rnew[i] = vect_r[i+1] + theta*dr;
		//printf("Rnew[%d] = %g\n", i, Rnew[i]);
	}

	/* Passage de t+1/2 a t+1 */
	//printf("Deuxieme phase\n");
	nouv[0] = l_border;
	//for(i=1;i<nr-1;i++){
	for(i=1;i<nr-1;i++){
		left = V(anc[i-1], R[i-1], vect_r[i]);
		right = V(anc[i], R[i], vect_r[i]);

		//printf("left = %g\n", left);
		//printf("right = %g\n", right);

		/*
		condition = left >= right;
		printf("condition = %d\n", condition);
		if(!condition){printf("left - right = %g\n", left - right);}
		*/

		/* *************** CAS D'UN CHOC **************** */

		//if(left>=right){
		if(left >= right || fabs(left-right) < 0.000000000001){
			sigma = runge_kutta4(vect_r[i], dt/2, dt/2,
                                    left, right, vect_r[i]);
			//fprintf(file, "%g %g\n", vect_r[i] + dr/2, sigma);
			//printf("Rnew[%d] = %g\n", i-1, Rnew[i-1]);
			//printf("sigma[%d] = %g\n", i-1, sigma);
			//printf("anc[%d] = %g\n", i, anc[i]);
			if(Rnew[i-1]<=sigma){
				nouv[i] = V(anc[i-1], R[i-1], Rnew[i-1]);
				//nouv[i+1] = V(anc[i], R[i], Rnew[i]);
				//printf("nouv[%d] = %g\n\n", i, nouv[i]);
			}
			else{
				nouv[i] = V(anc[i], R[i], Rnew[i-1]);
				//nouv[i+1] = V(anc[i+1], R[i+1], Rnew[i]);
				//printf("nouv[%d] = %g\n\n", i, nouv[i]);
			}
		}

		/* ************** CAS D'UNE RAREFACTION *********** */

		else{
			//printf("Rarefaction !\n");
			//printf("Pour i = %d et deuxieme phase\n\n", i);
			sigma_l = runge_kutta4(vect_r[i] + dr/2, 0, dt/2,
                                    left, left,
                                    vect_r[i] + dr/2);
			//printf("sigma_l[%d] = %g\n", i, sigma_l);
			sigma_r = runge_kutta4(vect_r[i] + dr/2, 0, dt/2,
                                    right, right,
                                    vect_r[i] + dr/2);
			//printf("sigma_r[%d] = %g\n", i, sigma_r);
			if(R[i] < sigma_l){
				nouv[i] = V(anc[i], vect_R[i], R[i]);
				//printf("nouv[%d] = %g\n\n", i, nouv[i]);
			}
			else if(R[i] > sigma_r){
				nouv[i] = V(anc[i+1],vect_R[i+1], R[i]);
				//printf("nouv[%d] = %g\n\n", i, nouv[i]);
			}
			else{
				//nouv[i] = 2*R[i]/dt;
				real sgn = (R[i] >= r_dis) ? 1. : -1;
				K = dichotomie(0.0000001, 0., 1., R[i], t+dt/2, r_dis);
				nouv[i] = sgn*sqrt(1 - K*K*(1-2*M/R[i]));
			}

		}
	}
	nouv[nr-1] = r_border;

	vect_R[0] = vect_r[0];
	for(i=1;i<nr-1;i++){
		vect_R[i] = Rnew[i-1];
	}
	vect_R[nr-1] = vect_r[nr-1];

	//fclose(file);
	free(R);
	free(Rnew);
}
