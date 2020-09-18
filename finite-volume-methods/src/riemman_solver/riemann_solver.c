/**********************************************************************
 * Nom ............ : riemann_solver.c
 * Role ........... : 
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "riemann_solver.h"
#include "io.h"
#include "parametre.h"


real V(real v_0, real r_dis, real r){

	real sgn = 0;
	if(v_0 >= 0)
		sgn = 1.;
	else
		sgn = -1.;

	return sgn * (real)sqrt( 
                       1./EPSILON/EPSILON
                      -(1./EPSILON/EPSILON -v_0*v_0)
                      /(1.-2*M/r_dis)
                      *(1.-2*M/r)
                     );
}

/* Second membre de l'EDO pour trouver la forme du choc */
real scd_mbr(real t, real y, real v_l, real v_r, real r_dis){
	return (1-2*M/y)
          *(V(v_l, r_dis, y)+V(v_r, r_dis, y))
          /2;
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


void riemann_solver(const char* FILENAME, const int nr, const int nt,
                    real v_l, real v_r, real dt, real t_0, real r_left,
                    real r_right, real r_dis, real* r){

	FILE* file = fopen("sigma.dat", "w");
	int i = 0;  /* Compteur pour l'espace */
	int j = 0;  /* Compteur pour le temps */

	real* sigma = NULL;  /* Vecteur contenant la courbe de choc */
	real* sigma_l = NULL;  /* Vecteurs contenant les courbes */
	real* sigma_r = NULL;  /* limites du domaine de rarefaction */
	real* v = NULL;   /* Vecteur contenant la solution pour un temps */

	/* Allocation memoire des tableaux */
	v = (real*)calloc((size_t)nr, sizeof(real));
	sigma = (real*)calloc((size_t)nt, sizeof(real));
	sigma_l = (real*)calloc((size_t)nt, sizeof(real));
	sigma_r = (real*)calloc((size_t)nt, sizeof(real));

	/* Remplissage du vecteur au temps n a partir de v_l et v_r */
	for(i=0;i<nr;i++){
		if(r[i]<=r_dis){
			v[i] = V(v_l, r_left, r[i]);
		}
		else{
			v[i] = V(v_r, r_right, r[i]);
		}
	}
	ecrit(FILENAME,nr, r, v);
	
	if(v_l>v_r){ /* Cas d'un choc */
		sigma[0] = r_dis;
		for(j=1;j<nt;j++){
			sigma[j] = runge_kutta4(sigma[j-1], t_0 + (real)(j-1)*dt,
                                    dt, v_l, v_r, r_dis);
		}
		for(j=0;j<nt;j++){
			fprintf(file, "%g %g\n", t_0 + (real)j*dt, sigma[j]);
		}

		for(j=1;j<nt;j++){
			for(i=0;i<nr;i++){
				if(r[i] <= sigma[j]){
					v[i] = V(v_l, r_left, r[i]);
				}
				else{
					v[i] = V(v_r, r_right, r[i]);
				}
			}
			ecrit(FILENAME, nr, r, v);
		}
	}
	else{ /* Cas d'une rarefaction */
		sigma_l[0] = r_dis;
		sigma_r[0] = r_dis;

		for(j=1;j<nt;j++){
			sigma_l[j] = runge_kutta4(sigma[j-1], t_0 + (real)(j-1)*dt,
                                    dt, v_l, v_l, r_dis);
		}
		for(j=1;j<nt;j++){
			sigma_r[j] = runge_kutta4(sigma[j-1], t_0 + (real)(j-1)*dt,
                                    dt, v_r, v_r, r_dis);
		}

		for(j=1;j<nt;j++){
			for(i=0;i<nr;i++){
				if(r[i] <= sigma_l[j]){
					v[i] = V(v_l, sigma[j], r[i]);
				}
				else if(r[i] > sigma_l[j] && r[i] < sigma_r[j]){
					v[i] = 0;
				}
				else if(r[i] > sigma_r[j]){
					v[i] = V(v_r, sigma[j], r[i]);
				}
				else{
					fprintf(stderr, "Erreur dans la rarefaction !\n");
					exit(EXIT_FAILURE);
				}
			}
			ecrit(FILENAME, nr, r, v);
		}
	}
	free(sigma_l);
	free(sigma_r);
	free(sigma);
	free(v);
}
