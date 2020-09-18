/**********************************************************************
 * Nom ............ : burgers.c
 * Role ........... : Resolution de l'equation de Burgers relativiste
 *                    par une methode de Lax-Friedrichs.
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "scan_real.h"
#include "methodes.h"
#include "io.h"

/* Nom du fichier ou on ecrit la solution pour w*/
#define FILENAME_W "lax_friedrichs_w.dat"

/* Nom du fichier ou on ecrit la solution pour la vitesse v*/
#define FILENAME_V "lax_friedrichs_v.dat"   

#define EPSILON 1.

int main(void){

	/* --------------- DECLARATION DES VARIABLES --------------- */
	
	int j = 0;                /* compteur pour l'espace */
	int n = 0;                /* compteur pour le temps */

	int nr = 0;               /* nombre de points en abscisses
                                     pour un t donne */
	int nt = 0;               /* nombre de points en temps
                                     pour un x donne */

	real M = 0.;              /* Masse de l'objet */
	real h = 0.;              /* pas pour l'abscisses */
	real k = 0.;              /* pas pour le temps */
	real alpha = 0.;          /* alpha = k/h */

	real r_left = 0.;         /* Borne inferieure en r */
	real r_right = 0.;        /* Borne superieure en r */
	real t_f = 0.;            /* Temps final */

	real r_dis = 0.;          /* Rayon de la discontinuite */
	real l_border = 0.;       /* Condition au bord gauche */
	real r_border = 0.;       /* Condition au bord droit */

	real u_left = 0.;         /* Valeur a gauche */
	real u_right = 0.;        /* Valeur a droite */

	real* vect_r = NULL;      /* Discretisation de l'axe r */
	real* unew = NULL;        /* vecteur contenant U^n+1 */
	real* uold = NULL;        /* vecteur contenant U^n */

	real* vect_v = NULL;      /* vecteur contenant la vitesse */

	/* ---------------- FIN DE LA DECLARATION ----------------- */


	/* Dialogue avec l'utilisateur */

	printf("Donner la borne inferieure en r: ");
	scan_real(&r_left);
	printf("r_left = %g\n", r_left);

	printf("Donner la borne superieure en r: ");
	scan_real(&r_right);
	printf("r_right = %g\n\n", r_right);

	printf("Donner le temps final: ");
	scan_real(&t_f);
	printf("t_f = %g\n", t_f);

	printf("Entrer la valeur de la masse M de l'objet: ");
	scan_real(&M);
	printf("M = %g\n", M);

	printf("Entrer la valeur du pas h en espace: ");
	scan_real(&h);
	printf("h = %g\n", h);

	printf("Entrer la valeur du pas k en temps: ");
	scan_real(&k);
	printf("k = %g\n", k);

	printf("Entrer la valeur de la donnee initiale a gauche u_left: ");
	scan_real(&u_left);
	printf("u_left = %g\n", u_left);

	printf("Entrer la valeur de la donnee initiale a droite u_right: ");
	scan_real(&u_right);
	printf("u_right = %g\n", u_right);

	printf("Entrer le rayon de la discontinuite: ");
	scan_real(&r_dis);
	printf("r_dis = %g\n", r_dis);

	/* Calcul du nombre de points en x */
	nr = (int) (ceil((r_right - r_left)/h + 1));
	printf("nr = %d\n\n", nr);

	/* Calcul du nombre de points en t */
	nt = (int) (ceil(t_f/k + 1));
	printf("nt = %d\n\n", nt);

	/* Calcul et affichage de alpha */
	alpha = k/h;
	printf("alpha = %g\n\n", alpha);

	/* Allocation memoire des tableaux et matrices */
	vect_v = (real*)calloc((size_t)nr, sizeof(real));
	vect_r = (real*)calloc((size_t)nr, sizeof(real));
	unew = (real*)calloc((size_t)nr, sizeof(real));
	uold = (real*)calloc((size_t)nr, sizeof(real));

	/* Remplissage du vecteur r */
	vect_r[0] = r_left;
	for(j=1;j<nr-1;j++){
		vect_r[j] = r_left + (real)j * h;
	}
	vect_r[nr-1] = r_right;

	/* Evaluation en r_left de la branche gauche */
	if(u_left<=0){
		l_border = -(real)sqrt(1./EPSILON/EPSILON 
                   - (1./EPSILON/EPSILON-u_left*u_left
				   * (1.-2*M/r_dis)*(1.-2*M/r_dis)
                   * (1.-2*M/r_dis)*(1.-2*M/r_dis))
                   / (1.-2*M/r_dis)
                   * (1.-2*M/r_left))
                   / (real)(1.-2*M/r_left)/(real)(1.-2*M/r_left);
	}
	else if(u_left>0){
		l_border = (real)sqrt(1./EPSILON/EPSILON 
                   - (1./EPSILON/EPSILON-u_left*u_left
				   * (1.-2*M/r_dis)*(1.-2*M/r_dis)
                   * (1.-2*M/r_dis)*(1.-2*M/r_dis))
                   / (1.-2*M/r_dis)
                   * (1.-2*M/r_left))
                   / (real)(1.-2*M/r_left)/(real)(1.-2*M/r_left);
	}
	else{
		fprintf(stderr, "Erreur !\n");
		exit(EXIT_FAILURE);
	}

	/* Evaluation en r_right de la branche droite */
	if(u_right<=0){
		r_border = -(real)sqrt(1./EPSILON/EPSILON 
                   - (1./EPSILON/EPSILON-u_right*u_right
				   * (1.-2*M/r_dis)*(1.-2*M/r_dis)
                   * (1.-2*M/r_dis)*(1.-2*M/r_dis))
                   / (1.-2*M/r_dis)
                   * (1.-2*M/r_right))
                   / (real)(1.-2*M/r_right)/(real)(1.-2*M/r_right);
	}
	else if(u_right>0){
		r_border = (real)sqrt(1./EPSILON/EPSILON 
                   - (1./EPSILON/EPSILON-u_right*u_right
				   * (1.-2*M/r_dis)*(1.-2*M/r_dis)
                   * (1.-2*M/r_dis)*(1.-2*M/r_dis))
                   / (1.-2*M/r_dis)
                   * (1.-2*M/r_right))
                   / (real)(1.-2*M/r_right)/(real)(1.-2*M/r_right);
	}
	else{
		fprintf(stderr, "Erreur !\n");
		exit(EXIT_FAILURE);
	}


	/* ---------- CALCUL + ECRITURE DE LA SOLUTION ------------- */

	remove(FILENAME_W); /* Si le fichier existe deja on le supprime */
	remove(FILENAME_V);

	/* Ecriture de l'entete contenant le nombre
	   de points en x et d'instant */
	entete(FILENAME_W, nr-2, nt, r_left, r_right, t_f, h, k);
	entete(FILENAME_V, nr-2, nt, r_left, r_right, t_f, h, k);

	/* Initialisation du probleme */
	uold[0] = l_border;
	for(j=1;j<nr;j++){
		if(vect_r[j] <= r_dis && u_left < 0){
			uold[j] = -(real)sqrt(1./EPSILON/EPSILON 
                           - (1./EPSILON/EPSILON-u_left*u_left
                           *(1-2*M/r_dis)*(1-2*M/r_dis)
                           *(1-2*M/r_dis)*(1-2*M/r_dis))
                           / (1.-2*M/r_dis)
                           * (1.-2*M/vect_r[j]))
                           / (1-2*M/vect_r[j])/(1-2*M/vect_r[j]);
		}
		else if(vect_r[j] <= r_dis && u_left >= 0){
			uold[j] = (real)sqrt(1./EPSILON/EPSILON 
                           - (1./EPSILON/EPSILON-u_left*u_left
                           *(1-2*M/r_dis)*(1-2*M/r_dis)
                           *(1-2*M/r_dis)*(1-2*M/r_dis))
                           / (1.-2*M/r_dis)
                           * (1.-2*M/vect_r[j]))
                           / (1-2*M/vect_r[j])/(1-2*M/vect_r[j]);
		}
		else if(vect_r[j] > r_dis && u_right < 0){
			uold[j] = - (real)sqrt(1./EPSILON/EPSILON 
                      - (1./EPSILON/EPSILON-u_right*u_right
                      * (1-2*M/r_dis)*(1-2*M/r_dis)
                      * (1-2*M/r_dis)*(1-2*M/r_dis))
                      / (1.-2*M/r_dis)
                      * (1.-2*M/vect_r[j]))
	              / (1-2*M/vect_r[j])/(1-2*M/vect_r[j]);
		}
		else if(vect_r[j] > r_dis && u_right >= 0){
			uold[j] = (real)sqrt(1./EPSILON/EPSILON 
                    - (1./EPSILON/EPSILON-u_right*u_right
                    * (1-2*M/r_dis)*(1-2*M/r_dis)
                    * (1-2*M/r_dis)*(1-2*M/r_dis))
                    / (1.-2*M/r_dis)
                    * (1.-2*M/vect_r[j]))
                    / (1-2*M/vect_r[j])/(1-2*M/vect_r[j]);
		}
		else{
			fprintf(stderr, "Erreur !\n");
			exit(EXIT_FAILURE);
		}
	}
	uold[nr-1] = r_border;

	for(j=0;j<nr;j++){
		vect_v[j] = (1-2*M/vect_r[j])*(1-2*M/vect_r[j])*uold[j];
	}


	ecrit(FILENAME_W, nr, uold);  /* Ecriture de la donnee initiale */
	ecrit(FILENAME_V, nr, vect_v);

	for(n=1;n<nt;n++){
		lax_friedrichs(nr, M, h, k, l_border, r_border,
                       vect_r, uold, unew);
		for(j=0;j<nr;j++){
			vect_v[j] = (1-2*M/vect_r[j])*(1-2*M/vect_r[j])*unew[j];
		}
		ecrit(FILENAME_W, nr, unew);
		ecrit(FILENAME_V, nr, vect_v);
		for(j=0;j<nr;j++){          /* Transfert de unew dans 
                                       uold pour la */
			uold[j] = unew[j];      /* prochaine iteration */
		}
	}

	/* --------------------------------------------------------- */


	/* Desallocation memoire des tableaux et matrices */
	free(unew);
	free(uold);
	free(vect_v);
	free(vect_r);

	exit(EXIT_SUCCESS);
}
