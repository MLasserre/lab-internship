/**********************************************************************
 * Nom ............ : burgers.c
 * Role ........... : Resolution de l'equation de Burgers relativiste
 *                    avec une methode de Glimm generalisee.
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <time.h>

#include "real_precision.h"
#include "functions.h"
#include "scan_real.h"
#include "methodes.h"
#include "io.h"

/* Nom du fichier ou on ecrit la solution */
#define FILENAME "glimm.dat"

#define DEFVAR
#include "parametres.h"

int main(void){

	EPSILON = 1.;

    /* Initialisation de la "graine" pour le tirage de la
       variable aleatoire*/
    srand((unsigned int)time(NULL));

	/* --------------- DECLARATION DES VARIABLES --------------- */

	int j = 0;                /* compteur pour l'espace */
	int n = 0;                /* compteur pour le temps */

	int nr = 0;               /* nombre de points en abscisses
                                 pour un t donne */
	int nt = 0;               /* nombre de points en temps
                                 pour un x donne */

	real h = 0.;              /* pas pour l'abscisses */
	real k = 0.;              /* pas pour le temps */

	real r_left = 0.;         /* Borne inferieure en x */
	real r_right = 0.;        /* Borne superieure en x */
	real t_0 = 0.;            /* Temps initial */
	real t_f = 0.;            /* Temps final */

	real r_dis = 0.;          /* Abscisse de la discontinuite */
	real l_border = 0.;       /* Condition au bord gauche */
	real r_border = 0.;       /* Condition au bord droit */

	real v_l = 0.;        /* Valeur a gauche de la discontinuite */
	real v_r = 0.;        /* Valeur a droite de la discontinuite */

	real* vnew = NULL;        /* vecteur contenant V^n+1 */
	real* vold = NULL;        /* vecteur contenant V^n */

	real* vect_r = NULL; /* Vecteur contenant l'intervalle discretise */
	real* vect_R = NULL;      /* Vecteur contenant le rayon pour lequel 
                                 est calculee la solution. Il ne 
                                 correspond plus avec la grille. */

	/* ---------------- FIN DE LA DECLARATION ----------------- */


	/* Dialogue avec l'utilisateur */

	printf("Donner la borne inferieure en r: ");
	scan_real(&r_left);
	printf("r_left = %g\n", r_left);

	printf("Donner la borne superieure en r: ");
	scan_real(&r_right);
	printf("r_right = %g\n\n", r_right);

	printf("Donner le temps initial: ");
	scan_real(&t_0);
	printf("t_0 = %g\n", t_0);

	printf("Donner le temps final: ");
	scan_real(&t_f);
	printf("t_f = %g\n", t_f);

	printf("Donner la valeure de la masse M: ");
	scan_real(&M);
	printf("M = %g\n", M);

	printf("Entrer la valeur du pas h en espace: ");
	scan_real(&h);
	printf("h = %g\n", h);

	printf("Entrer la valeur du pas k en temps: ");
	scan_real(&k);
	printf("k = %g\n", k);

	printf("Entrer la valeur a gauche de la discontinuite v_left: ");
	scan_real(&v_l);
	printf("v_l = %g\n", v_l);

	printf("Entrer la valeur a droit de la discontinuite v_right: ");
	scan_real(&v_r);
	printf("v_right = %g\n", v_r);

	printf("Entrer l'abscisse de la discontinuite: ");
	scan_real(&r_dis);
	printf("r_dis = %g\n", r_dis);

	/* Calcul du nombre de points en x */
	nr = (int) (ceil((r_right - r_left)/h + 1));
	printf("nr = %d\n\n", nr);

	/* Calcul du nombre de points en t */
	nt = (int) (ceil((t_f - t_0)/k + 1));
	printf("nt = %d\n\n", nt);


	/* Allocation memoire des tableaux et matrices */
	vnew = (real*)calloc((size_t)nr, sizeof(real));
	vold = (real*)calloc((size_t)nr, sizeof(real));
	vect_r = (real*)calloc((size_t)nr, sizeof(real));
	vect_R = (real*)calloc((size_t)nr, sizeof(real));


	/* ---------- CALCUL + ECRITURE DE LA SOLUTION ------------- */

	remove(FILENAME);  /* Si le fichier existe deja on le supprime */

	/* Ecriture de l'entete contenant le
       nombre de points en x et d'instant */
	entete(FILENAME, nr-2, nt, r_left, r_right, t_0, t_f, h, k);

	/* Remplissage de vect_r */
	vect_r[j] = r_left;
	for(j=1;j<nr-1;j++){
		vect_r[j] = r_left + (real)j*h;
	}
	vect_r[nr-1] = r_right;

	/* Calcul de la valeur aux bornes */
	l_border = V(v_l, r_dis, r_left);
	r_border = V(v_r, r_dis, r_right);

	/* Initialisation du probleme */
	vold[0] = l_border;
	for(j=1;j<nr-1;j++){
		if(vect_r[j]<=r_dis){
			vold[j] = V(v_l, r_dis, vect_r[j]);
		}
		else if(vect_r[j]>r_dis){
			vold[j] = V(v_r, r_dis, vect_r[j]);
		}
		else{
			fprintf(stderr, "Erreur !\n");
			exit(EXIT_FAILURE);
		}
	}
	vold[nr-1] = r_border;

	for(j=0;j<nr;j++){
		vect_R[j] = vect_r[j];
	}

	printf("Temps n = 0\n");
	for(j=0;j<nr;j++){
		printf("anc[%d] = %g\n", j, vold[j]);
	}
	printf("\n");

	ecrit(FILENAME, nr, vold); /* Ecriture de la solution */
                               /* au temps initial */

	for(n=1;n<nt;n++){
		printf("Temps n = %d\n", n);
		glimm(nr, k, h, r_dis, n*k, vect_R, vect_r, vold, vnew);
		ecrit(FILENAME, nr, vnew);
		for(j=0;j<nr;j++){          /* Transfert de unew dans 
                                       uold pour la */
			vold[j] = vnew[j];      /* prochaine iteration */
		}
	}

	/* --------------------------------------------------------- */

	/* Desallocation memoire des tableaux et matrices */
	free(vect_r);
	free(vect_R);
	free(vnew);
	free(vold);

	exit(EXIT_SUCCESS);
}
