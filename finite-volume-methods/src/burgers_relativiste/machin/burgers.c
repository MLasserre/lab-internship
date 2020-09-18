/**********************************************************************
 * Nom ............ : burgers.c
 * Role ........... : Resolution de l'equation de Burgers relativiste.
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "functions.h"
#include "scan_real.h"
#include "methodes.h"
#include "io.h"

#define FILENAME "euler.dat"   /* Nom du fichier ou on
                                           ecrit la solution */
#define DEFVAR
#include "parametres.h"

int main(void){

	EPSILON = 1.;

	/* --------------- DECLARATION DES VARIABLES --------------- */
	
	int j = 0;                /* compteur pour l'espace */
	int n = 0;                /* compteur pour le temps */

	int nr = 0;               /* nombre de points en abscisses
                                     pour un t donne */
	int nt = 0;               /* nombre de points en temps
                                     pour un x donne */

	real h = 0.;              /* pas pour l'abscisses */
	real k = 0.;              /* pas pour le temps */

	real r_left = 0.;         /* Borne inferieure en r */
	real r_right = 0.;        /* Borne superieure en r */
	real t_0 = 0.;            /* Temps initial */
	real t_f = 0.;            /* Temps final */

	real r_dis = 0.;          /* Rayon de la discontinuite */
	real l_border = 0.;       /* Condition au bord gauche */
	real r_border = 0.;       /* Condition au bord droit */

	real u_left = 0.;         /* Valeur a gauche */
	real u_right = 0.;        /* Valeur a droite */

	real* vect_r = NULL;      /* Discretisation de l'axe r */

	real* unew = NULL;        /* vecteur contenant U^n+1 */
	real* uold = NULL;        /* vecteur contenant U^n */

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

	printf("Entrer la valeur de la masse M de l'objet: ");
	scan_real(&M);
	printf("M = %g\n", M);

	printf("Entrer la valeur du pas h en espace: ");
	scan_real(&h);
	printf("h = %g\n", h);

	printf("Entrer la valeur du pas k en temps: ");
	scan_real(&k);
	printf("k = %g\n", k);

	/*
	printf("Entrer la condition au bord gauche l_border: ");
	scan_real(&l_border);
	printf("l_border = %g\n", l_border);

	printf("Entrer la condition au bord droit r_border: ");
	scan_real(&r_border);
	printf("r_border = %g\n", r_border);
	*/

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
	nt = (int) (ceil((t_f - t_0)/h + 1));
	printf("nt = %d\n\n", nt);

	/* Allocation memoire des tableaux et matrices */
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
	l_border = V(u_left, r_dis, r_left);
	r_border = V(u_right, r_dis, r_right);


	/* ---------- CALCUL + ECRITURE DE LA SOLUTION ------------- */

	remove(FILENAME); /* Si le fichier existe deja on le supprime */

    /* Ecriture de l'entete contenant le nombre
       de points en x et d'instant */
	entete(FILENAME, nr-2, nt, r_left, r_right, 0, t_f, h, k);

	/* Initialisation du probleme */
	uold[0] = l_border;
	for(j=1;j<nr-1;j++){
		if(vect_r[j] <= r_dis){
			uold[j] = V(u_left, r_dis, vect_r[j]);
		}
		else if(vect_r[j] > r_dis){
			uold[j] = V(u_right, r_dis, vect_r[j]);
		}
		else{
			fprintf(stderr, "Erreur !\n");
			exit(EXIT_FAILURE);
		}
	}
	uold[nr-1] = r_border;

	ecrit(FILENAME, nr, uold);

	for(n=1;n<nt;n++){
		backward_euler(nr, l_border, r_border, M,
                       k, h, vect_r, uold, unew);
		ecrit(FILENAME, nr, unew);
		for(j=0;j<nr;j++){          /* Transfert de unew dans 
                                       uold pour la */
			uold[j] = unew[j];      /* prochaine iteration */
		}
	}

	/* --------------------------------------------------------- */


	/* Desallocation memoire des tableaux et matrices */
	free(unew);
	free(uold);
	free(vect_r);

	exit(EXIT_SUCCESS);
}
