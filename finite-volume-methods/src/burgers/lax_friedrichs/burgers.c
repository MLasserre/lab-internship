/**********************************************************************
 * Nom ............ : burgers.c
 * Role ........... : Resolution de l'equation de Burgers avec des
 *                    methodes conservatives.
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "scan_real.h"
#include "methodes.h"
#include "io.h"

#define FILENAME "lax-friedrichs.dat"   /* Nom du fichier ou est
                                           ecrit la solution */


int main(void){

	/* --------------- DECLARATION DES VARIABLES --------------- */
	
	int j = 0;                /* compteur pour l'espace */
	int n = 0;                /* compteur pour le temps */

	int nx = 0;               /* nombre de points en abscisses
                                     pour un t donne */
	int nt = 0;               /* nombre de points en temps
                                     pour un x donne */

	real h = 0.;              /* pas pour l'abscisses */
	real k = 0.;              /* pas pour le temps */
	real alpha = 0.;          /* alpha = k/h */

	real x_left = 0.;         /* Borne inferieure en x */
	real x_right = 0.;        /* Borne superieure en x */
	real t_f = 0.;            /* Temps final */

    real x_dis = 0.;          /* Abscisse de la discontinuite */
    real l_border = 0.;       /* Condition au bord gauche */
    real r_border = 0.;       /* Condition au bord droit */

    real u_left = 0.;         /* Valeur a gauche de la discontinuite */
    real u_right = 0.;        /* Valeur a droite de la discontinuite */

	real* unew = NULL;        /* vecteur contenant U^n+1 */
	real* uold = NULL;        /* vecteur contenant U^n */

	/* ---------------- FIN DE LA DECLARATION ----------------- */


	/* Dialogue avec l'utilisateur */

	printf("Donner la borne inferieure en x: ");
	scan_real(&x_left);
	printf("x_left = %g\n", x_left);

	printf("Donner la borne superieure en x: ");
	scan_real(&x_right);
	printf("x_right = %g\n\n", x_right);

	printf("Donner le temps final: ");
	scan_real(&t_f);
	printf("t_f = %g\n", t_f);

	printf("Entrer la valeur du pas h en espace: ");
	scan_real(&h);
	printf("h = %g\n", h);

	printf("Entrer la valeur du pas k en temps: ");
	scan_real(&k);
	printf("k = %g\n", k);

	printf("Entrer la condition au bord gauche l_border: ");
	scan_real(&l_border);
	printf("l_border = %g\n", l_border);

	printf("Entrer la condition au bord droit r_border: ");
	scan_real(&r_border);
	printf("r_border = %g\n", r_border);
    
	printf("Entrer la valeur a gauche de la discontinuite u_left: ");
	scan_real(&u_left);
	printf("u_left = %g\n", u_left);

	printf("Entrer la valeur a droite de la discontinuite u_right: ");
	scan_real(&u_right);
	printf("u_right = %g\n", u_right);

	printf("Entrer l'abscisse de la discontinuite: ");
	scan_real(&x_dis);
	printf("x_dis = %g\n", x_dis);

	/* Calcul du nombre de points en x */
	nx = (int) (ceil((x_right - x_left)/h + 1)); 
	printf("nx = %d\n\n", nx);

	/* Calcul du nombre de points en t */
	nt = (int) (ceil(t_f/k + 1));
	printf("nt = %d\n\n", nt);

	/* Calcul et affichage de alpha */
	alpha = k/h;
	printf("alpha = %g\n\n", alpha);

	/* Allocation memoire des vecteurs */
	unew = (real*)calloc((size_t)nx, sizeof(real));
	uold = (real*)calloc((size_t)nx, sizeof(real));


	/* ---------- CALCUL + ECRITURE DE LA SOLUTION ------------- */

	remove(FILENAME);           /* Si le fichier existe
                                   deja on le supprime */

    /* Ecriture de l'entete contenant le nombre de points en x et
       d'instant */
	entete(FILENAME, nx-2, nt, x_left, x_right, t_f, h, k); 

    /* Initialisation du probleme */
    uold[0] = l_border;
	for(j=1;j<nx-1;j++){
		if((x_left + (real)j*h)<=x_dis){
			uold[j] = u_left;
		}
		else if((x_left + (real) j*h)>x_dis){
			uold[j] = u_right;
		}
		else{
			fprintf(stderr, "Erreur !\n");
			exit(EXIT_FAILURE);
		}
	}
    uold[nx-1] = r_border;

	ecrit(FILENAME, nx, uold); /* Ecriture de la solution */
                               /* au temps initial */
	for(n=1;n<nt;n++){
		lax_friedrichs(nx, alpha, l_border, r_border, uold, unew);
		ecrit(FILENAME, nx, unew);
		for(j=0;j<nx;j++){      /* Transfert de unew dans */
			uold[j] = unew[j];  /* uold pour la prochaine iteration */
		}
	}

	/* --------------------------------------------------------- */


	/* Desallocation memoire des tableaux */
	free(unew);
	free(uold);

	exit(EXIT_SUCCESS);
}
