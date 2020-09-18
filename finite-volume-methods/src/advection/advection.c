/**********************************************************************
 * Nom ............ : advection.c
 * Role ........... : Resolution de l'equation d'advection par des
 *                    methodes aux differences finies.
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "scan_real.h"
#include "methodes.h"
#include "io.h"

#define FILENAME "lax-friedrichs.dat"   /* Nom du fichier ou on
                                         ecrit la solution */


int main(void){

	/* --------------- DECLARATION DES VARIABLES --------------- */
	
	int j = 0;                /* compteur pour l'espace */
	int n = 0;                /* compteur pour le temps */

	int nx = 0;               /* nombre de points en abscisses
                                     pour un t donne */
	int nt = 0;               /* nombre de points en temps
                                     pour un x donne */

	real h = 0.;              /* Pas pour l'abscisses */
	real k = 0.;              /* Pas pour le temps */
	real a = 0.;              /* Vitesse */
	real alpha = 0.;          /* alpha = k*a/h */

	real u_left = 0.;         /* Condition au bord gauche */
	real u_right = 0.;        /* Condition au bord droit */
	real x_left = 0.;         /* Borne inferieure en x */
	real x_right = 0.;        /* Borne superieure en x */
	real t_f = 0.;            /* Temps final */

	real* unew = NULL;        /* Vecteur contenant U^n+1 */
	real* uold = NULL;        /* Vecteur contenant U^n */

	/* ---------------- FIN DE LA DECLARATION ----------------- */


	/* Dialogue avec l'utilisateur */

	printf("Donner la borne inferieure en x: ");
	scan_real(&x_left);
	printf("x_left = %g\n", x_left);

	printf("Donner la borne superieure en x: ");
	scan_real(&x_right);
	printf("x_right = %g\n", x_right);

	printf("Donner la condition au bord gauche u_left: ");
	scan_real(&u_left);
	printf("u_left = %g\n", u_left);

	printf("Donner la condition au bord droit u_right: ");
	scan_real(&u_right);
	printf("u_right = %g\n", u_right);

	printf("Donner le temps final: ");
	scan_real(&t_f);
	printf("t_f = %g\n", t_f);

	printf("Entrer la valeur de la vitesse a: ");
	scan_real(&a);
	printf("a = %g\n", a);

	printf("Entrer la valeur du pas h en espace: ");
	scan_real(&h);
	printf("h = %g\n", h);

	printf("Entrer la valeur du pas k en temps: ");
	scan_real(&k);
	printf("k = %g\n", k);

	/* Calcul du nombre de points en x */
	nx = (int) ((x_right - x_left)/h + 1); /* Ameliorer au cas ou
                                                  t_f est pas int */
	printf("nx = %d\n\n", nx);

	/* Calcul du nombre de points en t */
	nt = (int) (t_f/k + 1);  // Meme remarque qu'au dessus
	printf("nt = %d\n\n", nt);

	/* Calcul et affichage de alpha */
	alpha = k*a/h;
	printf("alpha = %g\n\n", alpha);

	/* Allocation memoire des tableaux et matrices */
	unew = (real*)calloc((size_t)nx, sizeof(real));
	uold = (real*)calloc((size_t)nx, sizeof(real));


	/* ---------- CALCUL + ECRITURE DE LA SOLUTION ------------- */

	remove(FILENAME);           /* Si le fichier existe
                                       deja on le supprime */
	entete(FILENAME, nx, nt); /* Ecriture de l'entete
                                       contenant le nombre
                                       de points en x et
                                       d'instant */

	for(j=0;j<nx;j++){    /* Initialisation du probleme */
		if((x_left + (real)j*h)<=0){
			uold[j] = 1.;
		}
		else if((x_left + (real) j*h)>0){
			uold[j] = 0.;
		}
		else{
			fprintf(stderr, "Erreur !\n");
			exit(EXIT_FAILURE);
		}
	}
	ecrit(FILENAME, nx, uold);

	for(n=0;n<nt-1;n++){
		lax_friedrichs(nx, u_left, u_right, alpha, uold, unew);
		ecrit(FILENAME, nx, unew);
		for(j=0;j<nx;j++){          /* Transfert de unew dans 
                                               uold pour la */
			uold[j] = unew[j];  /* prochaine iteration */
		}
	}

	/* --------------------------------------------------------- */


	/* Desallocation memoire des tableaux et matrices */
	free(unew);
	free(uold);

	exit(EXIT_SUCCESS);
}
