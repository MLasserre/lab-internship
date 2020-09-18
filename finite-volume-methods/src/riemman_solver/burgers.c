/**********************************************************************
 * Nom ............ : burgers.c
 * Role ........... : Resolution de l'equation de Burgers relativiste.
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "scan_real.h"
#include "riemann_solver.h"
#include "io.h"

#define FILENAME "solution_riemann.dat"   /* Nom du fichier ou on
                                           ecrit la solution */
#define DEFVAR
#include "parametre.h"

int main(void){

	EPSILON = 1.;             /* Definition de EPSILON a 1 */

	/* --------------- DECLARATION DES VARIABLES --------------- */

	int i = 0;                /* compteur pour l'espace */

	int nr = 0;               /* nombre de points en abscisses
                                     pour un t donne */
	int nt = 0;               /* nombre de points en temps
                                     pour un x donne */

	real dr = 0.;              /* pas pour l'abscisses */
	real dt = 0.;              /* pas pour le temps */

	real r_left = 0.;         /* Borne inferieure en r */
	real r_right = 0.;        /* Borne superieure en r */
	real t_0 = 0.;            /* Temps initial */
	real t_f = 0.;            /* Temps final */

	real r_dis = 0.;          /* Rayon de la discontinuite */

	real v_left = 0.;         /* Valeur a gauche */
	real v_right = 0.;        /* Valeur a droite */

	real* vect_r = NULL;      /* Discretisation de l'axe r */

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

	printf("Entrer la valeur du pas dr en espace: ");
	scan_real(&dr);
	printf("dr = %g\n", dr);

	printf("Entrer la valeur du pas dt en temps: ");
	scan_real(&dt);
	printf("dt = %g\n", dt);

	printf("Entrer la valeur de la donnee initiale a gauche v_left: ");
	scan_real(&v_left);
	printf("v_left = %g\n", v_left);

	printf("Entrer la valeur de la donnee initiale a droite v_right: ");
	scan_real(&v_right);
	printf("v_right = %g\n", v_right);

	printf("Entrer le rayon de la discontinuite: ");
	scan_real(&r_dis);
	printf("r_dis = %g\n", r_dis);

	/* Calcul du nombre de points en r */
	nr = (int) (ceil((r_right - r_left)/dr + 1));
	printf("nr = %d\n\n", nr);

	/* Calcul du nombre de points en t */
	nt = (int) (ceil((t_f - t_0)/dt + 1));
	printf("nt = %d\n\n", nt);

	/* Allocation memoire des tableaux et matrices */
	vect_r = (real*)calloc((size_t)nr, sizeof(real));

	/* Remplissage du vecteur r */
	vect_r[0] = r_left;
	for(i=1;i<nr-1;i++){
		vect_r[i] = r_left + (real)i * dr;
	}
	vect_r[nr-1] = r_right;

	v_left = V(v_left, r_dis, vect_r[0]);
	v_right = V(v_right, r_dis, vect_r[nr-1]);

	remove(FILENAME);
	entete(FILENAME, nr-2, nt, r_left, r_right, t_0, t_f, dr, dt);
	riemann_solver(FILENAME, nr, nt, v_left, v_right, dt, t_0,
                   r_left, r_right, r_dis, vect_r);


	/* Desallocation memoire des tableaux et matrices */
	free(vect_r);

	exit(EXIT_SUCCESS);
}
