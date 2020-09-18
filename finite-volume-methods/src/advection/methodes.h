/**********************************************************************
 * Nom ............ : methodes.h
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#ifndef METHODES
#define METHODES

void backward_euler(const int n, real u_gauche, real u_droite,
                    real alpha, real* anc, real* nouv);
void one_sided_left(const int n, real u_gauche, real u_droite,
                    real alpha, real* anc, real* nouv);
void one_sided_right(const int n, real u_gauche, real u_droite,
                     real alpha, real* anc, real* nouv);
void lax_friedrichs(const int n, real u_gauche, real u_droite,
                    real alpha, real* anc, real* nouv);
void lax_wendroff(const int n, real u_gauche, real u_droite,
                  real alpha, real* anc, real* nouv);
void beam_warming(const int n, real u_gauche_1, real u_gauche_2,
                  real alpha, real* anc, real* nouv);

#endif /* METHODES */
