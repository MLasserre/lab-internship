/**********************************************************************
 * Nom ............ : methodes.h
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#ifndef METHODES
#define METHODES

void backward_euler(const int n, real u_gauche, real u_droite, real M,
                    real k, real h, real* vect_r, real* anc, real* nouv);

#endif /* METHODES */
