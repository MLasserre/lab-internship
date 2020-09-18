/**********************************************************************
 * Nom ............ : methodes.h
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#ifndef METHODES
#define METHODES

void lax_friedrichs(const int n, real alpha, real l_border,
                    real r_border, real* anc, real* nouv);

#endif /* METHODES */
