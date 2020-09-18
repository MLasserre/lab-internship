/**********************************************************************
 * Nom ............ : methodes.h
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#ifndef METHODES
#define METHODES

void godunov(const int N, real alpha, real M, real l_border,
             real r_border, real r_0, real* anc, real* nouv);
			
#endif /* METHODES */
