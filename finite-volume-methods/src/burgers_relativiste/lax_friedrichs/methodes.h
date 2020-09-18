/**********************************************************************
 * Nom ............ : methodes.h
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#ifndef METHODES
#define METHODES

real num_flux_function(real alpha, real a, real b);

void lax_friedrichs(const int n, real M, real dr, real dt,
                    real l_border, real r_border, real* r, real* anc,
                    real* nouv);

#endif /* METHODES */
