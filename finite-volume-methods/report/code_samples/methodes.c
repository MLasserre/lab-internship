/**********************************************************************
 * Nom ............ : methodes.c
 * Role ........... : Implementation de la methode de Lax-Friedrichs
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "methodes.h"

/* Fonction de flux numerique */
real num_flux_function(real alpha, real a, real b){
	return (a*a + b*b)/4 - (b - a)/(2*alpha);
}

/* Methode de Lax-Friedrichs */
void lax_friedrichs(const int n, real M, real dr, real dt,
                    real l_border, real r_border, real* r, real* anc,
                    real* nouv){

	int j = 0;
	real alpha = dt/dr;

	for(j=1;j<n-1;j++){
		nouv[j] = anc[j]
                 +alpha*(1./(2*(1.-2*M/(r[j]+(real)0.5*dr)))
                        -1./(2*(1.-2*M/(r[j]-(real)0.5*dr))))
		         -alpha*((1.-2*M/(r[j]+(real)0.5*dr))
		                *(1.-2*M/(r[j]+(real)0.5*dr))
		                *(1.-2*M/(r[j]+(real)0.5*dr))
                        *num_flux_function(alpha, anc[j], anc[j+1])
                        -(1.-2*M/(r[j]-(real)0.5*dr))
                        *(1.-2*M/(r[j]-(real)0.5*dr))
		                *(1.-2*M/(r[j]-(real)0.5*dr))
                         *num_flux_function(alpha, anc[j-1], anc[j]));
	}
	nouv[0] = l_border; 
	nouv[n-1] = r_border;
}

