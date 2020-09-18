/**********************************************************************
 * Nom ............ : methodes.c
 * Role ........... : 
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "methodes.h"

real num_flux_function(real alpha, real a, real b){
	return (a*a + b*b)/4 - (b - a)/(2*alpha);
}

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
/*
	for(j=1;j<n-1;j++){
		nouv[j] = (anc[j+1] + anc[j-1])/2
                  - alpha/2 * (anc[j+1]*anc[j+1]-anc[j-1]*anc[j-1])/2
                  - alpha/(1-2*M/(r_left+((real)j-(real)0.5)*dr))/2
                  + alpha/(1-2*M/(r_left+((real)j+0.5)*dr))/2
                  + alpha*M/(r_left+((real)j+0.5)*dr)
                  * ((anc[j+1]*anc[j+1] + anc[j]*anc[j])/2
                     - 1/alpha*(anc[j+1]-anc[j]))
                  - alpha*M/(r_left+((real)j- (real)0.5)*dr)
                  * ((anc[j]*anc[j] + anc[j-1]*anc[j-1])/2
                     - 1/alpha*(anc[j]-anc[j-1]));

	}
*/
	nouv[0] = l_border; /* A ameliorer au cas ou l'intervalle est 
                        strictement positif ou negatif */
	nouv[n-1] = r_border;
}

