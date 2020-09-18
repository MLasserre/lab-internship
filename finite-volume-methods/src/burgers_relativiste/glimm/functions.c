/**********************************************************************
 * Nom .......... : functions.c
 * Role ......... : 
 * Auteur ....... : Marvin Lasserre
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "real_precision.h"
#include "functions.h"
#include "parametres.h"

real V(real v_0, real r_0, real r){
	real sgn = 0.;
	if(v_0 >= 0.){
		sgn = 1.;
	}
	else if(v_0 < 0.){
		sgn = -1.;
	}
	else{
		fprintf(stderr, "Erreur dans la determination de sgn !\n");
		exit(EXIT_FAILURE);
	}
	
	return sgn*(real)sqrt( 1./EPSILON/EPSILON
                          -(1/EPSILON/EPSILON-v_0*v_0)
                          /(1-2*M/r_0)*(1-2*M/r));
}
real R(real r, real K){
	real A = pow(1-K*K, 1.5);
	real B = sqrt(1-K*K);
	real C = sqrt(1 - K*K*(1-2*M/r));

	return (2*M*A*log(r - 2*M) - 2*M*A*log(2*r*C + (2*M - r)*K*K)
          +r*B*C + M*(2-3*K*K)*log(r*B*C + (M-r)*K*K + r))/A;
}

