/**********************************************************************
 * Nom ............ : methodes.h
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#ifndef METHODES
#define METHODES

real dichotomie(real epsilon, real inf, real sup, real r,
                real t, real r_dis);
real scd_mbr(real t, real y, real v_l, real v_r, real r_dis);
real runge_kutta4(real u_i, real t_i, real dt, real v_l, real v_r,
                  real r_dis);
void glimm(const int nr, real dt, real dx, real r_dis, real t,
           real* vect_R, real* vect_r, real* anc, real* nouv);
			
#endif /* METHODES */
