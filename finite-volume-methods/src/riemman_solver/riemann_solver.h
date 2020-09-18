/**********************************************************************
 * Nom ............ : methodes.h
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#ifndef METHODES
#define METHODES

real V(real v_0, real r_dis, real r);
real scd_mbr(real t, real y, real v_l, real v_r, real r_dis);
real runge_kutta4(real u_i, real t_i, real dt, real v_l, real v_r,
                  real r_dis);
void riemann_solver(const char* FILENAME, const int nr, const int nt,
                    real v_l, real v_r, real dt, real t_0, real r_left,
                    real r_right, real r_dis, real* r);

#endif /* METHODES */
