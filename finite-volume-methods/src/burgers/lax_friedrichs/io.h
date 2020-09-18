/**********************************************************************
 * Nom ............ : io.h
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#ifndef IO_H
#define IO_H

void entete(const char* filename, const int nb_x, const int nb_t,
            const real x_left, const real x_right, const real t_f,
            const real x_step, const real t_step);

void ecrit(const char* filename, const int size_x, real* mat);

#endif /* IO_H */
