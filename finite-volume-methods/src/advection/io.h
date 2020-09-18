/**********************************************************************
 * Nom ............ : io.h
 * Auteur ......... : Marvin Lasserre
 *********************************************************************/

#ifndef IO_H
#define IO_H

void entete(const char* filename, const int nb_x, const int nb_t);
void ecrit(const char* filename, const int size_x, real* mat);

#endif /* IO_H */
