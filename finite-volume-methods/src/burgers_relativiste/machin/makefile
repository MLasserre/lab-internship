# makefile  pour C99 avec  bibliotheques openblas/lapack
RM = /bin/rm -f
CC = gcc

# utiliser une variable d'environnement pour le repertoire openblas-lapack
LA_SYSTEM = ${REP_LAPACK}
#
# pour mandriva-2010 32 bits (machine # virtuelle)
# LA_SYSTEM = /home/lefrere/M1/mdv2010-32
#
# pour sappli 64 bits
# LA_SYSTEM = /home/lefrere/M1/sappli1-64
#
# partie commune 
#
LA_INCPATH=$(LA_SYSTEM)/include
LA_LIBPATH=$(LA_SYSTEM)/lib
#
#LD = gcc
LD = gfortran
#GFPATH = /usr/lib/gcc/i386-redhat-linux5E/4.1.2
#LDFLAGS = -L$(LA_LIBPATH) -L$(GFPATH) -static
LDFLAGS = -L$(LA_LIBPATH) -static
LIBS= -lopenblas -lpthread -lm
CFLAGS = -std=c99 -Wextra -Wall -Wunused -Wshadow -Wconversion -Wmissing-prototypes \
         -ffloat-store -pedantic  -I$(LA_INCPATH)

# sources-------------Editer les 2 lignes ci-dessous-------
MAIN= burgers.c
SUB= methodes.c scan_real.c io.c functions.c
#----------------------------------------------------------
HEAD = $(SUB:.c=.h)
OBJ_MAIN = $(MAIN:.c=.o)
OBJ_SUB = $(SUB:.c=.o)
EXEC = $(MAIN:.c=) 

$(EXEC): $(OBJ_MAIN) $(OBJ_SUB)
	$(LD) $(LDFLAGS) $(OBJ_MAIN) $(OBJ_SUB) $(LIBS) -o $@

$(OBJ_MAIN): $(MAIN) $(HEAD) real_precision.h
	$(CC) $(CFLAGS) -c  $< -o $@

# ajout des dependances exceptionnelles
burgers.o : parametres.h
functions.o : parametres.h

#matop.h : matrices.h


%.o: %.c %.h real_precision.h
	$(CC) $(CFLAGS) -c  $< -o $@

clean: 
	$(RM) *.o *~

superclean:
	$(RM) *.o *.out *~ $(EXEC)
   
go : $(EXEC)
	./$(EXEC) < entree.in
