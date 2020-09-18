RM = /bin/rm -f
CC = gcc

#LD = gcc
LD = gfortran
LDFLAGS = -static
LIBS= -lpthread -lm
CFLAGS = -std=c99 -Wextra -Wall -Wunused -Wshadow -Wconversion -Wmissing-prototypes \
         -ffloat-store -pedantic

# sources-------------Editer les 2 lignes ci-dessous-------
MAIN= burgers.c
SUB= methodes.c scan_real.c io.c
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

#matop.h : matrices.h


%.o: %.c %.h real_precision.h
	$(CC) $(CFLAGS) -c  $< -o $@

clean: 
	$(RM) *.o

superclean:
	$(RM) *.o *.out $(EXEC)
   
hyperclean:
	$(RM) *.o *.out *.dat $(EXEC)

go : $(EXEC)
	./$(EXEC) < entree.in
