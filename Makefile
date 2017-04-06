# Makefile for mandelbrot area code

#
# C compiler and options for Cray
#
# CC=	gcc-6
# LIB=

#
# C compiler and options for Intel
#
#CC=     icc -O3 -openmp
#LIB=    -lm

#
# C compiler and options for PGI
#
#CC=     pgcc -O3 -mp
#LIB=	-lm

#
# C compiler and options for GNU
#
CC=     gcc-6 -O3 -fopenmp
LIB=	-lm

#
# Object files
#
OBJ=    loops.o

#
# Compile
#
loops:   $(OBJ)
	$(CC) -o $@ $(OBJ) $(LIB)

.c.o:
	$(CC) -c $<

#
# Clean out object files and the executable.
#
clean:
	rm *.o loops
