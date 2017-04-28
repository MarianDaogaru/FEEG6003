# Makefile for the CaseStudy

#
# C compiler and options
#


CC=	mpicc
LIB=	-lm

#
# Object files
#
OBJ=	do_2d.o \
	pgmio.o \
	comms.o

#
# Compile
#
do_2d:	$(OBJ)
	$(CC) -o $@ $(OBJ) $(LIB)

.c.o:
	$(CC) -c $<

#
# Clean out object files and the executable.
#
clean:
	rm *.o do_2d
