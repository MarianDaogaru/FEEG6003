# Makefile for the CaseStudy

#
# C compiler and options
#


CC=	mpicc
LIB=	-lm

#
# Object files
#
OBJ=	do_serial.o \
	pgmio.o


#
# Compile
#
do_serial:	$(OBJ)
	$(CC) -o $@ $(OBJ) $(LIB)

.c.o:
	$(CC) -c $<

#
# Clean out object files and the executable.
#
clean:
	rm *.o do_serial
