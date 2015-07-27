######################################################################
#	Makefile for MEMSAT3                                         #
######################################################################

#### CC		Command to execute C compiler
#### CFLAGS	Flags to pass to C compiler.

INCDIR		= .
CC		= cc

CFLAGS		= -O
LIBS		= -lm

all:		globmem mem_pred nnsat

install:
		/bin/cp globmem mem_pred nnsat ../bin

clean:
		/bin/rm -f globmem mem_pred nnsat

globmem:	globmem.c globmem_net.h
		$(CC) $(CFLAGS) globmem.c $(LIBS) -o globmem

mem_pred:	mem_pred.c sigmem_net.h
		$(CC) $(CFLAGS) mem_pred.c $(LIBS) -o mem_pred

nnsat:		nnsat.c
		$(CC) $(CFLAGS) nnsat.c $(LIBS) -o nnsat
