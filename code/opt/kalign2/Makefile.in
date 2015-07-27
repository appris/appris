PREFIX			= /usr/local/bin
TEST				= test/
CC					= gcc  
CFLAGS			= -O9  -Wall 
DEBUGFLAGS		= -ggdb -Wall 

SOURCES			= kalign2_distance_calculation.c kalign2_dp.c kalign2_input.c kalign2_main.c kalign2_mem.c kalign2_inferface.c kalign2_misc.c kalign2_tree.c  kalign2_profile.c kalign2_alignment_types.c kalign2_feature.c kalign2_hirschberg.c  kalign2_advanced_gaps.c kalign2_hirschberg_dna.c  kalign2_output.c kalign2_string_matching.c kalign2_profile_alignment.c 
PROGS				= kalign
DEBUGPROGS	= kalign_debug


OBJECTS			= $(SOURCES:.c=.o)
DEBUGOBJECTS	= $(SOURCES:.c=_debug.o)


.PHONY: clean

all: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(PROGS)
	
%.o:	%.c
	$(CC) $(CFLAGS) -c $<
	
debug: $(DEBUGOBJECTS)
	$(CC) $(DEBUGFLAGS) $(DEBUGOBJECTS) -o $(DEBUGPROGS)	

%_debug.o:	%.c
	$(CC) $(DEBUGFLAGS) -c $< -o $@


install:
	cp $(PROGS) /usr/local/bin/
	
clean: 
	rm -f $(PROGS) $(OBJECTS)
	rm -f $(DEBUGPROGS) $(DEBUGOBJECTS)
	rm -f config.status config.log config.h Makefile
	rm -f test*
	rm -f *~
	
