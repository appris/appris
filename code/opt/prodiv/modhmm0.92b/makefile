COMPILER = gcc

O_FILE_OPTIONS = -w -c
EXEC_FILE_OPTIONS = -w -lm

IO_FUNCS_MULTIALPHA = readhmm_multialpha.o readseqs_multialpha.o readhmm.o
COMP_FUNCS_MULTIALPHA = core_algorithms_multialpha.o training_algorithms_multialpha.o
GEN_FUNCS = std_funcs.o debug_funcs.o std_calculation_funcs.o
GENGETOPT_FUNCS_SEARCH = cmdline_hmmsearch.o
HEADERS = structs.h funcs.h


BASIC_SEARCHING_PROGS = modhmms


all: $(BASIC_SEARCHING_PROGS)



#gengetopt functions
cmdline_hmmsearch.o: cmdline_hmmsearch.c cmdline_hmmsearch.h
	$(COMPILER) $(O_FILE_OPTIONS) cmdline_hmmsearch.c

#general functions
std_funcs.o: std_funcs.c structs.h funcs.h
	$(COMPILER) $(O_FILE_OPTIONS) std_funcs.c

std_calculation_funcs.o: std_calculation_funcs.c structs.h funcs.h
	$(COMPILER) $(O_FILE_OPTIONS) std_calculation_funcs.c

debug_funcs.o: debug_funcs.c $(HEADERS)
	$(COMPILER) $(O_FILE_OPTIONS) debug_funcs.c


#io functions
readhmm.o: readhmm.c $(HEADERS)
	$(COMPILER) $(O_FILE_OPTIONS) readhmm.c

readhmm_multialpha.o: readhmm_multialpha.c $(HEADERS)
	$(COMPILER) $(O_FILE_OPTIONS) readhmm_multialpha.c

readseqs_multialpha.o: readseqs_multialpha.c $(HEADERS)
	$(COMPILER) $(O_FILE_OPTIONS) readseqs_multialpha.c



#computational functions
core_algorithms_multialpha.o: core_algorithms_multialpha.c $(HEADERS)
	$(COMPILER) $(O_FILE_OPTIONS) core_algorithms_multialpha.c

training_algorithms_multialpha.o: training_algorithms_multialpha.c $(HEADERS)
	$(COMPILER) $(O_FILE_OPTIONS) training_algorithms_multialpha.c



#executables
modhmms: hmmsearch.o $(IO_FUNCS_MULTIALPHA) $(COMP_FUNCS_MULTIALPHA) $(GEN_FUNCS) $(GENGETOPT_FUNCS_SEARCH)
	$(COMPILER) $(EXEC_FILE_OPTIONS)  hmmsearch.o $(IO_FUNCS_MULTIALPHA) $(COMP_FUNCS_MULTIALPHA) $(GEN_FUNCS) $(GENGETOPT_FUNCS_SEARCH) -o modhmms
hmmsearch.o: hmmsearch.c $(HEADERS)
	$(COMPILER) $(O_FILE_OPTIONS) hmmsearch.c
