##########################################
# Makefile                               #
# Makefile for the code developed in TP1 #
#                                        #
# T. Dufaud                              #
##########################################
################################
# Variables for this makefile
################################
# 
# -- option to dedicated machine
#
# Rk: You should create a file such as your-machineName.mk
# Follow the example of the machine called "ambre" in the 
# file ambre.mk
#
HOSTNAME?=$(shell hostname)
include $(HOSTNAME).mk

# 
# -- Compiler Option
DEBUG_FLAGS = -Og -g
RELEASE_FLAGS = -O3
OPTC=${OPTCLOCAL} ${RELEASE_FLAGS}

#
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src
TPDIROBJ=$(TPDIR)/bin

#
# -- librairies
LIBS=${LIBSLOCAL}

# -- Include directories
INCLATLAS=${INCLUDEBLASLOCAL}
INCL= -I $(TPDIR)/include $(INCLATLAS)

#
#################################################################
# makefile
############
#
SOL?=
OBJENV= tp_env.o
OBJLIBPOISSON= lib_poisson1D$(SOL).o lib_poisson1D_writers.o lib_poisson1D_richardson$(SOL).o lib_poisson1D_csc_csr.o
OBJTP2ITER= $(OBJLIBPOISSON) tp_poisson1D_iter.o
OBJTP2DIRECT= $(OBJLIBPOISSON) tp_poisson1D_direct.o

TEST_FOLDERS = tests
TEST_FILES = $(wildcard $(TEST_FOLDERS)/*.c)

#
.PHONY: all

all: bin/tp_testenv bin/tpPoisson1D_iter bin/tpPoisson1D_direct
run: run_testenv run_tpPoisson1D_iter run_tpPoisson1D_direct

testenv: bin/tp_testenv

tpPoisson1D_iter: bin/tpPoisson1D_iter

tpPoisson1D_direct: bin/tpPoisson1D_direct

%.o : $(TPDIRSRC)/%.c
	$(CC) $(OPTC) -c $(INCL) $< -o $(TPDIROBJ)/$@

bin/tp_testenv: $(OBJENV)
	$(CC) -o bin/tp_testenv $(OPTC) $(addprefix $(TPDIROBJ)/, $(OBJENV)) $(LIBS)

bin/tpPoisson1D_iter: $(OBJTP2ITER)
	$(CC) -o bin/tpPoisson1D_iter $(OPTC) $(addprefix $(TPDIROBJ)/, $(OBJTP2ITER)) $(LIBS)

bin/tpPoisson1D_direct: $(OBJTP2DIRECT)
	$(CC) -o bin/tpPoisson1D_direct $(OPTC) $(addprefix $(TPDIROBJ)/, $(OBJTP2DIRECT)) $(LIBS)

# Pour compiler tous les tests
bin/%: $(TEST_FOLDERS)/%.c $(OBJLIBPOISSON)
	$(CC) -o $@ $(OPTC) $(INCL) $(addprefix $(TPDIROBJ)/, $(OBJLIBPOISSON)) $< $(LIBS)

bin/benchmark_direct_methods: benchmarks/benchmark_direct_methods.c $(OBJLIBPOISSON)
	$(CC) -o $@ $(OPTC) $(INCL) $(addprefix $(TPDIROBJ)/, $(OBJLIBPOISSON)) $< $(LIBS)

run_tests: bin/test_creation_poisson bin/test_forward_error bin/test_facto_LU bin/test_csr
	$(foreach test, $^, $(test);)

run_testenv: bin/tp_testenv
	bin/tp_testenv

run_benchmark: bin/benchmark_direct_methods
	bin/benchmark_direct_methods

run_tpPoisson1D_iter: bin/tpPoisson1D_iter
	bin/tpPoisson1D_iter
	#bin/tpPoisson1D_iter 1
	#bin/tpPoisson1D_iter 2

run_tpPoisson1D_direct: bin/tpPoisson1D_direct
	bin/tpPoisson1D_direct
	bin/tpPoisson1D_direct 1
	bin/tpPoisson1D_direct 2

clean:
	rm *.o bin/*
