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
OPTC=${OPTCLOCAL}

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
OBJLIBPOISSON= lib_poisson1D$(SOL).o lib_poisson1D_writers.o lib_poisson1D_richardson$(SOL).o
OBJTP2ITER= $(OBJLIBPOISSON) tp_poisson1D_iter.o
OBJTP2DIRECT= $(OBJLIBPOISSON) tp_poisson1D_direct.o

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
	$(CC) -o bin/tpPoisson1D_direct $(OPTC) $(addprefix $(TPDIROBJ)/, $(OBJTP2ITER)) $(LIBS)

run_testenv: bin/tp_testenv
	bin/tp_testenv

run_tpPoisson1D_iter: bin/tpPoisson1D_iter
	bin/tpPoisson1D_iter
	bin/tpPoisson1D_iter 1
	bin/tpPoisson1D_iter 2

run_tpPoisson1D_direct: bin/tpPoisson1D_direct
	bin/tpPoisson1D_direct
	bin/tpPoisson1D_direct 1
	bin/tpPoisson1D_direct 2

run_test:
	$(CC) -o bin/test $(OPTC) $(INCL) $(TPDIRSRC)/test.c $(LIBS)
	bin/test

clean:
	rm *.o bin/*
