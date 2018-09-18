# MAKEFILE src/programs/5th_ed/p121/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p121 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p121 executable as RELEASE edition
#   debug: Build p121 executable as DEBUG edition
#   install: Copy p121 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p121 executable
#   relink: Rebuild p121 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

MERO_DIR = $(HOME)/array_work/
MERO_LD_FLAGS = -L$(MERO_DIR) -lmeroarray

MODULE	= ../../../../include/mpi
LIB	= ../../../../lib
.SUFFIXES: .o .f90

all:	release

release:
	@echo
	@echo "Building P121 RELEASE"
	$(MAKE) p121 \
	FFLAGS="-c -I$(MODULE) $(P121_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P121_MPI_REL_LDFLAGS)"
	@echo "Done p121 RELEASE"
	@echo
	@echo "Building P121_MERO RELEASE"
	$(MAKE) p121_mero \
	FFLAGS="-c -I$(MODULE) $(P121_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) $(MERO_LD_FLAGS) -lParaFEM_mpi.$(VERSION) $(P121_MPI_REL_LDFLAGS)"
	@echo "Done p121_mero RELEASE"

debug:
	@echo
	@echo "Building p121 DEBUG"
	$(MAKE) p121 \
	FFLAGS="-c -I$(MODULE) $(P121_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P121_MPI_DBG_LDFLAGS)"
	@echo "Done p121 DEBUG"

p121: p121.o
	$(FC) p121.o -o p121 $(LDFLAGS)

p121.o: p121.f90
	$(FC) p121.f90 $(FFLAGS)

p121_mero: p121_mero.o
	$(FC) p121_mero.o -o p121_mero $(LDFLAGS) $(MERO_LD_FLAGS)

p121_mero.o: p121_mero.f90
	$(FC) p121_mero.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p121 p121_mero ../../../../bin/p121 ../../../../bin/p121_mero

relink: p121 

install:
	cp p121 ../../../../bin
	cp p121_mero ../../../../bin

test:
	./test.sh p121
