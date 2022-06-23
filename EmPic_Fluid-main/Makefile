#
# Directoriesi
HOMEDIR = $(shell pwd | sed -e 's/\.*//')
RESDIR = $(HOMEDIR)/results
#
# Files
CODE = run.o
#
MPI=mpicxx
MPIOPTS=-O3 
LINKOPTS=$(MPIOPTS)
HYPREFLAGS=$(LINKOPTS) $(LIBS) -lstdc++
FLAGS1=-Wl,--start-group
LIB=$$MKL_LINK/libmkl_intel_lp64.a $$MKL_LINK/libmkl_sequential.a $$MKL_LINK/libmkl_core.a
FLAGS2=-Wl,--end-group
#
#
# Targets
default: Makefile
	@+make $(CODE)

.SUFFIXES:.c
$(CODE): main_MPI.cpp
	$(MPI) main_MPI.cpp -o $(CODE) $(MPIFLAGS)  $(HYPREFLAGS1) 

run: main_MPI.cpp
	$(MPI) main_MPI.cpp -o $(CODE) $(MPIFLAGS)  $(HYPREFLAGS1) $(FLAGS1)  $(LIB) $(FLAGS2)

debug: ex3.c
	$(MPI) ex3.c -g -o $(CODE) $(MPIFLAGS)  $(HYPREFLAGS) $(FLAGS1)  $(LIB) $(FLAGS2)

cleancode:
	rm -f $(CODE)

clean:
	rm -f $(CODE)
	rm -f ./output/*.dat
	rm -f ./output/ion/*.dat
	rm -f ./output/electron/*.dat
	rm -rf ./output/restart
	cd ./output/results/ion; rm -f *.txt
	cd ./output/results/electron; rm -f *.txt
	cd ./output/results/fields; rm -f *.txt

cleandata:
	rm -f ./output/*.dat
	rm -f ./output/ion/*.dat
	rm -f ./output/electron/*.dat
distclean:
	@+make clean


