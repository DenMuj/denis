# DBEC-GP-OMP-CUDA-MPI programs are developed by:
#
# Vladimir Loncar, Antun Balaz
# (Scientific Computing Laboratory, Institute of Physics Belgrade, Serbia)
#
# Srdjan Skrbic
# (Department of Mathematics and Informatics, Faculty of Sciences, University of Novi Sad, Serbia)
#
# Paulsamy Muruganandam
# (Bharathidasan University, Tamil Nadu, India)
#
# Luis E. Young-S, Sadhan K. Adhikari
# (UNESP - Sao Paulo State University, Brazil)
#
#
# Public use and modification of these codes are allowed provided that the
# following papers are cited:
# [1] Vladimir Loncar et al., Comput. Phys. Commun. NNN (2016) NNN.
# [2] Vladimir Loncar et al., Comput. Phys. Commun. 200 (2016) 406.
# [3] R. Kishor Kumar et al., Comput. Phys. Commun. 195 (2015) 117.
#
# The authors would be grateful for all information and/or comments
# regarding the use of the programs.

compiler = intel
CC = mpicc
FFTW_INCS=-I$(FFTW3_INCLUDE_PATH)
FFTW_LIBS=-L$(FFTW3_LIBRARY_PATH)
OMPI_INCS=-I$(MPI_INCLUDE)
OMPI_LIBS=-L$(MPI_LIB)
CFLAGS = -O3 -Wall
OMPFLAGS = -fopenmp
LIBS = -lm -lrt -lfftw3_mpi -lfftw3
OMPLIBS = -lfftw3_omp

ifeq ($(compiler), intel)
	CFLAGS = -Ofast -xHOST -Wall
	OMPFLAGS = -qopenmp
endif

MPI_TARGETS=imre3d-ms-psi0 imre3d-ms-ddiX-rot-mpi-qf-grad-mu-self-Q35s imre3d-ddiX-rot-mpi-qf-grad-mu-self-Q35s imre3d-ddi0X-rot-mpi-qf-grad-mu-self-Q35s imre3d-ms-compgX-rot-mpi-qf-mu-self-Q35s imre3d-compgX-rot-mpi-qf-mu-self-Q35s imag3d-compfxv-rot-mpi-qf-mu-self-Q35s imag3d-compgX-rot-mpi-qf-mu-self-Q35s imag3d-compgX-mpi-qf-mu-self-Q35s imag3d-compfxv-cyl-mpi-qf-mu-self-Q35s imag3d-compfxv-mpi-qf-mu-self-Q35s imag3d-comp-mpi-qf-mu-self-Q35s imag3d-compg-mpi-qf-mu-self-Q35s imag3d-mpi-qf-mu-self-Q35s real3d-mpi-qf-mu-nono-Q35s imag3d-mpi imag3d-mpi-qf imag3d-mpi-qf-mu imag3d-mpi-qf-mu-self imag3d-comp-mpi-qf-mu-self imag3d-compg-mpi-qf-mu-self imag3d-mpi-qf-mu-self-Nd-Dd-rnd imag3d-mpi-qf-mu-self-Nd-Dd-rnd2 imag3d-mpi-qf-mu-self-Nd-Dd-rnd3 imag3d-mpi-qf-mu-self-Nd-Dd-rnd4 imag3d-mpi-qf-mu-self-rnd imag3d-mpi-qf-mu-self-Nd imag3d-mpi-qf-mu-self-Nd-Dd imag3d-mpi-qf-mu-self-Q5s imag3d-mpi-qf-mu-self-sin imag3d-mpi-qf-grad-mu-self imag3d-mpi-qf-grad0-mu-self imag3d-mpi-qf-mu-self-link imag3d-mpi-qf-bel-mu-self imag3d-mpi-ring-qf imag3d-mpi-ring-qf-mu imag3d-mpi-ring-qf-mu-self imag3d-mpi-ring-qf-mu-ND-self imag3d-mpi-ring-qf-mu-ND-shift-self imag3d-mpi-ring-qf-mu-ND-shift4-self imag3d-mpi-ring-qf-mu-ND-drop-self imag3d-mpi-ring-qf-mu-self-link imag3d-mpi-ring-qf-mu-prec real3d-mpi real3d-mpi-qf real3d-mpi-qf-mu real3d-mpi-qf-mu-nono real3d-mpi-qf-mu-nono-cont real3d-mpi-qf-mu-nono-K3 real3d-mpi-qf-mu-nono-K3-Q5s real3d-mpi-qf-mu-no-K3-Q5s real3d-mpi-qf-mu-hm real3d-mpi-qf-mu-hm-cont real3d-mpi-ring-qf real3d-mpi-ring-qf-mu
MPI_OBJECTS=$(MPI_TARGETS:=.o)

UTIL_TARGETS=diffint mem cfg
UTIL_OBJECTS=$(UTIL_TARGETS:=.o)

all: $(MPI_TARGETS)
	rm -rf $(UTIL_OBJECTS)

$(MPI_TARGETS): % : $(UTIL_OBJECTS) %.o
	$(CC) $(CFLAGS) $(OMPFLAGS) $(FFTW_INCS) $(OMPI_INCS) -o $@ $@.o $(UTIL_OBJECTS) $(LIBS) $(FFTW_LIBS) $(OMPLIBS)
	rm -rf $@.o

$(MPI_OBJECTS): %.o:
	$(CC) $(CFLAGS) $(OMPFLAGS) $(FFTW_INCS) $(OMPI_INCS) -c src/$*/$*.c -o $@

$(UTIL_OBJECTS): %.o:
	$(CC) $(CFLAGS) $(FFTW_INCS) $(OMPI_INCS) -c src/utils/$*.c -o $@

clean:
	rm -rf $(MPI_TARGETS) $(MPI_OBJECTS) $(UTIL_OBJECTS)

clear:
	rm -rf $(UTIL_OBJECTS)
