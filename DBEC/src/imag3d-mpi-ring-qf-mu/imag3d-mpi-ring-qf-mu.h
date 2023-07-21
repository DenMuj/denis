/**
 * DBEC-GP-OMP-CUDA-MPI programs are developed by:
 *
 * Vladimir Loncar, Antun Balaz
 * (Scientific Computing Laboratory, Institute of Physics Belgrade, Serbia)
 *
 * Srdjan Skrbic
 * (Department of Mathematics and Informatics, Faculty of Sciences, University of Novi Sad, Serbia)
 *
 * Paulsamy Muruganandam
 * (Bharathidasan University, Tamil Nadu, India)
 *
 * Luis E. Young-S, Sadhan K. Adhikari
 * (UNESP - Sao Paulo State University, Brazil)
 *
 *
 * Public use and modification of these codes are allowed provided that the
 * following papers are cited:
 * [1] V. Loncar et al., Comput. Phys. Commun. 209 (2016) 190.
 * [2] V. Loncar et al., Comput. Phys. Commun. 200 (2016) 406.
 * [3] R. Kishor Kumar et al., Comput. Phys. Commun. 195 (2015) 117.
 *
 * The authors would be grateful for all information and/or comments
 * regarding the use of the programs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <omp.h>
#include <sys/time.h>
#include <mpi.h>

#define MAX(a, b, c)       (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c)
#define MAX_FILENAME_SIZE  256
#define RMS_ARRAY_SIZE     4

#define BOHR_RADIUS        5.2917720859e-11

#define FFT_FLAG FFTW_MEASURE
/* #define FFT_FLAG FFTW_ESTIMATE */

static enum OutputFlags { DEN_X = 1 << 0, DEN_Y = 1 << 1, DEN_Z = 1 << 2, DEN_XY = 1 << 3, DEN_XZ = 1 << 4, DEN_YZ = 1 << 5, DEN_XY0 = 1 << 6, DEN_X0Z = 1 << 7, DEN_0YZ = 1 << 8, DEN_XYZ = 1 << 9 } outflags;

char *input, *output, *muoutput, *rmsout, *Niterout, *finalout;
long outstpx, outstpy, outstpz;

double edd, h2, h4, q3, q5, norm_psi2, norm_psi3, ux, uy, uz, ueta, Rring;
int QF;

int opt;
long Na;
long Niter, Nsnap;
long Nx, Ny, Nz;
int rank, nprocs;
long localNx, localNy, offsetNx, offsetNy;
long Nx2, Ny2, Nz2;
double g, gd;
double aho, tau, as, add;
double dx, dy, dz;
double dx2, dy2, dz2;
double dt;
double vnu, vlambda, vgamma, veta;
double par;
double pi, cutoff;

double *x, *y, *z;
double *x2, *y2, *z2;
double ***pot;
double ***potdd;

double Ax0, Ay0, Az0, Ax0r, Ay0r, Az0r, Ax, Ay, Az;
double *calphax, *calphay, *calphaz;
double *cgammax, *cgammay, *cgammaz;

fftw_plan plan_forward, plan_backward, plan_forward0, plan_backward0;
fftw_plan plan_transpose_x, plan_transpose_y, plan_transpose_dpsi;

void readpar(void);
void initpsi(double ***);
void initpot(void);
void gencoef(void);
void initpotdd(double *, double *, double *, double *, double *, double *);
void calcnorm(double *, double ***, double **, double **, double **);
void calcmu(double *, double ***, double ***, double ***, double ***, double ***, fftw_complex *, double **, double **, double **, double **, double **, double **);
void calcmudet(double *, double ***, double ***, double ***, double ***, double ***, double ***, fftw_complex *, double **, double **, double **, double **, double **, double **, double ***, double ***, double ***);
void calcpsidd2(double ***, double ***, fftw_complex *);
void calcpsidd20(double ***, double ***, fftw_complex *);
void calcrms(double *, double ***, double ***, double **, double **, double **);
void calcnu(double ***, double ***, fftw_complex *);
void calclux(double ***, double **);
void calcluy(double ***, double **);
void calcluz(double ***, double **);

void outdenx(double ***, double **, double *, double *, MPI_File);
void outdeny(double ***, double **, double *, double *, MPI_File);
void outdenz(double ***, double **, double *, double *, MPI_File);
void outdenxy(double ***, double ***, double *, MPI_File);
void outdenxz(double ***, double ***, double *, MPI_File);
void outdenyz(double ***, double ***, double *, MPI_File);
void outpsi2xy(double ***, double ***, MPI_File);
void outpsi2xz(double ***, double ***, MPI_File);
void outpsi2yz(double ***, double ***, MPI_File);
void outdenxyz(double ***, double ***, MPI_File);

extern double simpint(double, double *, long);
extern void diff(double, double *, double *, long);

extern int cfg_init(char *);
extern char *cfg_read(char *);

extern double *alloc_double_vector(long);
extern double **alloc_double_matrix(long, long);
extern double ***alloc_double_tensor(long, long, long);
extern fftw_complex *alloc_fftw_complex_vector(long);
extern void free_double_vector(double *);
extern void free_double_matrix(double **);
extern void free_double_tensor(double ***);
extern void free_fftw_complex_vector(fftw_complex *);
