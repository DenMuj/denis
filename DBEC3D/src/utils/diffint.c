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

#include "diffint.h"

/**
 *    Spatial 1D integration with Simpson's rule.
 *    h - space step
 *    f - array with the function values
 *    N - number of integration points
 */
double simpint(double h, double *f, long N) {
   /*int c;
   long cnti;
   double sum;

   sum = f[0];
   for (cnti = 1; cnti < N - 1; cnti ++) {
      c = 2 + 2 * (cnti % 2);
      sum += c * f[cnti];
   }
   sum += f[N - 1];

   return sum * h / 3.;*/

   long cnti;
   double sum, sumi, sumj, sumk;

   sumi = 0.; sumj = 0.; sumk = 0.;

   for(cnti = 1; cnti < N - 1; cnti += 2) {
      sumi += f[cnti];
      sumj += f[cnti - 1];
      sumk += f[cnti + 1];
   }

   sum = sumj + 4. * sumi + sumk;
   if(N % 2 == 0) sum += (5. * f[N - 1] + 8. * f[N - 2] - f[N - 3]) / 4.;

   return sum * h / 3.;
}

/**
 *    Richardson extrapolation formula for calculation of space
 *    derivatives of a real function.
 *    h  - space step
 *    f  - array with the function values
 *    df - array with the first derivatives of the function
 *    N  - number of space mesh points
 */
void diff(double h, double *f, double *df, long N) {
   long cnti;

   df[0] = 0.;
   df[1] = (f[2] - f[0]) / (2. * h);

   for (cnti = 2; cnti < N - 2; cnti ++) {
      df[cnti] = (f[cnti - 2] - 8. * f[cnti - 1] + 8. * f[cnti + 1] - f[cnti + 2]) / (12. * h);
   }

   df[N - 2] = (f[N - 1] - f[N - 3]) / (2. * h);
   df[N - 1] = 0.;

   return;
}

/**
 *    Richardson extrapolation formula for calculation of space
 *    derivatives of a complex function..
 *    h  - space step
 *    f  - array with the function values
 *    df - array with the first derivatives of the function
 *    N  - number of space mesh points
 */
void diffc(double h, double complex *f, double complex *df, long N) {
   long cnti;

   df[0] = 0.;
   df[1] = (f[2] - f[0]) / (2. * h);

   for (cnti = 2; cnti < N - 2; cnti ++) {
      df[cnti] = (f[cnti - 2] - 8. * f[cnti - 1] + 8. * f[cnti + 1] - f[cnti + 2]) / (12. * h);
   }

   df[N - 2] = (f[N - 1] - f[N - 3]) / (2. * h);
   df[N - 1] = 0.;

   return;
}
