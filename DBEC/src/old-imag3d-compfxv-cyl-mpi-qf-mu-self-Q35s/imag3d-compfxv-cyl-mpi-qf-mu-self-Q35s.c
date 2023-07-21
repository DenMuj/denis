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

#include "imag3d-compfxv-cyl-mpi-qf-mu-self-Q35s.h"

int main(int argc, char **argv) {
   FILE *out;
   FILE *filerms;
   FILE *muout;
   MPI_File mpifile;
   char filename[MAX_FILENAME_SIZE];
   int nthreads, rankNx2;
   long offsetNx2;
   long cnti, cntj, cntk, snap, nsteps, cnte, cntl;
   double norm, mu[6], mutot, mutotold;
   double *rms;
   double complex **cbeta;
   double complex ***psi, ***psi_t;
   double complex **tmpdpsi1, **tmpdpsi2;
   double ***dpsi, ***dpsi_t;
   double ***psidd2, ***psidd20;
   double **tmpxi, **tmpyi, **tmpzi, **tmpxj, **tmpyj, **tmpzj, ***tmpxmui, ***tmpymui, ***tmpzmui;
   double **outx, **outy, **outz;
   double ***outxy, ***outxz, ***outyz;
   double ***outxyz;
   fftw_complex *psidd2fft;

   int provided;
   MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   pi = 3.14159265358979;

   struct timeval start, stop, iter_start, iter_stop;
   double wall_time, init_time, iter_time;
   iter_time = 0.;
   gettimeofday(&start, NULL);

   if ((argc != 3) || (strcmp(*(argv + 1), "-i") != 0)) {
      fprintf(stderr, "Usage: %s -i <input parameter file> \n", *argv);
      exit(EXIT_FAILURE);
   }

   if (! cfg_init(argv[2])) {
      fprintf(stderr, "Wrong input parameter file.\n");
      exit(EXIT_FAILURE);
   }

   readpar();

   //---------------------------------------------------------------------------

   edd = (4. * pi / 3.) * gd / g;

   if (edd == 0.) {
      q3 = 1.;
      q5 = 1.;
   } else {
      if(edd == 1.) {
         q3 = 3. * sqrt(3.) / 4.;
         q5 = 3. * sqrt(3.) / 2.;
      } else {
           q3 = creal((-(sqrt(3.)*pow(-1. + edd,2)*clog(-3.*(-1. + edd)*edd)) + 2*(sqrt(edd)*(5. + edd)*sqrt(1. + 2.*edd) +
                sqrt(3.)*pow(-1. + edd,2)*log(3.*edd + sqrt(3.)*sqrt(edd)*sqrt(1. + 2.*edd))))/ (16.*sqrt(edd)));

           q5 = creal((6.*sqrt(1. + 2.*edd)*(11. + edd*(4. + 9.*edd)) + (5.*sqrt(3.)*pow(-1. + edd,3)*(clog(-3.*(-1. + edd)*edd) -
                2.*log(3.*edd + sqrt(3.)*sqrt(edd)*sqrt(1. + 2.*edd))))/sqrt(edd))/96.);
      }
   }

   q5 = 1. + 1.5 * edd * edd;
   q3 = 1. + 0.3 * edd * edd;
   q3 *= QF * QDEPL;
   q5 *= QF;

   h2 = 32. * sqrt(pi) * pow(as  * BOHR_RADIUS / aho, 2.5) * pow(Na, 1.5) * (4 * q5 + q3) / 3.;
   h4 = 8. / (3. * sqrt(pi)) * pow(as * BOHR_RADIUS / aho, 1.5) * sqrt(Na) * q3;

   //---------------------------------------------------------------------------

   assert(Nx % nprocs == 0);
   assert(Ny % nprocs == 0);

   localNx = Nx / nprocs;
   localNy = Ny / nprocs;
   offsetNx = rank * localNx;
   offsetNy = rank * localNy;

   Nx2 = Nx / 2; Ny2 = Ny / 2; Nz2 = Nz / 2;
   dx2 = dx * dx; dy2 = dy * dy; dz2 = dz * dz;

   rankNx2 = Nx2 / localNx;
   offsetNx2 = Nx2 % localNx;

   #pragma omp parallel
   #pragma omp master
   nthreads = omp_get_num_threads();

   // Allocation of memory ------------------------------------------

   rms = alloc_double_vector(RMS_ARRAY_SIZE);

   x = alloc_double_vector(Nx);
   y = alloc_double_vector(Ny);
   z = alloc_double_vector(Nz);

   x2 = alloc_double_vector(Nx);
   y2 = alloc_double_vector(Ny);
   z2 = alloc_double_vector(Nz);

   pot = alloc_double_tensor(localNx, Ny, Nz);
   potdd = alloc_double_tensor(localNy, Nx, Nz);
   psi = alloc_complex_tensor(localNx, Ny, Nz);
   psi_t = alloc_complex_tensor(Nx, localNy, Nz);
   psidd2 = alloc_double_tensor(localNx, Ny, 2 * (Nz2 + 1));
   if (muoutput != NULL) psidd20 = alloc_double_tensor(localNx, Ny, 2 * (Nz2 + 1));

   dpsi = alloc_double_tensor(localNx, Ny, Nz);
   dpsi_t = alloc_double_tensor(Nx, localNy, Nz);

   calphax = alloc_double_vector(Nx - 1);
   calphay = alloc_double_vector(Ny - 1);
   calphaz = alloc_double_vector(Nz - 1);
   cbeta = alloc_complex_matrix(nthreads, MAX(Nx, Ny, Nz) - 1);
   cgammax = alloc_double_vector(Nx - 1);
   cgammay = alloc_double_vector(Ny - 1);
   cgammaz = alloc_double_vector(Nz - 1);

   tmpxi = alloc_double_matrix(nthreads, Nx);
   tmpyi = alloc_double_matrix(nthreads, Ny);
   tmpzi = alloc_double_matrix(nthreads, Nz);
   tmpxj = alloc_double_matrix(nthreads, Nx);
   tmpyj = alloc_double_matrix(nthreads, Ny);
   tmpzj = alloc_double_matrix(nthreads, Nz);
   if (muoutput != NULL) tmpxmui = alloc_double_tensor(6, nthreads, Nx);
   if (muoutput != NULL) tmpymui = alloc_double_tensor(6, nthreads, Ny);
   if (muoutput != NULL) tmpzmui = alloc_double_tensor(6, nthreads, Nz);
   tmpdpsi1 = alloc_complex_matrix(nthreads, MAX(Nx, Ny, Nz));
   tmpdpsi2 = alloc_complex_matrix(nthreads, MAX(Nx, Ny, Nz));


   psidd2fft = alloc_fftw_complex_vector(Nx * localNy * (Nz2 + 1));

   outx = alloc_double_matrix(localNx, 2);
   outy = alloc_double_matrix(localNy, 2);
   outz = alloc_double_matrix(Nz, 2); // Because rank 0 will assemble outz
   outxy = alloc_double_tensor(localNx, Ny, 3);
   outxz = alloc_double_tensor(localNx, Nz, 3);
   outyz = alloc_double_tensor((rank == rankNx2) ? Ny : localNy, Nz, 3);
   outxyz = dpsi;

   // -----------------------------------------------------------------------

   if (opt == 2) par = 2.;
   else par = 1.;

   g = par * g;
   gd = par * gd;
   h2 = par * h2;
   tau = par * tau;

   if (input == NULL) {
      if (urho == 0.) {
         if (vgamma > 0.) {
            urho = 1. / sqrt(vgamma);
         } else urho = 1.;
      }

      if (uz == 0.) {
         if (vlambda > 0.) {
            uz = 1. / sqrt(vlambda);
         } else uz = 1.;
      }
   }

   if (rank == 0) {
      if (output != NULL) {
         sprintf(filename, "%s.txt", output);
         out = fopen(filename, "w");
      } else out = stdout;
   } else out = fopen("/dev/null", "w");

   if (rank == 0) {
      if (muoutput != NULL) {
         sprintf(filename, "%s.txt", muoutput);
         muout = fopen(filename, "w");
      } else muout = NULL;
   } else muout = fopen("/dev/null", "w");

   if (rank == 0) {
      if (rmsout != NULL) {
         sprintf(filename, "%s.txt", rmsout);
         filerms = fopen(filename, "w");
      } else filerms = NULL;
   } else filerms = fopen("/dev/null", "w");

   fprintf(out, "\n**********************************************\n");
   fprintf(out, "Imaginary-time propagation in 3D, OPTION = %d\nMPI nodes = %d, OMP threads = %d, cores = %d\n", opt, nprocs, nthreads, nprocs * nthreads);
   fprintf(out, "**********************************************\n\nInteractions\n");

   if (cfg_read("G") != NULL) {
      fprintf(out, "Contact: G = %.6le\n", g);
   } else {
      fprintf(out, "Contact: Natoms = %.6le, as = %.6le * a0, G = %.6le\n", (double) Na, as, g);
   }

   if (cfg_read("GDD") != NULL) {
      fprintf(out, "DDI: GD = %.6le, edd = %.6le\n", gd, edd);
   } else {
        fprintf(out, "DDI: add = %.6le * a0, GD = %.6le, edd = %.6le\n", add, gd, edd);
   }

   fprintf(out, "     Dipolar cutoff R = %.6le\n\n",  cutoff);

   if (QF == 1) {
      fprintf(out, "QF = 1, QDEPL = %i: h2 = %.6le, h4 = %.6le\n        q3 = %.6le, q5 = %.6le\n\n", QDEPL, h2, h4, q3, q5);
   } else  fprintf(out, "QF = 0\n\n");

   fprintf(out, "Trap parameters\nGAMMA = %.6le, NU = %.6le, LAMBDA = %.6le\n\n", vgamma, vnu, vlambda);
   fprintf(out, "Space discretization: NX = %li, NY = %li, NZ = %li\n", Nx, Ny, Nz);
   fprintf(out, "                      DX = %.6le, DY = %.6le, DZ = %.6le\n", dx, dy, dz);
   if (cfg_read("AHO") != NULL) fprintf(out, "      Unit of length: aho = %.6le m\n", aho);
   fprintf(out, "\nTime discretization: NITER = %li (NSNAP = %li)\n", Niter, Nsnap);
   fprintf(out, "                     DT = %.6le\n",  dt);
   if (cfg_read("TAU") != NULL) fprintf(out, "       Unit of time: tau = %.6le s\n", tau);
   fprintf(out, "\nInitial state: ");
   if (input != NULL) {
      fprintf(out, "file %s\n\n", input);
   } else {
      fprintf(out, "Gaussian\n               URHO = %.6le, UZ = %.6le\n\n", urho, uz);
   }

   fprintf(out, "-------------------------------------------------------------------\n");
   fprintf(out, "Snap      Norm0          Norm           mu             <r>\n");
   fprintf(out, "-------------------------------------------------------------------\n");
   fflush(out);

   if (muoutput != NULL) {
      fprintf(muout, "----------------------------------------------------------------------------------------------------------------\n");
      fprintf(muout, "Snap      mu             Kin            Pot            Contact        DDI            Dmu(Contact)   Dmu(DDI) \n");
      fprintf(muout, "----------------------------------------------------------------------------------------------------------------\n");
      fflush(muout);
   }

   if (rmsout != NULL) {
      fprintf(filerms, "\n**********************************************\n");
      fprintf(filerms, "Imaginary-time propagation in 3D, OPTION = %d\nMPI nodes = %d, OMP threads = %d, cores = %d\n", opt, nprocs, nthreads, nprocs * nthreads);
      fprintf(filerms, "**********************************************\n\nInteractions\n");

      if (cfg_read("G") != NULL) {
         fprintf(filerms, "Contact: G = %.6le\n", g);
      } else {
         fprintf(filerms, "Contact: Natoms = %.6le, as = %.6le * a0, G = %.6le\n", (double) Na, as, g);
      }

      if (cfg_read("GDD") != NULL) {
         fprintf(filerms, "DDI: GD = %.6le, edd = %.6le\n", gd, edd);
      } else {
           fprintf(filerms, "DDI: add = %.6le * a0, GD = %.6le, edd = %.6le\n", add, gd, edd);
      }

      fprintf(filerms, "     Dipolar cutoff R = %.6le\n\n",  cutoff);

      if (QF == 1) {
         fprintf(filerms, "QF = 1: h2 = %.6le, h4 = %.6le\n        q3 = %.6le, q5 = %.6le\n\n", h2, h4, q3, q5);
      } else  fprintf(filerms, "QF = 0\n\n");

      fprintf(filerms, "Trap parameters\nGAMMA = %.6le, NU = %.6le, LAMBDA = %.6le\n\n", vgamma, vnu, vlambda);       
      fprintf(filerms, "Space discretization: NX = %li, NY = %li, NZ = %li\n", Nx, Ny, Nz);
      fprintf(filerms, "                      DX = %.6le, DY = %.6le, DZ = %.6le\n", dx, dy, dz);
      if (cfg_read("AHO") != NULL) fprintf(filerms, "      Unit of length: aho = %.6le m\n", aho);
      fprintf(filerms, "\nTime discretization: NITER = %li (NSNAP = %li)\n", Niter, Nsnap);
      fprintf(filerms, "                     DT = %.6le\n",  dt);
      if (cfg_read("TAU") != NULL) fprintf(filerms, "       Unit of time: tau = %.6le s\n", tau);
      fprintf(filerms, "\nInitial state: ");
      if (input != NULL) {
         fprintf(filerms, "file %s\n\n", input);
      } else {
        fprintf(filerms, "Gaussian\n               URHO = %.6le, UZ = %.6le\n\n", urho, uz);
      }

      fprintf(filerms, "-------------------------------------------------------------------\n");
      fprintf(filerms, "Snap      <r>            <x>            <y>            <z>\n");
      fprintf(filerms, "-------------------------------------------------------------------\n");
      fflush(filerms);
   }

   fftw_init_threads();
   fftw_mpi_init();
   fftw_plan_with_nthreads(nthreads);

   long n[] = { Nx, Ny, Nz };
   plan_forward = fftw_mpi_plan_many_dft_r2c(3, n, 1, localNx, localNy, psidd2[0][0], psidd2fft, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_OUT);
   plan_backward = fftw_mpi_plan_many_dft_c2r(3, n, 1, localNy, localNx, psidd2fft, psidd2[0][0], MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_IN);

   if (muoutput != NULL) plan_forward0 = fftw_mpi_plan_many_dft_r2c(3, n, 1, localNx, localNy, psidd20[0][0], psidd2fft, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_OUT);
   if (muoutput != NULL) plan_backward0 = fftw_mpi_plan_many_dft_c2r(3, n, 1, localNy, localNx, psidd2fft, psidd20[0][0], MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_IN);

   plan_transpose_x = fftw_mpi_plan_many_transpose(Nx, Ny * Nz, 2, localNx, localNy * Nz, (double *) **psi, (double *) **psi_t, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_OUT);
   plan_transpose_y = fftw_mpi_plan_many_transpose(Ny * Nz, Nx, 2, localNy * Nz, localNx, (double *) **psi_t, (double *) **psi, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_IN);
   plan_transpose_dpsi = fftw_mpi_plan_many_transpose(Ny * Nz, Nx, 1, localNy * Nz, localNx, (double *) **dpsi_t, **dpsi, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_IN);

   initpsi(psi);
   initpot();
   gencoef();
   initpotdd(*tmpxi, *tmpyi, *tmpzi, *tmpxj, *tmpyj, *tmpzj);

   calcnorm(&norm, psi, tmpxi, tmpyi, tmpzi);
   norm_psi2 *= norm * norm;
   norm_psi3 *= norm * norm * norm;
   calcrms(rms, psi, psi_t, tmpxi, tmpyi, tmpzi);

   if (muoutput != NULL) {
      calcmudet(mu, psi, psi_t, dpsi, dpsi_t, psidd2, psidd20, psidd2fft, tmpxi, tmpyi, tmpzi, tmpxj, tmpyj, tmpzj, tmpdpsi1, tmpdpsi2, tmpxmui, tmpymui, tmpzmui);

      mutot = 0.;
      for (cnte = 0; cnte < 6; cnte ++) {
         mu[cnte] /= par * norm_psi2;
         mutot += mu[cnte];
      }

      fprintf(muout, "%-9d %-14.6le %-14.6le %-14.6le %-14.6le %-14.6le %-14.6le %-14.6le\n", 0, mutot, mu[0], mu[1], mu[2], mu[3], mu[4], mu[5]);
      fflush(muout);
   } else {
      calcmu(&mutot, psi, psi_t, dpsi, dpsi_t, psidd2, psidd2fft, tmpxi, tmpyi, tmpzi, tmpxj, tmpyj, tmpzj, tmpdpsi1, tmpdpsi2);
      mutot /= par * norm_psi2;
   }

   MPI_Bcast(&mutot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   fprintf(out, "%-9d %-14.6le %-14.6le %-14.6le %-14.6le\n", 0, norm_psi2, norm_psi2 + h4 * norm_psi3, mutot, *rms / sqrt(norm_psi2));
   fflush(out);

   mutotold = mutot;

   if (rmsout != NULL) {
      fprintf(filerms, "%-9d %-14.6le %-14.6le %-14.6le %-14.6le\n", 0, rms[0] / sqrt(norm_psi2), rms[1] / sqrt(norm_psi2), rms[2] / sqrt(norm_psi2), rms[3] / sqrt(norm_psi2));
      fflush(filerms);
   }

   if (Niterout != NULL) {

      char itername[10];
      sprintf(itername, "-%06d-", 0);

      if (outflags & DEN_XYZ) {
         sprintf(filename, "%s%s3d.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outdenxyz(psi, outxyz, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & DEN_X) {
         sprintf(filename, "%s%s1d_x.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outdenx(psi, outx, *tmpyi, *tmpzi, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & DEN_Y) {
         sprintf(filename, "%s%s1d_y.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outdeny(psi_t, outy, *tmpxi, *tmpzi, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & DEN_Z) {
         sprintf(filename, "%s%s1d_z.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outdenz(psi, outz, *tmpxi, *tmpyi, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & DEN_XY) {
         sprintf(filename, "%s%s2d_xy.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outdenxy(psi, outxy, *tmpzi, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & DEN_XZ) {
         sprintf(filename, "%s%s2d_xz.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outdenxz(psi, outxz, *tmpyi, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & DEN_YZ) {
          sprintf(filename, "%s%s2d_yz.bin", Niterout, itername);
          MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
          outdenyz(psi_t, outyz, *tmpxi, mpifile);
          MPI_File_close(&mpifile);
      }

      if (outflags & DEN_XY0) {
         sprintf(filename, "%s%s3d_xy0.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outpsi2xy(psi, outxy, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & DEN_X0Z) {
         sprintf(filename, "%s%s3d_x0z.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outpsi2xz(psi, outxz, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & DEN_0YZ) {
         sprintf(filename, "%s%s3d_0yz.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outpsi2yz(psi, outyz, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & ARG_XY0) {
         sprintf(filename, "%s-phase%s3d_xy0.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outargxy(psi, outxy, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & ARG_X0Z) {
         sprintf(filename, "%s-phase%s3d_x0z.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outargxz(psi, outxz, mpifile);
         MPI_File_close(&mpifile);
      }

      if (outflags & ARG_0YZ) {
         sprintf(filename, "%s-phase%s3d_0yz.bin", Niterout, itername);
         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
         outargyz(psi, outyz, mpifile);
         MPI_File_close(&mpifile);
      }
   }

   //  -----------------------------------------------------------NITER

   gettimeofday(&iter_start, NULL);

   nsteps = Niter / Nsnap;
   for (snap = 1; snap <= Nsnap; snap ++) {
      for (cntl = 0; cntl < nsteps; cntl ++) {
         calcnu(psi, psidd2, psidd2fft);
         calcluy(psi, cbeta);
         calcluz(psi, cbeta);

         #pragma omp for private(cnti, cntj, cntk)
         for (cnti = 0; cnti < localNx; cnti ++) {
            for (cntj = 0; cntj < Ny; cntj ++) {
               for (cntk = 0; cntk < cntj; cntk ++) {
                  psi[cnti][cntj][cntk] = 0.5 * (psi[cnti][cntj][cntk] + psi[cnti][cntk][cntj]);
                  psi[cnti][cntk][cntj] = psi[cnti][cntj][cntk];
               }
            }
         }

         calclux(psi_t, cbeta);
         calcnorm(&norm, psi, tmpxi, tmpyi, tmpzi);
      }

      norm_psi2 *= norm * norm;
      norm_psi3 *= norm * norm * norm;
      calcrms(rms, psi, psi_t, tmpxi, tmpyi, tmpzi);

      if (muoutput != NULL) {
         calcmudet(mu, psi, psi_t, dpsi, dpsi_t, psidd2, psidd20, psidd2fft, tmpxi, tmpyi, tmpzi, tmpxj, tmpyj, tmpzj, tmpdpsi1, tmpdpsi2, tmpxmui, tmpymui, tmpzmui);

         mutot = 0.;
         for (cnte = 0; cnte < 6; cnte ++) {
            mu[cnte] /= par * norm_psi2;
            mutot += mu[cnte];
         }

         fprintf(muout, "%-9li %-14.6le %-14.6le %-14.6le %-14.6le %-14.6le %-14.6le %-14.6le\n", snap, mutot, mu[0], mu[1], mu[2], mu[3], mu[4], mu[5]);
         fflush(muout);
      } else {
         calcmu(&mutot, psi, psi_t, dpsi, dpsi_t, psidd2, psidd2fft, tmpxi, tmpyi, tmpzi, tmpxj, tmpyj, tmpzj, tmpdpsi1, tmpdpsi2);
         mutot /= par * norm_psi2;
      }

      MPI_Bcast(&mutot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      fprintf(out, "%-9li %-14.6le %-14.6le %-14.6le %-14.6le\n", snap, norm_psi2, norm_psi2 + h4 * norm_psi3, mutot, *rms / sqrt(norm_psi2));
      fflush(out);

      if (rmsout != NULL) {
         fprintf(filerms, "%-9li %-14.6le %-14.6le %-14.6le %-14.6le\n", snap, rms[0] / sqrt(norm_psi2), rms[1] / sqrt(norm_psi2), rms[2] / sqrt(norm_psi2), rms[3] / sqrt(norm_psi2));
         fflush(filerms);
      }

      if (Niterout != NULL) {

         char itername[10];
         sprintf(itername, "-%06li-", snap);

         if (outflags & DEN_XYZ) {
            sprintf(filename, "%s%s3d.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenxyz(psi, outxyz, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & DEN_X) {
            sprintf(filename, "%s%s1d_x.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenx(psi, outx, *tmpyi, *tmpzi, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & DEN_Y) {
            sprintf(filename, "%s%s1d_y.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdeny(psi_t, outy, *tmpxi, *tmpzi, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & DEN_Z) {
            sprintf(filename, "%s%s1d_z.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenz(psi, outz, *tmpxi, *tmpyi, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & DEN_XY) {
            sprintf(filename, "%s%s2d_xy.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenxy(psi, outxy, *tmpzi, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & DEN_XZ) {
            sprintf(filename, "%s%s2d_xz.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenxz(psi, outxz, *tmpyi, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & DEN_YZ) {
             sprintf(filename, "%s%s2d_yz.bin", Niterout, itername);
             MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
             outdenyz(psi_t, outyz, *tmpxi, mpifile);
             MPI_File_close(&mpifile);
         }

         if (outflags & DEN_XY0) {
            sprintf(filename, "%s%s3d_xy0.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outpsi2xy(psi, outxy, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & DEN_X0Z) {
            sprintf(filename, "%s%s3d_x0z.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outpsi2xz(psi, outxz, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & DEN_0YZ) {
            sprintf(filename, "%s%s3d_0yz.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outpsi2yz(psi, outyz, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & ARG_XY0) {
            sprintf(filename, "%s-phase%s3d_xy0.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outargxy(psi, outxy, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & ARG_X0Z) {
            sprintf(filename, "%s-phase%s3d_x0z.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outargxz(psi, outxz, mpifile);
            MPI_File_close(&mpifile);
         }

         if (outflags & ARG_0YZ) {
            sprintf(filename, "%s-phase%s3d_0yz.bin", Niterout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outargyz(psi, outyz, mpifile);
            MPI_File_close(&mpifile);
         }
      }

      if (fabs((mutotold - mutot) / mutot) < murel) break;
      mutotold = mutot;
   }

   gettimeofday(&iter_stop, NULL);
   iter_time += (double) (((iter_stop.tv_sec - iter_start.tv_sec) * 1000 + (iter_stop.tv_usec - iter_start.tv_usec)/1000.0) + 0.5);

   fprintf(out, "-------------------------------------------------------------------\n");
   fflush(out);
   fprintf(muout, "----------------------------------------------------------------------------------------------------------------\n\n");
   fflush(muout);

   if (rmsout != NULL) {
      fprintf(filerms, "-------------------------------------------------------------------\n\n");
      fclose(filerms);
   }

   if (finalpsi != NULL) {
       sprintf(filename, "%s.bin", finalpsi);
       MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
       outpsi(psi, mpifile);
       MPI_File_close(&mpifile);
   }

   // Free all dynamically allocated memory. ---------------

   free_double_vector(rms);

   free_double_vector(x);
   free_double_vector(y);
   free_double_vector(z);

   free_double_vector(x2);
   free_double_vector(y2);
   free_double_vector(z2);

   free_double_tensor(pot);
   free_double_tensor(potdd);
   free_complex_tensor(psi);
   free_complex_tensor(psi_t);
   free_double_tensor(psidd2);
   if (muoutput != NULL) free_double_tensor(psidd20);

   free_double_tensor(dpsi);
   free_double_tensor(dpsi_t);

   free_double_vector(calphax);
   free_double_vector(calphay);
   free_double_vector(calphaz);
   free_complex_matrix(cbeta);
   free_double_vector(cgammax);
   free_double_vector(cgammay);
   free_double_vector(cgammaz);

   free_double_matrix(tmpxi);
   free_double_matrix(tmpyi);
   free_double_matrix(tmpzi);
   free_double_matrix(tmpxj);
   free_double_matrix(tmpyj);
   free_double_matrix(tmpzj);
   if (muoutput != NULL) free_double_tensor(tmpxmui);
   if (muoutput != NULL) free_double_tensor(tmpymui);
   if (muoutput != NULL) free_double_tensor(tmpzmui);
   free_complex_matrix(tmpdpsi1);
   free_complex_matrix(tmpdpsi2);

   fftw_destroy_plan(plan_forward);
   fftw_destroy_plan(plan_backward);
   if (muoutput != NULL) fftw_destroy_plan(plan_forward0);
   if (muoutput != NULL) fftw_destroy_plan(plan_backward0);
   fftw_destroy_plan(plan_transpose_x);
   fftw_destroy_plan(plan_transpose_y);
   fftw_destroy_plan(plan_transpose_dpsi);

   free_fftw_complex_vector(psidd2fft);

   free_double_matrix(outx);
   free_double_matrix(outy);
   free_double_matrix(outz);
   free_double_tensor(outxy);
   free_double_tensor(outxz);
   free_double_tensor(outyz);

   fftw_mpi_cleanup();

   // ----------------------------------------------------

   MPI_Finalize();

   gettimeofday(&stop, NULL);
   wall_time = (double) (((stop.tv_sec - start.tv_sec) * 1000 + (stop.tv_usec - start.tv_usec)/1000.0) + 0.5);
   init_time = wall_time - iter_time;
   fprintf(out, "\nInitialization/allocation wall-clock time: %.3f seconds\n", init_time / 1000.);
   fprintf(out, "Calculation (iterations) wall-clock time : %.3f seconds\n\n", iter_time / 1000.);

   if(output != NULL) fclose(out);
   if(muoutput != NULL) fclose(muout);

   return (EXIT_SUCCESS);
}

/**
 *    Reading input parameters from the configuration file.
 */
void readpar(void) {
   char *cfg_tmp;

   if ((cfg_tmp = cfg_read("OPTION")) == NULL) {
      fprintf(stderr, "OPTION is not defined in the configuration file\n");
      exit(EXIT_FAILURE);
   }
   opt = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("NATOMS")) == NULL) {
      fprintf(stderr, "NATOMS is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Na = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("AHO")) == NULL) {
      fprintf(stderr, "AHO is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   aho = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("G")) == NULL) {
      if ((cfg_tmp = cfg_read("AS")) == NULL) {
         fprintf(stderr, "AS is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      as = atof(cfg_tmp);

      g = 4. * pi * as * Na * BOHR_RADIUS / aho;
   } else {
      g = atof(cfg_tmp);
   }

   if ((cfg_tmp = cfg_read("GDD")) == NULL) {
      if ((cfg_tmp = cfg_read("ADD")) == NULL) {
         fprintf(stderr, "ADD is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      add = atof(cfg_tmp);

      gd = 3. * add * Na * BOHR_RADIUS / aho;
   } else {
      gd = atof(cfg_tmp);
   }

   if ((cfg_tmp = cfg_read("TAU")) != NULL) tau = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("QF")) == NULL) {
      QF = 0;
   } else QF = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("QDEPL")) == NULL) {
      QDEPL = 0;
   } else QDEPL = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("NX")) == NULL) {
      fprintf(stderr, "NX is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Nx = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("NY")) == NULL) {
      fprintf(stderr, "NY is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Ny = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("NZ")) == NULL) {
      fprintf(stderr, "Nz is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Nz = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("DX")) == NULL) {
      fprintf(stderr, "DX is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dx = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("DY")) == NULL) {
      fprintf(stderr, "DY is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dy = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("DZ")) == NULL) {
      fprintf(stderr, "DZ is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dz = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("DT")) == NULL) {
      fprintf(stderr, "DT is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dt = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("MUREL")) == NULL) {
      fprintf(stderr, "MUREL is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   murel = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("GAMMA")) == NULL) {
      fprintf(stderr, "GAMMA is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   vgamma = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("NU")) == NULL) {
      fprintf(stderr, "NU is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   vnu = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("LAMBDA")) == NULL) {
      fprintf(stderr, "LAMBDA is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   vlambda = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("NITER")) == NULL) {
      fprintf(stderr, "NITER is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Niter = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("NSNAP")) == NULL) Nsnap = 1;
   else Nsnap = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("CUTOFF")) == NULL) {
      fprintf(stderr, "CUTOFF is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   cutoff = atof(cfg_tmp);

   input = cfg_read("INPUT");
   input_type = cfg_read("INPUT_TYPE");
   output = cfg_read("OUTPUT");
   muoutput = cfg_read("MUOUTPUT");
   rmsout = cfg_read("RMSOUT");
   Niterout = cfg_read("NITEROUT");
   finalpsi = cfg_read("FINALPSI");

   if ((cfg_tmp = cfg_read("URHO")) == NULL) urho = 0.;
   else urho = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("UZ")) == NULL) uz = 0.;
   else uz = atof(cfg_tmp);

   if ((cfg_tmp = cfg_read("XI")) == NULL) xi = 0.;
   else xi = atof(cfg_tmp);

   if (Niterout != NULL) {
      if ((cfg_tmp = cfg_read("OUTFLAGS")) == NULL) {
         fprintf(stderr, "OUTFLAGS is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outflags = atoi(cfg_tmp);           
   } else outflags = 0;
 
   if ((Niterout != NULL) || (finalpsi != NULL)) {
      if ((cfg_tmp = cfg_read("OUTSTPX")) == NULL) {
         fprintf(stderr, "OUTSTPX is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpx = atol(cfg_tmp);

      if ((cfg_tmp = cfg_read("OUTSTPY")) == NULL) {
         fprintf(stderr, "OUTSTPY is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpy = atol(cfg_tmp);

      if ((cfg_tmp = cfg_read("OUTSTPZ")) == NULL) {
         fprintf(stderr, "OUTSTPZ is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpz = atol(cfg_tmp);
   }

   return;
}

/**
 *    Initialization of the space mesh and the initial wave function.
 *    psi - array with the wave function values
 *    tmpz - temporary array
 */
void initpsi(double complex ***psi) {
   long cnti, cntj, cntk;
   double cpsi;
   double tmp, *tmpr, *tmpc;
   MPI_Offset fileoffset;
   MPI_File file = MPI_FILE_NULL;

   tmpr = alloc_double_vector(Nz);
   tmpc = alloc_double_vector(2 * Nz);

   #pragma omp parallel for private(cnti)
   for (cnti = 0; cnti < Nx; cnti ++) {
      x[cnti] = (cnti - Nx2) * dx;
      x2[cnti] = x[cnti] * x[cnti];
   }

   #pragma omp parallel for private(cntj)
   for (cntj = 0; cntj < Ny; cntj ++) {
      y[cntj] = (cntj - Ny2) * dy;
      y2[cntj] = y[cntj] * y[cntj];
   }

   #pragma omp parallel for private(cntk)
   for (cntk = 0; cntk < Nz; cntk ++) {
      z[cntk] = (cntk - Nz2) * dz;
      z2[cntk] = z[cntk] * z[cntk];
   }

   if (input != NULL) {
      MPI_File_open(MPI_COMM_WORLD, input, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

      if (file == MPI_FILE_NULL) {
         fprintf(stderr, "Specify the proper input file with the initial wave function.\n");
         MPI_Finalize();
         exit(EXIT_FAILURE);
      }

      if (strcmp(input_type, "DEN") == 0) {
         fileoffset = rank * sizeof(double) * localNx * Ny * Nz;
         for (cnti = 0; cnti < localNx; cnti ++) {
            for (cntj = 0; cntj < Ny; cntj ++) {
               MPI_File_read_at(file, fileoffset, tmpr, Nz, MPI_DOUBLE, MPI_STATUS_IGNORE);
               for (cntk = 0; cntk < Nz; cntk ++) {
                  psi[cnti][cntj][cntk] = sqrt(tmpr[cntk]);
               }
               fileoffset += Nz * sizeof(double);
            }
         }
      } else {
         if (strcmp(input_type, "PSI") == 0) {
            fileoffset = rank * sizeof(double complex) * localNx * Ny * Nz;
            for (cnti = 0; cnti < localNx; cnti ++) {
               for (cntj = 0; cntj < Ny; cntj ++) {
                  MPI_File_read_at(file, fileoffset, tmpc, 2 * Nz, MPI_DOUBLE, MPI_STATUS_IGNORE);
                  for (cntk = 0; cntk < Nz; cntk ++) {
                     psi[cnti][cntj][cntk] = tmpc[2 * cntk] + I * tmpc[2 * cntk + 1];
                  }
                  fileoffset += Nz * sizeof(double complex);
               }
            }
         } else {
            fprintf(stderr, "Specify the proper input_type for the file with the initial wave function.\n");
            MPI_Finalize();
            exit(EXIT_FAILURE);
         }
      }

      MPI_File_close(&file);
      MPI_File_close(&file);
   } else {

      cpsi = sqrt(32.) / (pow(pi, 0.75) * urho * urho * sqrt(uz));

      #pragma omp parallel for private(cnti, cntj, cntk, tmp)
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               tmp = exp(- 2.0 * ((z2[cntk] + y2[cntj]) / (urho * urho) + x2[offsetNx + cnti] / (uz * uz)));
               psi[cnti][cntj][cntk] = cpsi * (y[cntj] + I * z[cntk]) * tmp;
            }
         }
      }
   }

   #pragma omp parallel for private(cnti, cntj, cntk)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            psi[cnti][cntj][cntk] *= (y[cntj] + I * z[cntk]) / sqrt(y2[cntj] + z2[cntk] + xi * xi);
         }
      }
   }

   free_double_vector(tmpr);
   free_double_vector(tmpc);

   return;
}

/**
 *    Initialization of the potential.
 */
void initpot() {
   long cnti, cntj, cntk;
   double vgamma2, vnu2, vlambda2;

   vgamma2 = vgamma * vgamma;
   vnu2 = vnu * vnu;
   vlambda2 = vlambda * vlambda;

   #pragma omp parallel for private(cnti, cntj, cntk)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            pot[cnti][cntj][cntk] = 0.5 * par * (vgamma2 * x2[offsetNx + cnti] + vnu2 * y2[cntj] + vlambda2 * z2[cntk]);
         }
      }
   }

   return;
}

/**
 *    Crank-Nicolson scheme coefficients generation.
 */
void gencoef(void) {
   long cnti;

   Ax0 = 1. + dt / dx2 / (3. - par);
   Ay0 = 1. + dt / dy2 / (3. - par);
   Az0 = 1. + dt / dz2 / (3. - par);

   Ax0r = 1. - dt / dx2 / (3. - par);
   Ay0r = 1. - dt / dy2 / (3. - par);
   Az0r = 1. - dt / dz2 / (3. - par);

   Ax = - 0.5 * dt / dx2 / (3. - par);
   Ay = - 0.5 * dt / dy2 / (3. - par);
   Az = - 0.5 * dt / dz2 / (3. - par);

   calphax[Nx - 2] = 0.;
   cgammax[Nx - 2] = - 1. / Ax0;
   for (cnti = Nx - 2; cnti > 0; cnti --) {
      calphax[cnti - 1] = Ax * cgammax[cnti];
      cgammax[cnti - 1] = - 1. / (Ax0 + Ax * calphax[cnti - 1]);
   }

   calphay[Ny - 2] = 0.;
   cgammay[Ny - 2] = - 1. / Ay0;
   for (cnti = Ny - 2; cnti > 0; cnti --) {
      calphay[cnti - 1] = Ay * cgammay[cnti];
      cgammay[cnti - 1] = - 1. / (Ay0 + Ay * calphay[cnti - 1]);
   }

   calphaz[Nz - 2] = 0.;
   cgammaz[Nz - 2] = - 1. / Az0;
   for (cnti = Nz - 2; cnti > 0; cnti --) {
      calphaz[cnti - 1] = Az * cgammaz[cnti];
      cgammaz[cnti - 1] = - 1. / (Az0 + Az * calphaz[cnti - 1]);
   }

   return;
}

/**
 *    Initialization of the dipolar potential.
 *    kx  - array with the space mesh values in the x-direction in the K-space
 *    ky  - array with the space mesh values in the y-direction in the K-space
 *    kz  - array with the space mesh values in the z-direction in the K-space
 *    kx2 - array with the squared space mesh values in the x-direction in the
 *          K-space
 *    ky2 - array with the squared space mesh values in the y-direction in the
 *          K-space
 *    kz2 - array with the squared space mesh values in the z-direction in the
 *          K-space
 */
void initpotdd(double *kx, double *ky, double *kz, double *kx2, double *ky2, double *kz2) {
   long cnti, cntj, cntk;
   double dkx, dky, dkz, xk, tmp;

   dkx = 2. * pi / (Nx * dx);
   dky = 2. * pi / (Ny * dy);
   dkz = 2. * pi / (Nz * dz);

   for (cnti = 0; cnti < Nx2; cnti ++) kx[cnti] = cnti * dkx;
   for (cnti = Nx2; cnti < Nx; cnti ++) kx[cnti] = (cnti - Nx) * dkx;
   for (cntj = 0; cntj < Ny2; cntj ++) ky[cntj] = cntj * dky;
   for (cntj = Ny2; cntj < Ny; cntj ++) ky[cntj] = (cntj - Ny) * dky;
   for (cntk = 0; cntk < Nz2; cntk ++) kz[cntk] = cntk * dkz;
   for (cntk = Nz2; cntk < Nz; cntk ++) kz[cntk] = (cntk - Nz) * dkz;

   for (cnti = 0; cnti < Nx; cnti ++) kx2[cnti] = kx[cnti] * kx[cnti];
   for (cntj = 0; cntj < localNy; cntj ++) ky2[cntj] = ky[offsetNy + cntj] * ky[offsetNy + cntj];
   for (cntk = 0; cntk < Nz; cntk ++) kz2[cntk] = kz[cntk] * kz[cntk];

   #pragma omp parallel for private(cnti, cntj, cntk, xk, tmp)
   for (cntj = 0; cntj < localNy; cntj ++) {
      for (cnti = 0; cnti < Nx; cnti ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            xk = sqrt(kz2[cntk] + kx2[cnti] + ky2[cntj]);
            tmp = 1. + 3. * cos(xk * cutoff) / (xk * xk * cutoff * cutoff) - 3. * sin(xk * cutoff) / (xk * xk * xk * cutoff * cutoff * cutoff);
            potdd[cntj][cnti][cntk] = (4. * pi * (3. * kx2[cnti] / (kx2[cnti] + ky2[cntj] + kz2[cntk]) - 1.) / 3.) * tmp;
         }
      }
   }

   if (rank == 0) {
      potdd[0][0][0] = 0.;
   }

   return;
}



/**
 *    Calculation of the wave function norm and normalization.
 *    norm - wave function norm
 *    psi  - array with the wave function values
 *    tmpx - temporary array
 *    tmpy - temporary array
 *    tmpz - temporary array
 */
void calcnorm(double *norm, double complex ***psi, double **tmpx, double **tmpy, double **tmpz) {
   int threadid, quick;
   long cnti, cntj, cntk;
   double tmp, alpha;
   void *sendbuf;

   #pragma omp parallel private(threadid, cnti, cntj, cntk, tmp)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               tmp = cabs(psi[cnti][cntj][cntk]);
               tmpz[threadid][cntk] = tmp * tmp;
            }
            tmpy[threadid][cntj] = simpint(dz, tmpz[threadid], Nz);
         }
         (*tmpx)[cnti] = simpint(dy, tmpy[threadid], Ny);
      }
   }

   sendbuf = (rank == 0) ? MPI_IN_PLACE : *tmpx;
   MPI_Gather(sendbuf, localNx, MPI_DOUBLE, *tmpx, localNx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if (rank == 0) {
      norm_psi2 = simpint(dx, *tmpx, Nx);
      alpha = 1. / norm_psi2;
      *norm = sqrt(alpha);
   }

   if (QDEPL == 1) {
   #pragma omp parallel private(threadid, cnti, cntj, cntk, tmp)
   {
       threadid = omp_get_thread_num();

       #pragma omp for
       for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               tmp = cabs(psi[cnti][cntj][cntk]);
               tmpz[threadid][cntk] = tmp * tmp * tmp;
            }
            tmpy[threadid][cntj] = simpint(dz, tmpz[threadid], Nz);
         }
         (*tmpx)[cnti] = simpint(dy, tmpy[threadid], Ny);
      }
   }

   sendbuf = (rank == 0) ? MPI_IN_PLACE : *tmpx;
   MPI_Gather(sendbuf, localNx, MPI_DOUBLE, *tmpx, localNx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if (rank == 0) {
      norm_psi3 = simpint(dx, *tmpx, Nx);

      for (quick = 0; quick < 50; quick ++) {
         alpha = (1. - pow(alpha, 1.5) * h4 * norm_psi3) / norm_psi2;
      }

      *norm = sqrt(alpha);
   }
   }

   MPI_Bcast(norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   #pragma omp for private(cnti, cntj, cntk)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            psi[cnti][cntj][cntk] *= *norm;
         }
      }
   }

   return;
}

/**
 *    Calculation of the chemical potential and energy.
 *    mu        - chemical potential
 *    en        - energy
 *    psi       - array with the wave function values
 *    psi_t       - array with the transposed wave function values
 *    dpsi      - temporary array
 *    dpsi_t    - temporary array
 *    psidd2    - array with the squared wave function values
 *    psidd2fft - array with the squared wave function fft values
 *    tmpxi     - temporary array
 *    tmpyi     - temporary array
 *    tmpzi     - temporary array
 *    tmpxj     - temporary array
 *    tmpyj     - temporary array
 *    tmpzj     - temporary array
 */
void calcmudet(double *mu, double complex ***psi, double complex ***psi_t, double ***dpsi, double ***dpsi_t, double ***psidd2, double ***psidd20, fftw_complex *psidd2fft, double **tmpxi, double **tmpyi, double **tmpzi, double **tmpxj, double **tmpyj, double **tmpzj, double complex **tmpdpsi1, double complex **tmpdpsi2, double ***tmpxmui, double ***tmpymui, double ***tmpzmui) {
   int threadid;
   long cnti, cntj, cntk, cnte;
   double psi2, psi3, kpar, tmp;
   void *sendbuf;

   kpar = 0.5 * par;

   //fftw_execute(plan_transpose_x);

   #pragma omp parallel private(threadid, cnti, cntj, cntk)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cntj = 0; cntj < localNy; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            for (cnti = 0; cnti < Nx; cnti ++) {
               tmpdpsi1[threadid][cnti] = psi_t[cnti][cntj][cntk];
            }
            diffc(dx, tmpdpsi1[threadid], tmpdpsi2[threadid], Nx);
            for (cnti = 0; cnti < Nx; cnti ++) {
               dpsi_t[cnti][cntj][cntk] = tmpdpsi2[threadid][cnti] * conj(tmpdpsi2[threadid][cnti]);
            }
         }
      }
   }

   fftw_execute(plan_transpose_dpsi);

   #pragma omp parallel private(threadid, cnti, cntj, cntk)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            for (cntj = 0; cntj < Ny; cntj ++) {
               tmpdpsi1[threadid][cntj] = psi[cnti][cntj][cntk];
            }
            diffc(dy, tmpdpsi1[threadid], tmpdpsi2[threadid], Ny);
            for (cntj = 0; cntj < Ny; cntj ++) {
               dpsi[cnti][cntj][cntk] += tmpdpsi2[threadid][cntj] * conj(tmpdpsi2[threadid][cntj]);
            }
         }
      }
      #pragma omp barrier

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               tmpdpsi1[threadid][cntk] = psi[cnti][cntj][cntk];
            }
            diffc(dz, tmpdpsi1[threadid], tmpdpsi2[threadid], Nz);
            for (cntk = 0; cntk < Nz; cntk ++) {
               dpsi[cnti][cntj][cntk] += tmpdpsi2[threadid][cntk] * conj(tmpdpsi2[threadid][cntk]);
            }
         }
      }
   }

   calcpsidd2(psi, psidd2, psidd2fft);
   calcpsidd20(psi, psidd20, psidd2fft);

   #pragma omp parallel private(threadid, cnti, cntj, cntk, cnte, psi2, psi3, tmp)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               tmp = cabs(psi[cnti][cntj][cntk]);
               psi2 = tmp * tmp;
               psi3 = psi2 * tmp;
               tmpzmui[0][threadid][cntk] = kpar * dpsi[cnti][cntj][cntk];
               tmpzmui[1][threadid][cntk] = pot[cnti][cntj][cntk] * psi2;
               tmpzmui[2][threadid][cntk] = g * psi2 * psi2;
               tmpzmui[3][threadid][cntk] = gd * psidd20[cnti][cntj][cntk] * psi2;
               tmpzmui[4][threadid][cntk] = h2 * psi3 * psi2;
               tmpzmui[5][threadid][cntk] = gd * (psidd2[cnti][cntj][cntk] - psidd20[cnti][cntj][cntk]) * psi2;
            }
            for (cnte = 0; cnte < 6; cnte ++) tmpymui[cnte][threadid][cntj] = simpint(dz, tmpzmui[cnte][threadid], Nz);
         }
         for (cnte = 0; cnte < 6; cnte ++) tmpxmui[cnte][0][cnti] = simpint(dy, tmpymui[cnte][threadid], Ny);
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);

   for (cnte = 0; cnte < 6; cnte ++) {
      sendbuf = (rank == 0) ? MPI_IN_PLACE : tmpxmui[cnte][0];
      MPI_Gather(sendbuf, localNx, MPI_DOUBLE, tmpxmui[cnte][0], localNx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   }

   if (rank == 0) for (cnte = 0; cnte < 6; cnte ++) mu[cnte] = simpint(dx, tmpxmui[cnte][0], Nx);

   return;
}

void calcmu(double *mu, double complex ***psi, double complex ***psi_t, double ***dpsi, double ***dpsi_t, double ***psidd2, fftw_complex *psidd2fft, double **tmpxi, double **tmpyi, double **tmpzi, double **tmpxj, double **tmpyj, double **tmpzj, double complex **tmpdpsi1, double complex **tmpdpsi2) {
   int threadid;
   long cnti, cntj, cntk;
   double dpsi2, psi2, psi3, psi2lin, psidd2lin, kpar, tmp;
   void *sendbuf;

   kpar = 0.5 * par;
   calcpsidd2(psi, psidd2, psidd2fft);

   //fftw_execute(plan_transpose_x);

   #pragma omp parallel private(threadid, cnti, cntj, cntk)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cntj = 0; cntj < localNy; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            for (cnti = 0; cnti < Nx; cnti ++) {
               tmpdpsi1[threadid][cnti] = psi_t[cnti][cntj][cntk];
            }
            diffc(dx, tmpdpsi1[threadid], tmpdpsi2[threadid], Nx);
            for (cnti = 0; cnti < Nx; cnti ++) {
               dpsi_t[cnti][cntj][cntk] = tmpdpsi2[threadid][cnti] * conj(tmpdpsi2[threadid][cnti]);
            }
         }
      }
   }

   fftw_execute(plan_transpose_dpsi);

   #pragma omp parallel private(threadid, cnti, cntj, cntk, dpsi2, psi2, psi3, psi2lin, psidd2lin, tmp)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            for (cntj = 0; cntj < Ny; cntj ++) {
               tmpdpsi1[threadid][cntj] = psi[cnti][cntj][cntk];
            }
            diffc(dy, tmpdpsi1[threadid], tmpdpsi2[threadid], Ny);
            for (cntj = 0; cntj < Ny; cntj ++) {
               dpsi[cnti][cntj][cntk] += tmpdpsi2[threadid][cntj] * conj(tmpdpsi2[threadid][cntj]);
            }
         }
      }
      #pragma omp barrier

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               tmpdpsi1[threadid][cntk] = psi[cnti][cntj][cntk];
            }
            diffc(dz, tmpdpsi1[threadid], tmpdpsi2[threadid], Nz);
            for (cntk = 0; cntk < Nz; cntk ++) {
               dpsi[cnti][cntj][cntk] += tmpdpsi2[threadid][cntk] * conj(tmpdpsi2[threadid][cntk]);
            }
         }
      }
      #pragma omp barrier

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               dpsi2 = kpar * dpsi[cnti][cntj][cntk];
               tmp = cabs(psi[cnti][cntj][cntk]);
               psi2 = tmp * tmp;
               psi3 = psi2 * tmp;
               psi2lin = g * psi2;
               psidd2lin = gd * psidd2[cnti][cntj][cntk];
               tmpzi[threadid][cntk] = dpsi2 + (pot[cnti][cntj][cntk] + psi2lin + psidd2lin + h2 * psi3) * psi2;
            }
            tmpyi[threadid][cntj] = simpint(dz, tmpzi[threadid], Nz);
         }
         (*tmpxi)[cnti] = simpint(dy, tmpyi[threadid], Ny);
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);

   sendbuf = (rank == 0) ? MPI_IN_PLACE : *tmpxi;
   MPI_Gather(sendbuf, localNx, MPI_DOUBLE, *tmpxi, localNx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if (rank == 0) *mu = simpint(dx, *tmpxi, Nx);

   return;
}


/**
 *    Calculation of squared wave function values for dipole-dipole
 *    interaction.
 *    psi       - array with the wave function values
 *    psidd2    - array with the squared wave function values
 *    psidd2fft - array with the squared wave function fft values
 */
void calcpsidd2(double complex ***psi, double ***psidd2, fftw_complex *psidd2fft) {
   long cnti, cntj, cntk;
   long last = 0;
   double *psidd2tmp = (double *) psidd2fft;
   double tmp;

   #pragma omp parallel for private(cnti, cntj, cntk, tmp)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            tmp = cabs(psi[cnti][cntj][cntk]);
            psidd2[cnti][cntj][cntk] = tmp * tmp * (1. + h4 * tmp);
         }
      }
   }

   fftw_execute(plan_forward);

   #pragma omp parallel for private(cnti, cntj, cntk)
   for (cntj = 0; cntj < localNy; cntj ++) {
      for (cnti = 0; cnti < Nx; cnti ++) {
         for (cntk = 0; cntk < Nz2 + 1; cntk ++) {
            psidd2fft[cntj * Nx * (Nz2 + 1) + cnti * (Nz2 + 1) + cntk][0] *= potdd[cntj][cnti][cntk];
            psidd2fft[cntj * Nx * (Nz2 + 1) + cnti * (Nz2 + 1) + cntk][1] *= potdd[cntj][cnti][cntk];
         }
      }
   }

   fftw_execute(plan_backward);

   #pragma omp parallel for private(cnti, cntj, cntk)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            psidd2[cnti][cntj][cntk] /= (Nx * Ny * Nz);
         }
      }
   }

   if (nprocs > 1) {
      if (rank == 0) {
         MPI_Send(psidd2[0][0], Ny * Nz, MPI_DOUBLE, nprocs - 1, 0, MPI_COMM_WORLD);
      } else if (rank == nprocs - 1) {
         MPI_Recv(psidd2tmp, Ny * Nz, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         last = 1;
      }
   } else {
      psidd2tmp = psidd2[0][0];
      last = 1;
   }

   if (rank == nprocs - 1) {
      #pragma omp parallel for private(cntj, cntk)
      for (cntj = 0; cntj < Ny - 1; cntj ++) {
         for (cntk = 0; cntk < Nz - 1; cntk ++) {
            psidd2[localNx - 1][cntj][cntk] = psidd2tmp[cntj * Nz + cntk];
         }
      }
   }

   #pragma omp parallel for private(cnti, cntk)
   for (cnti = 0; cnti < localNx - last; cnti ++) {
      for (cntk = 0; cntk < Nz - 1; cntk ++) {
         psidd2[cnti][Ny - 1][cntk] = psidd2[cnti][0][cntk];
      }
   }

   #pragma omp parallel for private(cnti, cntj)
   for (cnti = 0; cnti < localNx - last; cnti ++) {
      for (cntj = 0; cntj < Ny - 1; cntj ++) {
         psidd2[cnti][cntj][Nz - 1] = psidd2[cnti][cntj][0];
      }
   }

   return;
}

void calcpsidd20(double complex ***psi, double ***psidd20, fftw_complex *psidd2fft) {            
   long cnti, cntj, cntk;
   long last = 0;
   double *psidd2tmp = (double *) psidd2fft;
   double tmp;

   #pragma omp parallel for private(cnti, cntj, cntk, tmp)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            tmp = cabs(psi[cnti][cntj][cntk]);
            psidd20[cnti][cntj][cntk] = tmp * tmp;
         }
      }
   }

   fftw_execute(plan_forward0);

   #pragma omp parallel for private(cnti, cntj, cntk)
   for (cntj = 0; cntj < localNy; cntj ++) {
      for (cnti = 0; cnti < Nx; cnti ++) {
         for (cntk = 0; cntk < Nz2 + 1; cntk ++) {
            psidd2fft[cntj * Nx * (Nz2 + 1) + cnti * (Nz2 + 1) + cntk][0] *= potdd[cntj][cnti][cntk];
            psidd2fft[cntj * Nx * (Nz2 + 1) + cnti * (Nz2 + 1) + cntk][1] *= potdd[cntj][cnti][cntk];
         }
      }
   }

   fftw_execute(plan_backward0);

   #pragma omp parallel for private(cnti, cntj, cntk)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            psidd20[cnti][cntj][cntk] /= (Nx * Ny * Nz);
         }
      }
   }

   if (nprocs > 1) {
      if (rank == 0) {
         MPI_Send(psidd20[0][0], Ny * Nz, MPI_DOUBLE, nprocs - 1, 0, MPI_COMM_WORLD);
      } else if (rank == nprocs - 1) {
         MPI_Recv(psidd2tmp, Ny * Nz, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         last = 1;
      }
   } else {
      psidd2tmp = psidd20[0][0];
      last = 1;
   }

   if (rank == nprocs - 1) {
      #pragma omp parallel for private(cntj, cntk)
      for (cntj = 0; cntj < Ny - 1; cntj ++) {
         for (cntk = 0; cntk < Nz - 1; cntk ++) {
            psidd20[localNx - 1][cntj][cntk] = psidd2tmp[cntj * Nz + cntk];
         }
      }
   }

   #pragma omp parallel for private(cnti, cntk)
   for (cnti = 0; cnti < localNx - last; cnti ++) {
      for (cntk = 0; cntk < Nz - 1; cntk ++) {
         psidd20[cnti][Ny - 1][cntk] = psidd20[cnti][0][cntk];
      }
   }

   #pragma omp parallel for private(cnti, cntj)
   for (cnti = 0; cnti < localNx - last; cnti ++) {
      for (cntj = 0; cntj < Ny - 1; cntj ++) {
         psidd20[cnti][cntj][Nz - 1] = psidd20[cnti][cntj][0];
      }
   }

   return;
}

/**
 *    Calculation of the root mean square radius.
 *    rms   - root mean square radius
 *    psi   - array with the wave function values
 *    psi_t - array with the transposed wave function values
 *    tmpx  - temporary array
 *    tmpy  - temporary array
 *    tmpz  - temporary array
 */
void calcrms(double *rms, double complex ***psi, double complex ***psi_t, double **tmpx, double **tmpy, double **tmpz) {
   int threadid;
   long cnti, cntj, cntk;
   double psi2;
   void *sendbuf;

   fftw_execute(plan_transpose_x);

   #pragma omp parallel private(threadid, cnti, cntj, cntk, psi2)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cntj = 0; cntj < localNy; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            for (cnti = 0; cnti < Nx; cnti ++) {
               psi2 = cabs(psi_t[cnti][cntj][cntk]);
               psi2 *= psi2;
               tmpx[threadid][cnti] = x2[cnti] * psi2;
            }
            tmpz[threadid][cntk] = simpint(dx, tmpx[threadid], Nx);
         }
         (*tmpy)[cntj] = simpint(dz, tmpz[threadid], Nz);
      }
   }

   sendbuf = (rank == 0) ? MPI_IN_PLACE : *tmpy;
   MPI_Gather(sendbuf, localNy, MPI_DOUBLE, *tmpy, localNy, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if (rank == 0) {
      rms[1] = sqrt(simpint(dy, *tmpy, Ny));
   }

   #pragma omp parallel private(threadid, cnti, cntj, cntk, psi2)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            for (cntj = 0; cntj < Ny; cntj ++) {
               psi2 = cabs(psi[cnti][cntj][cntk]);
               psi2 *= psi2;
               tmpy[threadid][cntj] = y2[cntj] * psi2;
            }
            tmpz[threadid][cntk] = simpint(dy, tmpy[threadid], Ny);
         }
         (*tmpx)[cnti] = simpint(dz, tmpz[threadid], Nz);
      }
   }

   sendbuf = (rank == 0) ? MPI_IN_PLACE : *tmpx;
   MPI_Gather(sendbuf, localNx, MPI_DOUBLE, *tmpx, localNx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if (rank == 0) {
      rms[2] = sqrt(simpint(dx, *tmpx, Nx));
   }

   #pragma omp parallel private(threadid, cnti, cntj, cntk, psi2)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               psi2 = cabs(psi[cnti][cntj][cntk]);
               psi2 *= psi2;
               tmpz[threadid][cntk] = z2[cntk] * psi2;
            }
            tmpy[threadid][cntj] = simpint(dz, tmpz[threadid], Nz);
         }
         (*tmpx)[cnti] = simpint(dy, tmpy[threadid], Ny);
      }
   }

   sendbuf = (rank == 0) ? MPI_IN_PLACE : *tmpx;
   MPI_Gather(sendbuf, localNx, MPI_DOUBLE, *tmpx, localNx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if (rank == 0) {
      rms[3] = sqrt(simpint(dx, *tmpx, Nx));
      rms[0] = sqrt(rms[1] * rms[1] + rms[2] * rms[2] + rms[3] * rms[3]);
   }

   return;
}
/**
 *    Time propagation with respect to H1 (part of the Hamiltonian without
 *    spatial derivatives).
 *    psi       - array with the wave function values
 *    psidd2    - array with the squared wave function values
 *    psidd2fft - array with the squared wave function fft values
 */
void calcnu(double complex ***psi, double ***psidd2, fftw_complex *psidd2fft) {
   long cnti, cntj, cntk;
   double psi2, psi3, psi2lin, psidd2lin, tmp;

   calcpsidd2(psi, psidd2, psidd2fft);

   #pragma omp parallel for private(cnti, cntj, cntk, psi2, psi3, psi2lin, psidd2lin, tmp)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            tmp = cabs(psi[cnti][cntj][cntk]);
            psi2 = tmp * tmp;
            psi3 = psi2 * tmp;
            psi2lin = g * psi2;
            psidd2lin = gd * psidd2[cnti][cntj][cntk];
            tmp = dt * (pot[cnti][cntj][cntk] + psi2lin + psidd2lin + h2 * psi3);
            psi[cnti][cntj][cntk] *= exp(- tmp);
         }
      }
   }

   return;
}

/**
 *    Time propagation with respect to H2 (x-part of the Laplacian).
 *    psi_t - array with the wave function values (transposed)
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calclux(double complex ***psi_t, double complex **cbeta) {
   int threadid;
   long cnti, cntj, cntk;
   double complex c;

   fftw_execute(plan_transpose_x);

   #pragma omp parallel private(threadid, cnti, cntj, cntk, c)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cntj = 0; cntj < localNy; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            cbeta[threadid][Nx - 2] = psi_t[Nx - 1][cntj][cntk];
            for (cnti = Nx - 2; cnti > 0; cnti --) {
               c = - Ax * psi_t[cnti + 1][cntj][cntk] + Ax0r * psi_t[cnti][cntj][cntk] - Ax * psi_t[cnti - 1][cntj][cntk];
               cbeta[threadid][cnti - 1] =  cgammax[cnti] * (Ax * cbeta[threadid][cnti] - c);
            }
            psi_t[0][cntj][cntk] = 0.;
            for (cnti = 0; cnti < Nx - 2; cnti ++) {
               psi_t[cnti + 1][cntj][cntk] = calphax[cnti] * psi_t[cnti][cntj][cntk] + cbeta[threadid][cnti];
            }
            psi_t[Nx - 1][cntj][cntk] = 0.;
         }
      }
   }

   fftw_execute(plan_transpose_y);

   return;
}

/**
 *    Time propagation with respect to H3 (y-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calcluy(double complex ***psi, double complex **cbeta) {
   int threadid;
   long cnti, cntj, cntk;
   double complex c;

   #pragma omp parallel private(threadid, cnti, cntj, cntk, c)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            cbeta[threadid][Ny - 2] = psi[cnti][Ny - 1][cntk];
            for (cntj = Ny - 2; cntj > 0; cntj --) {
               c = - Ay * psi[cnti][cntj + 1][cntk] + Ay0r * psi[cnti][cntj][cntk] - Ay * psi[cnti][cntj - 1][cntk];
               cbeta[threadid][cntj - 1] =  cgammay[cntj] * (Ay * cbeta[threadid][cntj] - c);
            }
            psi[cnti][0][cntk] = 0.;
            for (cntj = 0; cntj < Ny - 2; cntj ++) {
               psi[cnti][cntj + 1][cntk] = calphay[cntj] * psi[cnti][cntj][cntk] + cbeta[threadid][cntj];
            }
            psi[cnti][Ny - 1][cntk] = 0.;
         }
      }
   }

   return;
}

/**
 *    Time propagation with respect to H4 (z-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calcluz(double complex ***psi, double complex **cbeta) {
   int threadid;
   long cnti, cntj, cntk;
   double complex c;

   #pragma omp parallel private(threadid, cnti, cntj, cntk, c)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            cbeta[threadid][Nz - 2] = psi[cnti][cntj][Nz - 1];  // mozda ima neke razlike za real
            for (cntk = Nz - 2; cntk > 0; cntk --) {
               c = - Az * psi[cnti][cntj][cntk + 1] + Az0r * psi[cnti][cntj][cntk] - Az * psi[cnti][cntj][cntk - 1];
               cbeta[threadid][cntk - 1] =  cgammaz[cntk] * (Az * cbeta[threadid][cntk] - c);
            }
            psi[cnti][cntj][0] = 0.;
            for (cntk = 0; cntk < Nz - 2; cntk ++) {
               psi[cnti][cntj][cntk + 1] = calphaz[cntk] * psi[cnti][cntj][cntk] + cbeta[threadid][cntk];
            }
            psi[cnti][cntj][Nz - 1] = 0.;
         }
      }
   }

   return;
}
void outdenx(double complex ***psi, double **outx, double *tmpy, double *tmpz, MPI_File file) {
   long cnti, cntj, cntk;
   MPI_Offset fileoffset;

   fileoffset = rank * 2 * sizeof(double) * (localNx / outstpx);

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            tmpz[cntk] = cabs(psi[cnti][cntj][cntk]);
            tmpz[cntk] *= tmpz[cntk];
         }
         tmpy[cntj] = simpint(dz, tmpz, Nz);
      }
      outx[cnti / outstpx][0] = x[offsetNx + cnti];
      outx[cnti / outstpx][1] = simpint(dy, tmpy, Ny);
   }

   MPI_File_write_at_all(file, fileoffset, *outx, (localNx / outstpx) * 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void outdeny(double complex ***psi_t, double **outy, double *tmpx, double *tmpz, MPI_File file) {
   long cnti, cntj, cntk;
   MPI_Offset fileoffset;

   fileoffset = rank * 2 * sizeof(double) * (localNy / outstpy);

   fftw_execute(plan_transpose_x);

   for (cntj = 0; cntj < localNy; cntj += outstpy) {
      for (cnti = 0; cnti < Nx; cnti ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            tmpz[cntk] = cabs(psi_t[cnti][cntj][cntk]);
            tmpz[cntk] *= tmpz[cntk];
         }
         tmpx[cnti] = simpint(dz, tmpz, Nz);
      }
      outy[cntj / outstpy][0] = y[offsetNy + cntj];
      outy[cntj / outstpy][1] = simpint(dx, tmpx, Nx);
   }

   MPI_File_write_at_all(file, fileoffset, *outy, (localNy / outstpy) * 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void outdenz(double complex ***psi, double **outz, double *tmpx, double *tmpy, MPI_File file) {
   long cnti, cntj, cntk;
   MPI_Offset fileoffset;
   void *sendbuf;

   sendbuf = (rank == 0) ? MPI_IN_PLACE : tmpx;

   for (cntk = 0; cntk < Nz; cntk += outstpz) {
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            tmpy[cntj] = cabs(psi[cnti][cntj][cntk]);
            tmpy[cntj] *= tmpy[cntj];
         }
         tmpx[cnti] = simpint(dy, tmpy, Ny);
      }

      MPI_Gather(sendbuf, localNx, MPI_DOUBLE, tmpx, localNx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (rank == 0) {
         outz[cntk / outstpz][0] = z[cntk];
         outz[cntk / outstpz][1] = simpint(dx, tmpx, Nx);
      }
   }

   if (rank == 0) {
      fileoffset = 0;
      MPI_File_write_at(file, fileoffset, *outz, (Nz / outstpz) * 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
   }
}

void outdenxy(double complex ***psi, double ***outxy, double *tmpz, MPI_File file) {
   long cnti, cntj, cntk;
   MPI_Offset fileoffset;

   fileoffset = rank * 3 * sizeof(double) * (localNx / outstpx) * (Ny / outstpy);

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      for (cntj = 0; cntj < Ny; cntj += outstpy) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            tmpz[cntk] = cabs(psi[cnti][cntj][cntk]);
            tmpz[cntk] *= tmpz[cntk];
         }
         outxy[cnti / outstpx][cntj / outstpy][0] = x[offsetNx + cnti];
         outxy[cnti / outstpx][cntj / outstpy][1] = y[cntj];
         outxy[cnti / outstpx][cntj / outstpy][2] = simpint(dz, tmpz, Nz);
      }
   }

   MPI_File_write_at_all(file, fileoffset, **outxy, (localNx / outstpx) * (Ny / outstpy) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);

   return;
}

void outdenxz(double complex ***psi, double ***outxz, double *tmpy, MPI_File file) {
   long cnti, cntj, cntk;
   MPI_Offset fileoffset;

   fileoffset = rank * 3 * sizeof(double) * (localNx / outstpx) * (Nz / outstpz);

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      for (cntk = 0; cntk < Nz; cntk += outstpz) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            tmpy[cntj] = cabs(psi[cnti][cntj][cntk]);
            tmpy[cntj] *= tmpy[cntj];
         }
         outxz[cnti / outstpx][cntk / outstpz][0] = x[offsetNx + cnti];
         outxz[cnti / outstpx][cntk / outstpz][1] = z[cntk];
         outxz[cnti / outstpx][cntk / outstpz][2] = simpint(dy, tmpy, Ny);
      }
   }

   MPI_File_write_at_all(file, fileoffset, **outxz, (localNx / outstpx) * (Nz / outstpz) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);

   return;
}

void outdenyz(double complex ***psi_t, double ***outyz, double *tmpx, MPI_File file) {
   long cnti, cntj, cntk;
   MPI_Offset fileoffset;

   fileoffset = rank * 3 * sizeof(double) * (localNy / outstpy) * (Nz / outstpz);

   fftw_execute(plan_transpose_x);

   for (cntj = 0; cntj < localNy; cntj += outstpy) {
      for (cntk = 0; cntk < Nz; cntk += outstpz) {
         for (cnti = 0; cnti < Nx; cnti ++) {
            tmpx[cnti] = cabs(psi_t[cnti][cntj][cntk]);
            tmpx[cnti] *= tmpx[cnti];
         }
         outyz[cntj / outstpy][cntk / outstpz][0] = y[offsetNy + cntj];
         outyz[cntj / outstpy][cntk / outstpz][1] = z[cntk];
         outyz[cntj / outstpy][cntk / outstpz][2] = simpint(dx, tmpx, Nx);
      }
   }

   MPI_File_write_at_all(file, fileoffset, **outyz, (localNy / outstpy) * (Nz / outstpz) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void outpsi2xy(double complex ***psi, double ***outxy, MPI_File file) {
   long cnti, cntj;
   MPI_Offset fileoffset;

   fileoffset = rank * 3 * sizeof(double) * (localNx / outstpx) * (Ny / outstpy);

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      for (cntj = 0; cntj < Ny; cntj += outstpy) {
         outxy[cnti / outstpx][cntj / outstpy][0] = x[offsetNx + cnti];
         outxy[cnti / outstpx][cntj / outstpy][1] = y[cntj];
         outxy[cnti / outstpx][cntj / outstpy][2] = cabs(psi[cnti][cntj][Nz2]) * cabs(psi[cnti][cntj][Nz2]);
      }
   }

   MPI_File_write_at_all(file, fileoffset, **outxy, (localNx / outstpx) * (Ny / outstpy) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);

   return;
}

void outpsi2xz(double complex ***psi, double ***outxz, MPI_File file) {
   long cnti, cntk;
   MPI_Offset fileoffset;

   fileoffset = rank * 3 * sizeof(double) * (localNx / outstpx) * (Nz / outstpz);

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      for (cntk = 0; cntk < Nz; cntk += outstpz) {
         outxz[cnti / outstpx][cntk / outstpz][0] = x[offsetNx + cnti];
         outxz[cnti / outstpx][cntk / outstpz][1] = z[cntk];
         outxz[cnti / outstpx][cntk / outstpz][2] = cabs(psi[cnti][Ny2][cntk]) * cabs(psi[cnti][Ny2][cntk]);;
      }
   }

   MPI_File_write_at_all(file, fileoffset, **outxz, (localNx / outstpx) * (Nz / outstpz) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);

   return;
}

void outpsi2yz(double complex ***psi, double ***outyz, MPI_File file) {
   long cntj, cntk;
   int rankNx2, offsetNx2;
   MPI_Offset fileoffset;

   rankNx2 = Nx2 / localNx;
   offsetNx2 = Nx2 % localNx;

   fileoffset = 0;

   if (rank == rankNx2) {
      for (cntj = 0; cntj < Ny; cntj += outstpy) {
         for (cntk = 0; cntk < Nz; cntk += outstpz) {
            outyz[cntj / outstpy][cntk / outstpz][0] = y[cntj];
            outyz[cntj / outstpy][cntk / outstpz][1] = z[cntk];
            outyz[cntj / outstpy][cntk / outstpz][2] = cabs(psi[offsetNx2][cntj][cntk]) * cabs(psi[offsetNx2][cntj][cntk]);;
         }
      }

      MPI_File_write_at(file, fileoffset, **outyz, (Ny / outstpy) * (Nz / outstpz) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);
   }
}

void outdenxyz(double complex ***psi, double ***outxyz, MPI_File file) {
   long cnti, cntj, cntk;
   MPI_Offset fileoffset;

   // MPI IO returns error if the array is too large. As a workaround, we write just Ny * Nz at a time.

   fileoffset = rank * sizeof(double) * localNx * Ny * Nz;

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      for (cntj = 0; cntj < Ny; cntj += outstpy) {
         for (cntk = 0; cntk < Nz; cntk += outstpz) {
            outxyz[0][cntj][cntk] = cabs(psi[cnti][cntj][cntk]) * cabs(psi[cnti][cntj][cntk]);
         }
      }

      MPI_File_write_at_all(file, fileoffset, **outxyz, (Ny / outstpy) * (Nz / outstpz), MPI_DOUBLE, MPI_STATUS_IGNORE);
      fileoffset += (Ny / outstpy) * (Nz / outstpz) * sizeof(double);
   }

   return;
}

void outpsi(double complex ***psi, MPI_File file) {
   long cnti;
   MPI_Offset fileoffset;

   // MPI IO returns error if the array is too large. As a workaround, we write just Ny * Nz at a time.

   fileoffset = rank * sizeof(double complex) * localNx * Ny * Nz;

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      MPI_File_write_at_all(file, fileoffset, *(psi[cnti]), 2 * (Ny / outstpy) * (Nz / outstpz), MPI_DOUBLE, MPI_STATUS_IGNORE);
      fileoffset += (Ny / outstpy) * (Nz / outstpz) * sizeof(double complex);
   }

   return;
}

void outargxy(double complex ***psi, double ***outxy, MPI_File file) {
   long cnti, cntj;
   MPI_Offset fileoffset;

   fileoffset = rank * 3 * sizeof(double) * (localNx / outstpx) * (Ny / outstpy);

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      for (cntj = 0; cntj < Ny; cntj += outstpy) {
         outxy[cnti / outstpx][cntj / outstpy][0] = x[offsetNx + cnti];
         outxy[cnti / outstpx][cntj / outstpy][1] = y[cntj];
         outxy[cnti / outstpx][cntj / outstpy][2] = carg(psi[cnti][cntj][Nz2]);
      }
   }

   MPI_File_write_at_all(file, fileoffset, **outxy, (localNx / outstpx) * (Ny / outstpy) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);

   return;
}

void outargxz(double complex ***psi, double ***outxz, MPI_File file) {
   long cnti, cntk;
   MPI_Offset fileoffset;

   fileoffset = rank * 3 * sizeof(double) * (localNx / outstpx) * (Nz / outstpz);

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      for (cntk = 0; cntk < Nz; cntk += outstpz) {
         outxz[cnti / outstpx][cntk / outstpz][0] = x[offsetNx + cnti];
         outxz[cnti / outstpx][cntk / outstpz][1] = z[cntk];
         outxz[cnti / outstpx][cntk / outstpz][2] = carg(psi[cnti][Ny2][cntk]);
      }
   }

   MPI_File_write_at_all(file, fileoffset, **outxz, (localNx / outstpx) * (Nz / outstpz) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);

   return;
}

void outargyz(double complex ***psi, double ***outyz, MPI_File file) {
   long cntj, cntk;
   int rankNx2, offsetNx2;
   MPI_Offset fileoffset;

   rankNx2 = Nx2 / localNx;
   offsetNx2 = Nx2 % localNx;

   fileoffset = 0;

   if (rank == rankNx2) {
      for (cntj = 0; cntj < Ny; cntj += outstpy) {
         for (cntk = 0; cntk < Nz; cntk += outstpz) {
            outyz[cntj / outstpy][cntk / outstpz][0] = y[cntj];
            outyz[cntj / outstpy][cntk / outstpz][1] = z[cntk];
            outyz[cntj / outstpy][cntk / outstpz][2] = carg(psi[offsetNx2][cntj][cntk]);
         }
      }

      MPI_File_write_at(file, fileoffset, **outyz, (Ny / outstpy) * (Nz / outstpz) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);
   }
}
