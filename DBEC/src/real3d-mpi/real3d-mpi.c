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

#include "real3d-mpi.h"

int main(int argc, char **argv) {
   FILE *out;
   FILE *filerms;
   MPI_File mpifile;
   char filename[MAX_FILENAME_SIZE];
   int nthreads, rankNx2;
   long offsetNx2;
   long cnti;
   double norm, mu, en;
   double *rms;
   double complex **cbeta;
   double complex ***psi, ***psi_t;
   double ***dpsi, ***dpsi_t;
   double ***psidd2;
   double **tmpxi, **tmpyi, **tmpzi, **tmpxj, **tmpyj, **tmpzj;
   double **outx, **outy, **outz;
   double ***outxy, ***outxz, ***outyz;
   double ***outxyz;
   fftw_complex *psidd2fft;
   double psi2;

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

   dpsi = alloc_double_tensor(localNx, Ny, Nz);
   dpsi_t = alloc_double_tensor(Nx, localNy, Nz);

   calphax = alloc_complex_vector(Nx - 1);
   calphay = alloc_complex_vector(Ny - 1);
   calphaz = alloc_complex_vector(Nz - 1);
   cbeta = alloc_complex_matrix(nthreads, MAX(Nx, Ny, Nz) - 1);
   cgammax = alloc_complex_vector(Nx - 1);
   cgammay = alloc_complex_vector(Ny - 1);
   cgammaz = alloc_complex_vector(Nz - 1);

   tmpxi = alloc_double_matrix(nthreads, Nx);
   tmpyi = alloc_double_matrix(nthreads, Ny);
   tmpzi = alloc_double_matrix(nthreads, Nz);
   tmpxj = alloc_double_matrix(nthreads, Nx);
   tmpyj = alloc_double_matrix(nthreads, Ny);
   tmpzj = alloc_double_matrix(nthreads, Nz);

   psidd2fft = alloc_fftw_complex_vector(localNx * Ny * (Nz2 + 1));

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

   g = par * g0;
   gd = par * gd0;


   if (rank == 0) {
      if (output != NULL) {
         sprintf(filename, "%s.txt", output);
         out = fopen(filename, "w");
      } else out = stdout;
   } else out = fopen("/dev/null", "w");

   if (rank == 0) {
      if (rmsout != NULL) {
         sprintf(filename, "%s.txt", rmsout);
         filerms = fopen(filename, "w");
      } else filerms = NULL;
   } else filerms = fopen("/dev/null", "w");


   fprintf(out, "  Real time propagation 3D, OPTION = %d, OMP_NUM_THREADS = %d, MPI_NUM_PROCS = %d\n\n", opt, nthreads, nprocs);
   fprintf(out, "  Number of Atoms N = %li, Unit of length AHO = %.8f m\n", Na, aho);
   fprintf(out, "  Scattering length a = %.2f*a0, Dipolar ADD = %.2f*a0\n", as, add);
   fprintf(out, "  Nonlinearity G_3D = %.4f, Strength of DDI GD_3D = %.5f\n", g0, gd0);
   fprintf(out, "  Dipolar Cut off:   R = %.3f\n",  cutoff);
   fprintf(out, "  Parameters of trap: GAMMA = %.2f, NU = %.2f, LAMBDA = %.2f\n\n", vgamma, vnu, vlambda);
   fprintf(out, "  Space Step: NX = %li, NY = %li, NZ = %li\n", Nx, Ny, Nz);
   fprintf(out, "              DX = %.6f, DY = %.6f, DZ = %.6f\n", dx, dy, dz);
   fprintf(out, "  Time Step : NPAS = %li, ITER = %li\n", Npas, Npas/iter);
   fprintf(out, "              DT = %.6f\n\n",  dt);
   fprintf(out, "            ---------------------------------------------------------------------\n");
   fprintf(out, "              Iter     Norm       Chem       Ener/N       <r>    |Psi(0,0,0)|^2\n");
   fprintf(out, "            ---------------------------------------------------------------------\n");
   fflush(out);


   fftw_init_threads();
   fftw_mpi_init();
   fftw_plan_with_nthreads(nthreads);

   long n[] = { Nx, Ny, Nz };
   plan_forward = fftw_mpi_plan_many_dft_r2c(3, n, 1, localNx, localNy, psidd2[0][0], psidd2fft, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_OUT);
   plan_backward = fftw_mpi_plan_many_dft_c2r(3, n, 1, localNy, localNx, psidd2fft, psidd2[0][0], MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_IN);

   // Alternately, FFT plans can be created in this way if the block sizes are the same.
   //plan_forward = fftw_mpi_plan_dft_r2c_3d(Nx, Ny, Nz, **psidd2, psidd2fft, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_OUT);
   //plan_backward = fftw_mpi_plan_dft_c2r_3d(Nx, Ny, Nz, psidd2fft, **psidd2, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_IN);

   plan_transpose_x = fftw_mpi_plan_many_transpose(Nx, Ny * Nz, 2, localNx, localNy * Nz, (double *) **psi, (double *) **psi_t, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_OUT);
   plan_transpose_y = fftw_mpi_plan_many_transpose(Ny * Nz, Nx, 2, localNy * Nz, localNx, (double *) **psi_t, (double *) **psi, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_IN);
   plan_transpose_dpsi = fftw_mpi_plan_many_transpose(Ny * Nz, Nx, 1, localNy * Nz, localNx, (double *) **dpsi_t, **dpsi, MPI_COMM_WORLD, FFT_FLAG | FFTW_MPI_TRANSPOSED_IN);

   initpsi(psi, *tmpzi);
   initpot();
   gencoef();
   initpotdd(*tmpxi, *tmpyi, *tmpzi, *tmpxj, *tmpyj, *tmpzj);


   calcnorm(&norm, psi, tmpxi, tmpyi, tmpzi);
   calcmuen(&mu, &en, psi, psi_t, dpsi, dpsi_t, psidd2, psidd2fft, tmpxi, tmpyi, tmpzi, tmpxj, tmpyj, tmpzj);
   calcrms(rms, psi, psi_t, tmpxi, tmpyi, tmpzi);

   psi2 = cabs(psi[offsetNx2][Ny2][Nz2]) * cabs(psi[offsetNx2][Ny2][Nz2]);
   MPI_Bcast(&psi2, 1, MPI_DOUBLE, rankNx2, MPI_COMM_WORLD);

   fprintf(out, "Initial: %19.4f %11.5f %11.5f %10.5f %10.5f\n", norm, mu / par, en / par, *rms, psi2);
   fflush(out);

   if (initout != NULL) {
      /*
      sprintf(filename, "%s.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outdenxyz(psi, outxyz, mpifile);
      MPI_File_close(&mpifile);
      */

      sprintf(filename, "%s1d_x.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outdenx(psi, outx, *tmpyi, *tmpzi, mpifile);
      MPI_File_close(&mpifile);

      sprintf(filename, "%s1d_y.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outdeny(psi_t, outy, *tmpxi, *tmpzi, mpifile);
      MPI_File_close(&mpifile);

      sprintf(filename, "%s1d_z.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outdenz(psi, outz, *tmpxi, *tmpyi, mpifile);
      MPI_File_close(&mpifile);

      sprintf(filename, "%s2d_xy.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outdenxy(psi, outxy, *tmpzi, mpifile);
      MPI_File_close(&mpifile);

      sprintf(filename, "%s2d_xz.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outdenxz(psi, outxz, *tmpyi, mpifile);
      MPI_File_close(&mpifile);

      sprintf(filename, "%s2d_yz.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outdenyz(psi_t, outyz, *tmpxi, mpifile);
      MPI_File_close(&mpifile);
      /*
      sprintf(filename, "%s3d_xy0.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outpsi2xy(psi, outxy, mpifile);
      MPI_File_close(&mpifile);

      sprintf(filename, "%s3d_x0z.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outpsi2xz(psi, outxz, mpifile);
      MPI_File_close(&mpifile);

      sprintf(filename, "%s3d_0yz.bin", initout);
      MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
      outpsi2yz(psi, outyz, mpifile);
      MPI_File_close(&mpifile);
      */
   }

   if (rmsout != NULL) {
      fprintf(filerms, "  Real time propagation 3D,   OPTION = %d\n\n", opt);
      fprintf(filerms, "           ------------------------------------------------------------\n");
      fprintf(filerms, "RMS size:    Iter      <r>          <x>          <y>          <z>\n");
      fprintf(filerms, "           ------------------------------------------------------------\n");
      fprintf(filerms, "Initial: %21.5f %12.5f %12.5f %12.5f\n", rms[0], rms[1], rms[2], rms[3]);
      fflush(filerms);
   }

   // ------------------------------------------------------------------- NPAS

   gettimeofday(&iter_start, NULL);

   int step = 1;
   for (cnti = 0; cnti < Npas; cnti ++) {
      calcnu(psi, psidd2, psidd2fft);
      calcluy(psi, cbeta);
      calcluz(psi, cbeta);
      calclux(psi_t, cbeta);

      if ((cnti + 1) % iter == 0) {
         calcnorm(&norm, psi, tmpxi, tmpyi, tmpzi);
         calcmuen(&mu, &en, psi, psi_t, dpsi, dpsi_t, psidd2, psidd2fft, tmpxi, tmpyi, tmpzi, tmpxj, tmpyj, tmpzj);
         calcrms(rms, psi, psi_t, tmpxi, tmpyi, tmpzi);
         psi2 = cabs(psi[offsetNx2][Ny2][Nz2]) * cabs(psi[offsetNx2][Ny2][Nz2]);
         MPI_Bcast(&psi2, 1, MPI_DOUBLE, rankNx2, MPI_COMM_WORLD);

         //fprintf(out, "%8le   %8le   %8le   %8le   %8le\n", cnti * dt * par, rms[0], rms[1], rms[2], rms[3]);
         fprintf(out, "NPAS iter.: %5d %10.4f %11.5f %11.5f %10.5f %10.5f\n", step, norm, mu / par, en / par, *rms, psi2);
         fflush(out);

         if (Npasout != NULL) {
            char itername[10];
            sprintf(itername, "-%06d-", step);
            /*
            sprintf(filename, "%s.bin", Npasout);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenxyz(psi, outxyz, mpifile);
            MPI_File_close(&mpifile);
            */

            sprintf(filename, "%s%s1d_x.bin", Npasout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenx(psi, outx, *tmpyi, *tmpzi, mpifile);
            MPI_File_close(&mpifile);

            sprintf(filename, "%s%s1d_y.bin", Npasout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdeny(psi_t, outy, *tmpxi, *tmpzi, mpifile);
            MPI_File_close(&mpifile);

            sprintf(filename, "%s%s1d_z.bin", Npasout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenz(psi, outz, *tmpxi, *tmpyi, mpifile);
            MPI_File_close(&mpifile);

            sprintf(filename, "%s%s2d_xy.bin", Npasout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenxy(psi, outxy, *tmpzi, mpifile);
            MPI_File_close(&mpifile);

            sprintf(filename, "%s%s2d_xz.bin", Npasout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenxz(psi, outxz, *tmpyi, mpifile);
            MPI_File_close(&mpifile);

            sprintf(filename, "%s%s2d_yz.bin", Npasout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outdenyz(psi_t, outyz, *tmpxi, mpifile);
            MPI_File_close(&mpifile);

            sprintf(filename, "%s%s3d_xy0.bin", Npasout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outpsi2xy(psi, outxy, mpifile);
            MPI_File_close(&mpifile);

            sprintf(filename, "%s%s3d_x0z.bin", Npasout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outpsi2xz(psi, outxz, mpifile);
            MPI_File_close(&mpifile);

            sprintf(filename, "%s%s3d_0yz.bin", Npasout, itername);
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
            outpsi2yz(psi, outyz, mpifile);
            MPI_File_close(&mpifile);
        }

        if (rmsout != NULL) {
            fprintf(filerms, "NPAS iter.: %5d %12.5f %12.5f %12.5f %12.5f\n", step, rms[0], rms[1], rms[2], rms[3]);
            fflush(filerms);
        }

        step ++;
      }
   }

   gettimeofday(&iter_stop, NULL);
   iter_time += (double) (((iter_stop.tv_sec - iter_start.tv_sec) * 1000 + (iter_stop.tv_usec - iter_start.tv_usec)/1000.0) + 0.5);

   if (rmsout != NULL) {
      fprintf(filerms, "           ------------------------------------------------------------\n");
      fclose(filerms);
   }

   fprintf(out, "            ---------------------------------------------------------------------\n\n");

   if (finalout != NULL) {
       sprintf(filename, "%s.bin", finalout);
       MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifile);
       outdenxyz(psi, outxyz, mpifile);
       MPI_File_close(&mpifile);
   }

   // Free all dynamically allocated memory. ------------------------

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

   free_double_tensor(dpsi);
   free_double_tensor(dpsi_t);

   free_complex_vector(calphax);
   free_complex_vector(calphay);
   free_complex_vector(calphaz);
   free_complex_matrix(cbeta);
   free_complex_vector(cgammax);
   free_complex_vector(cgammay);
   free_complex_vector(cgammaz);

   free_double_matrix(tmpxi);
   free_double_matrix(tmpyi);
   free_double_matrix(tmpzi);
   free_double_matrix(tmpxj);
   free_double_matrix(tmpyj);
   free_double_matrix(tmpzj);

   fftw_destroy_plan(plan_forward);
   fftw_destroy_plan(plan_backward);
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
   fprintf(out, " Initialization/allocation wall-clock time: %.3f seconds\n", init_time / 1000.);
   fprintf(out, " Calculation (iterations) wall-clock time: %.3f seconds\n", iter_time / 1000.);

   if(output != NULL) fclose(out);

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

   if ((cfg_tmp = cfg_read("G0")) == NULL) {
      if ((cfg_tmp = cfg_read("AHO")) == NULL) {
         fprintf(stderr, "AHO is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      aho = atof(cfg_tmp);

      if ((cfg_tmp = cfg_read("AS")) == NULL) {
         fprintf(stderr, "AS is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      as = atof(cfg_tmp);

      g0 = 4. * pi * as * Na * BOHR_RADIUS / aho;
   } else {
      g0 = atof(cfg_tmp);
   }

   if ((cfg_tmp = cfg_read("GDD0")) == NULL) {
      if ((cfg_tmp = cfg_read("AHO")) == NULL) {
         fprintf(stderr, "AHO is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      aho = atof(cfg_tmp);

      if ((cfg_tmp = cfg_read("ADD")) == NULL) {
         fprintf(stderr, "ADD is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      add = atof(cfg_tmp);

      gd0 = 3. * add * Na * BOHR_RADIUS / aho;
   } else {
      gd0 = atof(cfg_tmp);
   }

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

   if ((cfg_tmp = cfg_read("NPAS")) == NULL) {
      fprintf(stderr, "NPAS is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Npas = atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("ITER")) == NULL)
      iter = Npas;
   else iter = Npas / atol(cfg_tmp);

   if ((cfg_tmp = cfg_read("CUTOFF")) == NULL) {
      fprintf(stderr, "CUTOFF is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   cutoff = atof(cfg_tmp);

   input = cfg_read("INPUT");
   output = cfg_read("OUTPUT");
   rmsout = cfg_read("RMSOUT");
   initout = cfg_read("INITOUT");
   Npasout = cfg_read("NPASOUT");
   finalout = cfg_read("FINALOUT");

   if ((initout != NULL) || (Npasout != NULL)) {
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
 *    psi  - array with the wave function values
 *    tmpz - temporary array
 */
void initpsi(double complex ***psi, double *tmpz) {
   long cnti, cntj, cntk;
   double cpsi;
   double tmp;
   MPI_Offset fileoffset;
   MPI_File file = MPI_FILE_NULL;

   cpsi = sqrt(pi * sqrt(pi / (vgamma * vnu * vlambda)));

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
      fileoffset = rank * sizeof(double) * localNx * Ny * Nz;
      MPI_File_open(MPI_COMM_WORLD, input, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

      if (file == MPI_FILE_NULL) {
         printf("Specify the proper input file with the initial wave function.\n");
         MPI_Finalize();
         exit(EXIT_FAILURE);
      }

      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            MPI_File_read_at(file, fileoffset, tmpz, Nz, MPI_DOUBLE, MPI_STATUS_IGNORE);
            for (cntk = 0; cntk < Nz; cntk ++) {
               psi[cnti][cntj][cntk] = sqrt(tmpz[cntk]);
            }
            fileoffset += Nz * sizeof(double);

            /*for (cntk = 0; cntk < Nz; cntk ++) {
               MPI_File_read_at(file, fileoffset, &tmp, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
               psi[cnti][cntj][cntk] = sqrt(tmp[0]);
               fileoffset += sizeof(double);
            }*/
         }
      }
      MPI_File_close(&file);
   }
   else {
      #pragma omp parallel for private(cnti, cntj, cntk, tmp)
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               tmp = exp(- 0.5 * (vgamma * x2[offsetNx + cnti] + vnu * y2[cntj] + vlambda * z2[cntk]));
               psi[cnti][cntj][cntk] = tmp / cpsi;
            }
         }
      }
   }

   return;
}

/**
 *    Initialization of the potential.
 */
void initpot() {
   long cnti, cntj, cntk;
   double vnu2, vlambda2, vgamma2;

   vnu2 = vnu * vnu;
   vlambda2 = vlambda * vlambda;
   vgamma2 = vgamma * vgamma;

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

   Ax0 = 1. + I * dt / dx2 / (3. - par);
   Ay0 = 1. + I * dt / dy2 / (3. - par);
   Az0 = 1. + I * dt / dz2 / (3. - par);

   Ax0r = 1. - I * dt / dx2 / (3. - par);
   Ay0r = 1. - I * dt / dy2 / (3. - par);
   Az0r = 1. - I * dt / dz2 / (3. - par);

   Ax = - 0.5 * I * dt / dx2 / (3. - par);
   Ay = - 0.5 * I * dt / dy2 / (3. - par);
   Az = - 0.5 * I * dt / dz2 / (3. - par);

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
            potdd[cntj][cnti][cntk] = (4. * pi * (3. * kz2[cntk] / (kx2[cnti] + ky2[cntj] + kz2[cntk]) - 1.) / 3.) * tmp;
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
   int threadid;
   long cnti, cntj, cntk;
   double tmp;
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
      *norm = sqrt(simpint(dx, *tmpx, Nx));
   }

   MPI_Bcast(norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   /*   tmp = 1. / *norm;

      #pragma omp for private(cnti, cntj, cntk)
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               psi[cnti][cntj][cntk] *= tmp;
            }
         }
      }
   */

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
void calcmuen(double *mu, double *en, double complex ***psi, double complex ***psi_t, double ***dpsi, double ***dpsi_t, double ***psidd2, fftw_complex *psidd2fft, double **tmpxi, double **tmpyi, double **tmpzi, double **tmpxj, double **tmpyj, double **tmpzj) {
   int threadid;
   long cnti, cntj, cntk;
   double psi2, psi2lin, psidd2lin, dpsi2;
   void *sendbuf;

   calcpsidd2(psi, psidd2, psidd2fft);

   fftw_execute(plan_transpose_x);

   #pragma omp parallel private(threadid, cnti, cntj, cntk)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cntj = 0; cntj < localNy; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            for (cnti = 0; cnti < Nx; cnti ++) {
               tmpxi[threadid][cnti] = cabs(psi_t[cnti][cntj][cntk]);
            }
            diff(dx, tmpxi[threadid], tmpxj[threadid], Nx);
            for (cnti = 0; cnti < Nx; cnti ++) {
               dpsi_t[cnti][cntj][cntk] = tmpxj[threadid][cnti] * tmpxj[threadid][cnti];
            }
         }
      }
   }

   fftw_execute(plan_transpose_dpsi);

   #pragma omp parallel private(threadid, cnti, cntj, cntk, psi2, psi2lin, psidd2lin, dpsi2)
   {
      threadid = omp_get_thread_num();

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            for (cntj = 0; cntj < Ny; cntj ++) {
               tmpyi[threadid][cntj] = cabs(psi[cnti][cntj][cntk]);
            }
            diff(dy, tmpyi[threadid], tmpyj[threadid], Ny);
            for (cntj = 0; cntj < Ny; cntj ++) {
               dpsi[cnti][cntj][cntk] += tmpyj[threadid][cntj] * tmpyj[threadid][cntj];
            }
         }
      }

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               tmpzi[threadid][cntk] = cabs(psi[cnti][cntj][cntk]);
            }
            diff(dz, tmpzi[threadid], tmpzj[threadid], Nz);
            for (cntk = 0; cntk < Nz; cntk ++) {
               dpsi[cnti][cntj][cntk] += tmpzj[threadid][cntk] * tmpzj[threadid][cntk];
            }
         }
      }
      #pragma omp barrier

      #pragma omp for
      for (cnti = 0; cnti < localNx; cnti ++) {
         for (cntj = 0; cntj < Ny; cntj ++) {
            for (cntk = 0; cntk < Nz; cntk ++) {
               psi2 = cabs(psi[cnti][cntj][cntk]);
               psi2 *= psi2;
               psi2lin = psi2 * g;
               psidd2lin = psidd2[cnti][cntj][cntk] * gd;
               dpsi2 = dpsi[cnti][cntj][cntk] / (3. - par);
               tmpzi[threadid][cntk] = (pot[cnti][cntj][cntk] + psi2lin + psidd2lin) * psi2 + dpsi2;
               tmpzj[threadid][cntk] = (pot[cnti][cntj][cntk] + 0.5 * psi2lin + 0.5 * psidd2lin) * psi2 + dpsi2;
            }
            tmpyi[threadid][cntj] = simpint(dz, tmpzi[threadid], Nz);
            tmpyj[threadid][cntj] = simpint(dz, tmpzj[threadid], Nz);
         }
         (*tmpxi)[cnti] = simpint(dy, tmpyi[threadid], Ny);
         (*tmpxj)[cnti] = simpint(dy, tmpyj[threadid], Ny);
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);

   sendbuf = (rank == 0) ? MPI_IN_PLACE : *tmpxi;
   MPI_Gather(sendbuf, localNx, MPI_DOUBLE, *tmpxi, localNx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   sendbuf = (rank == 0) ? MPI_IN_PLACE : *tmpxj;
   MPI_Gather(sendbuf, localNx, MPI_DOUBLE, *tmpxj, localNx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if (rank == 0) {
      *mu = simpint(dx, *tmpxi, Nx);
      *en = simpint(dx, *tmpxj, Nx);
   }

   return;
}

/**
 *    Calculation of squared wave function values for dipole-dipole
 *    interaction.
 *    psi       - array with the wave function values
 *    psidd2    - array with the squared wave function values
 *    psidd2fft - array with the squared wave function fft values
 */
void calcpsidd2(double complex***psi, double ***psidd2, fftw_complex *psidd2fft) {
   long cnti, cntj, cntk;
   long last = 0;
   double tmp;
   double *psidd2tmp = (double *) psidd2fft;

   #pragma omp parallel for private(cnti, cntj, cntk, tmp)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            tmp = cabs(psi[cnti][cntj][cntk]);
            psidd2[cnti][cntj][cntk] = tmp * tmp;
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
   double psi2, psi2lin, psidd2lin, tmp;

   calcpsidd2(psi, psidd2, psidd2fft);

   #pragma omp parallel for private(cnti, cntj, cntk, psi2, psi2lin, psidd2lin, tmp)
   for (cnti = 0; cnti < localNx; cnti ++) {
      for (cntj = 0; cntj < Ny; cntj ++) {
         for (cntk = 0; cntk < Nz; cntk ++) {
            psi2 = cabs(psi[cnti][cntj][cntk]);
            psi2 *= psi2;
            psi2lin = psi2 * g;
            psidd2lin = psidd2[cnti][cntj][cntk] * gd;
            tmp = dt * (pot[cnti][cntj][cntk] + psi2lin + psidd2lin);
            psi[cnti][cntj][cntk] *= cexp(- I * tmp);
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
            cbeta[threadid][Nz - 2] = psi[cnti][cntj][Nz - 1];
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
         outxz[cnti / outstpx][cntk / outstpz][2] = cabs(psi[cnti][Ny2][cntk]) * cabs(psi[cnti][Ny2][cntk]);
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
            outyz[cntj / outstpy][cntk / outstpz][2] = cabs(psi[offsetNx2][cntj][cntk]) * cabs(psi[offsetNx2][cntj][cntk]);
         }
      }

      MPI_File_write_at(file, fileoffset, **outyz, (Ny / outstpy) * (Nz / outstpz) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);
   }
}

void outdenxyz(double complex ***psi, double ***outxyz, MPI_File file) {
   long cnti, cntj, cntk;
   MPI_Offset fileoffset;

   // MPI IO returns error if the array is too large. As a workaround, we write just Ny * Nz at a time.

   //fileoffset = rank * sizeof(double) * localNx * Ny * Nz;
   fileoffset = rank * sizeof(double) * localNx * Ny * Nz;

   for (cnti = 0; cnti < localNx; cnti += outstpx) {
      for (cntj = 0; cntj < Ny; cntj += outstpy) {
         for (cntk = 0; cntk < Nz; cntk += outstpz) {
            //outxyz[cnti][cntj][cntk] = cabs(psi[cnti][cntj][cntk]) * cabs(psi[cnti][cntj][cntk]);
            outxyz[0][cntj][cntk] = cabs(psi[cnti][cntj][cntk]) * cabs(psi[cnti][cntj][cntk]);
         }
      }

      MPI_File_write_at_all(file, fileoffset, **outxyz, (Ny / outstpy) * (Nz / outstpz), MPI_DOUBLE, MPI_STATUS_IGNORE);
      fileoffset += (Ny / outstpy) * (Nz / outstpz) * sizeof(double);
   }

   //MPI_File_write_at_all(file, fileoffset, **outxyz, (localNx / outstpx) * (Ny / outstpy) * (Nz / outstpz), MPI_DOUBLE, MPI_STATUS_IGNORE);

   return;
}
