#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <fftw3-mpi.h>

#include "PAN_FFTW3.h"
#include "panphasia_functions.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

extern const int Nbasis;
extern const int irank_p[3][84];

extern size_t descriptor_order;
extern size_t descriptor_kk_limit;
extern size_t descriptor_base_size;

// Sum of the squared weights of all basis functions for one Fourier mode,
// computed on the fly instead of being stored for the whole grid
static FFTW_REAL compute_mode_weighting(const double complex *sph_bessel_coeff, size_t n4dimen, size_t ncopy,
                                        size_t ix_abs, size_t iy, size_t iz)
{
  FFTW_REAL mode_weighting = 0.0;
  for (size_t m = 0; m < ncopy; m++)
  {
    complex weight = sph_bessel_coeff[n4dimen * irank_p[0][m] + ix_abs] *
                     sph_bessel_coeff[n4dimen * irank_p[1][m] + iy] *
                     sph_bessel_coeff[n4dimen * irank_p[2][m] + iz];
    mode_weighting += cabs(weight) * cabs(weight);
  }
  return mode_weighting;
}

typedef struct
{
  int peer;
  int tag;
  size_t offset_planes;
  int count_planes;
} plane_xfer_t;

int PANPHASIA_compute_kspace_field_cropped_(size_t relative_level, ptrdiff_t N0_fourier_grid,
                                            ptrdiff_t N0_target_grid, ptrdiff_t N0_keep_grid,
                                            ptrdiff_t local_0_target, ptrdiff_t local_0_start_target,
                                            FFTW_COMPLEX *return_field)
{

  size_t copy_list[Nbasis];

  int pmax = 6;

  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  ptrdiff_t N0_pan_grid = descriptor_base_size << relative_level;

  if (N0_fourier_grid % N0_pan_grid != 0)
    return (100015);
  if (N0_target_grid <= 0 || N0_target_grid > N0_fourier_grid ||
      N0_keep_grid <= 0 || N0_keep_grid > N0_target_grid)
    return (100016);

  int fdim = N0_fourier_grid / N0_pan_grid;
  size_t nfft_dim = N0_fourier_grid;
  size_t npan_dim = N0_pan_grid;
  size_t ntarget_dim = (size_t)N0_target_grid;
  size_t nkeep_dim = (size_t)N0_keep_grid;
  size_t half_keep = nkeep_dim / 2;
  size_t target_z_stride = ntarget_dim / 2 + 1;

  int USE_SMALL_FFT = (N0_pan_grid % 2 == 0) && (sizeof(PAN_REAL) == sizeof(FFTW_REAL));
  int SHARED_FOUR_PAN_SPACE = USE_SMALL_FFT && (mpi_size == 1);

  if (!USE_SMALL_FFT)
  {
    if (N0_target_grid != N0_fourier_grid)
      return (100017);
    if (mpi_rank == 0)
    {
      printf("PANPHASIA: falling back to zero-stuffed %td^3 FFT (N_P = %td)\n",
             N0_fourier_grid, N0_pan_grid);
    }
  }

  int nsubdivide = SHARED_FOUR_PAN_SPACE ? 1 : 21; //(pmax%2==0)?pmax+1:pmax+2;

  size_t ncopy = (pmax + 1) * (pmax + 2) * (pmax + 3) / 6;

  if (ncopy % nsubdivide != 0)
    return (100010);
  int nchunk = ncopy / nsubdivide;

  int verbose = 1;
  int flag_output_mode = 2;
  int error;
  ptrdiff_t size_to_alloc_fourier = 0;
  ptrdiff_t size_to_alloc_pan = 0;
  ptrdiff_t size_to_alloc_pan_fft = 0;
  ptrdiff_t local_fourier_x_start = 1, local_fourier_x_end = 0;
  FFTW_PLAN output_coeff_forward_plan;

  ////////////////////////////////////////////////////////////////////////////////////

  if (pmax > descriptor_order)
    return (100020);

  for (size_t i = 0; i < Nbasis; i++)
    copy_list[i] = i;

  //printf("Dimensions of FT   (%td,%td,%td)\n",N0_fourier_grid,N0_fourier_grid,N0_fourier_grid);
  //printf("Dimensions of PG   (%td,%td,%td)\n",N0_pan_grid,N0_pan_grid,N0_pan_grid);
  //printf("local_no %td local_0_start_fourier %td\n",local_n0_fourier_return, local_0_start_fourier_return);

  //  Compute 1-D Spherical Bessel coefficients for each order //////////////////
  //  These are needed for the convolutions below              //////////////////
  size_t n4dimen;

  n4dimen = (nfft_dim % 4 == 0) ? 4 * (nfft_dim / 4) + 4 : 4 * (nfft_dim / 4) + 5;

  double complex *sph_bessel_coeff = FFTW_MALLOC(sizeof(double complex) * n4dimen * (pmax + 1));

  if (sph_bessel_coeff == NULL)
    return (100030);

  compute_sph_bessel_coeffs(nfft_dim, pmax, n4dimen, fdim, sph_bessel_coeff);

  //printf("Reached here! ndimen4 %ld\n",n4dimen);
  ///////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  // Determine sizes of Fourier and Panphasia coefficient 3-D arrays //

  ptrdiff_t local_n0_pan;
  ptrdiff_t local_0_start_pan;
  ptrdiff_t local_n1_pan;
  ptrdiff_t local_1_start_pan;
  {

    int rank = 3;

    if (USE_SMALL_FFT)
    {
      const ptrdiff_t ndimens_target[] = {N0_target_grid, N0_target_grid, N0_target_grid / 2 + 1};
      ptrdiff_t t_n0, t_s0, t_n1, t_s1;
      FFTW_MPI_LOCAL_SIZE_MANY_TRANSPOSED(rank, ndimens_target, 1,
                                          FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                          MPI_COMM_WORLD, &t_n0, &t_s0, &t_n1, &t_s1);
      if (t_s1 != local_0_start_target)
      {
        printf("Error local_1_start_target mismatch\n");
        return (100032);
      }
      if (t_n1 != local_0_target)
      {
        printf("Error local_1_size_target mismatch\n");
        return (100033);
      }

      const ptrdiff_t ndimens_alloc_pan[] = {N0_pan_grid, N0_pan_grid, N0_pan_grid + 2};
      ptrdiff_t local_n0_pan_fft, local_0_start_pan_fft;
      ptrdiff_t local_n1_pan_fft, local_1_start_pan_fft;
      if (SHARED_FOUR_PAN_SPACE)
      {
        size_to_alloc_pan_fft = FFTW_MPI_LOCAL_SIZE_MANY(rank, ndimens_alloc_pan, nchunk,
                                                         FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD,
                                                         &local_n0_pan_fft, &local_0_start_pan_fft);
        if ((local_0_start_pan_fft != 0) || (local_n0_pan_fft != N0_pan_grid))
          return (100037);
        local_0_start_pan = 0;
        local_n0_pan = N0_pan_grid;
        local_1_start_pan = 0;
        local_n1_pan = N0_pan_grid;
        size_to_alloc_pan = size_to_alloc_pan_fft;
      }
      else
      {
        size_to_alloc_pan_fft = FFTW_MPI_LOCAL_SIZE_MANY_TRANSPOSED(
            rank, ndimens_alloc_pan, nchunk,
            FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD,
            &local_n0_pan_fft, &local_0_start_pan_fft,
            &local_n1_pan_fft, &local_1_start_pan_fft);
        local_0_start_pan = local_0_start_pan_fft;
        local_n0_pan = local_n0_pan_fft;
        local_1_start_pan = local_1_start_pan_fft;
        local_n1_pan = local_n1_pan_fft;
        size_to_alloc_pan = (ptrdiff_t)ncopy * local_n0_pan * N0_pan_grid * (N0_pan_grid + 2);
      }
    }
    else
    {
      ptrdiff_t local_n0_fourier;
      ptrdiff_t local_0_start_fourier;
      const ptrdiff_t ndimens_alloc_fourier[] = {N0_fourier_grid, N0_fourier_grid, N0_fourier_grid + 2};
      ptrdiff_t howmany = nchunk;
      size_to_alloc_fourier = FFTW_MPI_LOCAL_SIZE_MANY(rank, ndimens_alloc_fourier, howmany,
                                                       FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD,
                                                       &local_n0_fourier, &local_0_start_fourier);

      if (local_0_start_fourier != local_0_start_target)
      {
        printf("Error local_0_start_fourier!=local_0_start_fourier_return\n");
        return (100032);
      };

      if (local_n0_fourier != local_0_target)
      {
        printf("Error local_n0_fourier!=local_n0_fourier_return\n");
        return (100033);
      };

      ptrdiff_t abs_fourier_x_start, abs_fourier_x_end;

      if (local_0_start_target % fdim == 0)
      {
        abs_fourier_x_start = local_0_start_target;
      }
      else
      {
        abs_fourier_x_start = local_0_start_target + fdim - (local_0_start_target % fdim);
      };

      if ((local_0_start_target + local_0_target - 1) % fdim == 0)
      {
        abs_fourier_x_end = local_0_start_target + local_0_target - 1;
      }
      else
      {
        abs_fourier_x_end = (local_0_start_target + local_0_target - 1) -
                            ((local_0_start_target + local_0_target - 1) % fdim);
      };

      if (local_0_target <= 0 || abs_fourier_x_end < abs_fourier_x_start)
      {
        local_0_start_pan = 0;
        local_n0_pan = 0;
        local_fourier_x_start = 1;
        local_fourier_x_end = 0;
      }
      else
      {
        if ((abs_fourier_x_end - abs_fourier_x_start) % fdim != 0)
          return (100036);

        local_0_start_pan = abs_fourier_x_start / fdim;
        local_n0_pan = 1 + (abs_fourier_x_end - abs_fourier_x_start) / fdim;

        local_fourier_x_start = abs_fourier_x_start - local_0_start_target;
        local_fourier_x_end = abs_fourier_x_end - local_0_start_target;
      }

      local_1_start_pan = local_0_start_pan;
      local_n1_pan = local_n0_pan;
      howmany = ncopy;
      size_to_alloc_pan = howmany * local_n0_pan * N0_pan_grid * (N0_pan_grid + 2);
    };

    //printf("size_to_alloc_fourier = %td\n",size_to_alloc_fourier);
    //printf("size_to_alloc_pan     = %td\n",size_to_alloc_pan);
  };
  /////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Allocate array to store Panphasia coefficients and compute them before allocating
  // FFT/exchange workspace, so the octree workspace is freed before FFT buffers are allocated.

  void *panphasia_coefficients = FFTW_MALLOC(sizeof(PAN_REAL) * (size_to_alloc_pan > 0 ? size_to_alloc_pan : 1));
  void *fourier_grids = NULL;

  if (panphasia_coefficients == NULL)
    return (100040);

  FFTW_REAL *ptr_real_fourier_grid;
  FFTW_REAL *ptr_panphasia_coefficients = panphasia_coefficients;

  FFTW_COMPLEX *ptr_cmplx_fourier_grid;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////// Compute the Panphasia coefficients   ////////////////////////////////

  if (local_n0_pan > 0)
  {

    size_t xorigin = local_0_start_pan, yorigin = 0, zorigin = 0;

    size_t xextent = local_n0_pan, yextent = N0_pan_grid, zextent = N0_pan_grid;

    if ((error = PANPHASIA_compute_coefficients_(&xorigin, &yorigin, &zorigin, &xextent, &yextent,
                                                 &zextent, copy_list, &ncopy, panphasia_coefficients,
                                                 &flag_output_mode, &verbose)))
      return (100100 + error);
  };

  if (SHARED_FOUR_PAN_SPACE)
  {
    ptr_real_fourier_grid = panphasia_coefficients;
    ptr_cmplx_fourier_grid = panphasia_coefficients;
  }
  else
  {
    ptrdiff_t fft_alloc = USE_SMALL_FFT ? size_to_alloc_pan_fft : size_to_alloc_fourier;
    fourier_grids = FFTW_MALLOC(sizeof(FFTW_REAL) * (fft_alloc > 0 ? fft_alloc : 2));

    if (fourier_grids == NULL)
      return (100041);

    ptr_real_fourier_grid = fourier_grids;
    ptr_cmplx_fourier_grid = fourier_grids;
  };
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////

  //////////  Output diagnostics for small grids only

  {

    if (N0_pan_grid < -1)
    {

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      char filename[100];
      sprintf(filename, "output_real_space_field.%d", rank);

      FILE *fp;
      printf("local_n0_pan %ld npan_dim %ld N0_pan_grid %ld\n", local_n0_pan,
             npan_dim, N0_pan_grid);

      fp = fopen(filename, "w");

      for (int ix = 0; ix < local_n0_pan; ix++)
        for (int iy = 0; iy < npan_dim; iy++)
          for (int iz = 0; iz < npan_dim; iz++)
          {
            int index = ix * N0_pan_grid * (N0_pan_grid + 2) + iy * (N0_pan_grid + 2) + iz;
            fprintf(fp, "%6lu%6d%6d %14.8lf %d\n", ix + local_0_start_pan, iy, iz,
                    ptr_panphasia_coefficients[index], index);
          };

      fclose(fp);
    };
  };

  //----------------------------------------------------------------------------------
  ////////////// Set up FTTW plan //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  // printf("Making the plan ... \n");

  //////////////////// Make plan for ncopy interleaved FTs ///////////////////////////
  //////////////////////////////////////////////////////////////////////////

  {
    int rank = 3;
    ptrdiff_t n_fft_plan = USE_SMALL_FFT ? N0_pan_grid : N0_fourier_grid;
    const ptrdiff_t ndimens[3] = {n_fft_plan, n_fft_plan, n_fft_plan};
    ptrdiff_t howmany = nchunk;
    ptrdiff_t block = FFTW_MPI_DEFAULT_BLOCK;
    ptrdiff_t tblock = FFTW_MPI_DEFAULT_BLOCK;
    unsigned flags = FFTW_ESTIMATE;
    if (USE_SMALL_FFT && !SHARED_FOUR_PAN_SPACE)
      flags |= FFTW_MPI_TRANSPOSED_OUT;

    output_coeff_forward_plan = FFTW_MPI_PLAN_MANY_DTF_R2C(rank, ndimens,
                                                           howmany, block, tblock,
                                                           ptr_real_fourier_grid, ptr_cmplx_fourier_grid,
                                                           MPI_COMM_WORLD, flags);
    if (output_coeff_forward_plan == NULL)
    {
      printf("Null plan\n");
      return (100051);
    };
  };

  // Precompute kept mode indices on the target grid (skipping modes discarded by
  // FourierInterpolateCopyTo as well as the target Nyquist planes)
  size_t *kept_y_local = NULL, *kept_y_f = NULL;
  size_t *kept_x_local = NULL, *kept_x_f = NULL;
  size_t n_kept_y = 0, n_kept_x = 0;

  if (USE_SMALL_FFT)
  {
    size_t alloc_y = (size_t)(local_0_target > 0 ? local_0_target : 1);
    size_t alloc_x = ntarget_dim > 0 ? ntarget_dim : 1;
    kept_y_local = malloc(sizeof(size_t) * alloc_y);
    kept_y_f = malloc(sizeof(size_t) * alloc_y);
    kept_x_local = malloc(sizeof(size_t) * alloc_x);
    kept_x_f = malloc(sizeof(size_t) * alloc_x);
    if (!kept_y_local || !kept_y_f || !kept_x_local || !kept_x_f)
      return (100047);

    for (ptrdiff_t i = 0; i < local_0_target; i++)
    {
      size_t iy_t = (size_t)(i + local_0_start_target);
      if (iy_t < half_keep)
      {
        kept_y_local[n_kept_y] = (size_t)i;
        kept_y_f[n_kept_y] = iy_t;
        n_kept_y++;
      }
      else if (iy_t > ntarget_dim - half_keep)
      {
        kept_y_local[n_kept_y] = (size_t)i;
        kept_y_f[n_kept_y] = iy_t + nfft_dim - ntarget_dim;
        n_kept_y++;
      }
    }

    for (size_t j = 0; j < ntarget_dim; j++)
    {
      if (j < half_keep)
      {
        kept_x_local[n_kept_x] = j;
        kept_x_f[n_kept_x] = j;
        n_kept_x++;
      }
      else if (j > ntarget_dim - half_keep)
      {
        kept_x_local[n_kept_x] = j;
        kept_x_f[n_kept_x] = j + nfft_dim - ntarget_dim;
        n_kept_x++;
      }
    }
  }

  // Set up plane pointers and MPI plane exchange schedule for the N_P FFT path
  size_t plane_cmplx_count = (size_t)nchunk * npan_dim * (npan_dim / 2 + 1);
  FFTW_COMPLEX **plane_ptr = NULL;
  FFTW_COMPLEX *recv_buf = NULL;
  plane_xfer_t *recv_runs = NULL;
  plane_xfer_t *send_runs = NULL;
  size_t n_recv_runs = 0, n_send_runs = 0;
  MPI_Request *xfer_reqs = NULL;
  MPI_Datatype row_mpi_type = MPI_DATATYPE_NULL, plane_mpi_type = MPI_DATATYPE_NULL;
  int have_mpi_types = 0;

  if (USE_SMALL_FFT && !SHARED_FOUR_PAN_SPACE)
  {
    plane_ptr = calloc(npan_dim, sizeof(FFTW_COMPLEX *));
    if (plane_ptr == NULL)
      return (100042);

    for (ptrdiff_t lp = 0; lp < local_n1_pan; lp++)
    {
      plane_ptr[local_1_start_pan + lp] = ptr_cmplx_fourier_grid + (size_t)lp * plane_cmplx_count;
    }

    long long my_start = (long long)local_1_start_pan;
    long long my_count = (long long)local_n1_pan;
    long long *pan_starts = malloc(sizeof(long long) * (size_t)mpi_size);
    long long *pan_counts = malloc(sizeof(long long) * (size_t)mpi_size);
    unsigned char *need_remote = calloc(npan_dim, sizeof(unsigned char));
    unsigned char *all_need = malloc((size_t)mpi_size * npan_dim);
    size_t *recv_slot_of_plane = malloc(sizeof(size_t) * npan_dim);

    if (!pan_starts || !pan_counts || !need_remote || !all_need || !recv_slot_of_plane)
      return (100043);

    MPI_Allgather(&my_start, 1, MPI_LONG_LONG, pan_starts, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&my_count, 1, MPI_LONG_LONG, pan_counts, 1, MPI_LONG_LONG, MPI_COMM_WORLD);

    int need_conj_z = (half_keep > npan_dim / 2 + 1);
    if (n_kept_x > 0 && half_keep > 0)
    {
      for (size_t yi = 0; yi < n_kept_y; yi++)
      {
        size_t iy_p = kept_y_f[yi] % npan_dim;
        need_remote[iy_p] = 1;
        if (need_conj_z)
        {
          size_t iy_p_conj = (iy_p == 0) ? 0 : (npan_dim - iy_p);
          need_remote[iy_p_conj] = 1;
        }
      }
    }

    for (ptrdiff_t lp = 0; lp < local_n1_pan; lp++)
      need_remote[local_1_start_pan + lp] = 0;

    MPI_Allgather(need_remote, (int)npan_dim, MPI_UNSIGNED_CHAR,
                  all_need, (int)npan_dim, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);

    recv_runs = malloc(sizeof(plane_xfer_t) * (npan_dim > 0 ? npan_dim : 1));
    send_runs = malloc(sizeof(plane_xfer_t) * ((size_t)mpi_size * (size_t)(local_n1_pan > 0 ? local_n1_pan : 1)));
    if (!recv_runs || !send_runs)
      return (100044);

    size_t total_recv_planes = 0;
    for (int owner = 0; owner < mpi_size; owner++)
    {
      if (owner == mpi_rank || pan_counts[owner] <= 0)
        continue;
      size_t s0 = (size_t)pan_starts[owner];
      size_t s1 = s0 + (size_t)pan_counts[owner];
      int tag = 0;
      size_t p = s0;
      while (p < s1)
      {
        while (p < s1 && !need_remote[p])
          p++;
        if (p >= s1)
          break;
        size_t p_start = p;
        while (p < s1 && need_remote[p])
        {
          recv_slot_of_plane[p] = total_recv_planes + (p - p_start);
          p++;
        }
        int len = (int)(p - p_start);
        recv_runs[n_recv_runs].peer = owner;
        recv_runs[n_recv_runs].tag = tag++;
        recv_runs[n_recv_runs].offset_planes = total_recv_planes;
        recv_runs[n_recv_runs].count_planes = len;
        n_recv_runs++;
        total_recv_planes += (size_t)len;
      }
    }

    for (int dest = 0; dest < mpi_size; dest++)
    {
      if (dest == mpi_rank || local_n1_pan <= 0)
        continue;
      const unsigned char *dest_need = all_need + (size_t)dest * npan_dim;
      int tag = 0;
      ptrdiff_t lp = 0;
      while (lp < local_n1_pan)
      {
        while (lp < local_n1_pan && !dest_need[local_1_start_pan + lp])
          lp++;
        if (lp >= local_n1_pan)
          break;
        ptrdiff_t lp_start = lp;
        while (lp < local_n1_pan && dest_need[local_1_start_pan + lp])
          lp++;
        send_runs[n_send_runs].peer = dest;
        send_runs[n_send_runs].tag = tag++;
        send_runs[n_send_runs].offset_planes = (size_t)lp_start;
        send_runs[n_send_runs].count_planes = (int)(lp - lp_start);
        n_send_runs++;
      }
    }

    if (total_recv_planes > 0)
    {
      recv_buf = FFTW_MALLOC(sizeof(FFTW_COMPLEX) * total_recv_planes * plane_cmplx_count);
      if (recv_buf == NULL)
        return (100045);
      for (size_t p = 0; p < npan_dim; p++)
      {
        if (need_remote[p])
          plane_ptr[p] = recv_buf + recv_slot_of_plane[p] * plane_cmplx_count;
      }
    }

    size_t total_reqs = n_recv_runs + n_send_runs;
    xfer_reqs = malloc(sizeof(MPI_Request) * (total_reqs > 0 ? total_reqs : 1));
    if (xfer_reqs == NULL)
      return (100046);

    MPI_Datatype real_mpi_type = (sizeof(FFTW_REAL) == sizeof(double)) ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Type_contiguous((int)(2 * (size_t)nchunk * (npan_dim / 2 + 1)), real_mpi_type, &row_mpi_type);
    MPI_Type_commit(&row_mpi_type);
    MPI_Type_contiguous((int)npan_dim, row_mpi_type, &plane_mpi_type);
    MPI_Type_commit(&plane_mpi_type);
    have_mpi_types = 1;

    free(recv_slot_of_plane);
    free(all_need);
    free(need_remote);
    free(pan_counts);
    free(pan_starts);
  }

  //printf("Plan completed ... \n");

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //----------------------------------------------------------------------------------

  {
    size_t total_target_cmplx = (size_t)local_0_target * ntarget_dim * target_z_stride;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t idx = 0; idx < total_target_cmplx; idx++)
      return_field[idx] = 0.0;
  }

  for (int iter = 0; iter < nsubdivide; iter++)
  {

    int moffset = iter * nchunk;

    if (!SHARED_FOUR_PAN_SPACE)
    {
      if (USE_SMALL_FFT)
      {
#ifdef USE_OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (ptrdiff_t ix_p = 0; ix_p < local_n0_pan; ix_p++)
          for (size_t iy_p = 0; iy_p < npan_dim; iy_p++)
          {
            size_t row_base = (size_t)ix_p * npan_dim * (npan_dim + 2) + iy_p * (npan_dim + 2);
            for (size_t iz_p = 0; iz_p < npan_dim; iz_p++)
            {
              size_t idx = row_base + iz_p;
              for (int m = 0; m < nchunk; m++)
              {
                ptr_real_fourier_grid[nchunk * idx + m] =
                    ptr_panphasia_coefficients[ncopy * idx + moffset + m];
              }
            }
            for (size_t iz_p = npan_dim; iz_p < npan_dim + 2; iz_p++)
            {
              size_t idx = row_base + iz_p;
              for (int m = 0; m < nchunk; m++)
                ptr_real_fourier_grid[nchunk * idx + m] = 0.0;
            }
          }
      }
      else
      {
        memset(ptr_real_fourier_grid, 0, sizeof(FFTW_REAL) * size_to_alloc_fourier);

        // Copy Panphasia coefficients to Fourier grid with appropriate stride - fdim

        size_t index_p, index_f;
        int m;

        for (int ix_f = local_fourier_x_start, ix_p = 0; ix_f <= local_fourier_x_end; ix_f += fdim, ix_p++)
#ifdef USE_OPENMP
#pragma omp parallel for collapse(2) private(index_p, index_f, m)
#endif

          for (size_t iy_p = 0; iy_p < npan_dim; iy_p++)
            for (size_t iz_p = 0; iz_p < npan_dim; iz_p++)
            {
              size_t iy_f = iy_p * fdim;
              size_t iz_f = iz_p * fdim;
              index_p = ix_p * N0_pan_grid * (N0_pan_grid + 2) + iy_p * (N0_pan_grid + 2) + iz_p;
              index_f = ix_f * N0_fourier_grid * (N0_fourier_grid + 2) + iy_f * (N0_fourier_grid + 2) + iz_f;

              for (m = 0; m < nchunk; m++)
              {
                ptr_real_fourier_grid[nchunk * index_f + m] =
                    ptr_panphasia_coefficients[ncopy * index_p + moffset + m];
              };
            };
      }
    };

    FFTW_MPI_EXECUTE_DFT_R2C(output_coeff_forward_plan,
                             ptr_real_fourier_grid, ptr_cmplx_fourier_grid);

    if (USE_SMALL_FFT && !SHARED_FOUR_PAN_SPACE)
    {
      int nreq = 0;
      for (size_t r = 0; r < n_recv_runs; r++)
      {
        MPI_Irecv(recv_buf + recv_runs[r].offset_planes * plane_cmplx_count,
                  recv_runs[r].count_planes, plane_mpi_type,
                  recv_runs[r].peer, recv_runs[r].tag, MPI_COMM_WORLD, &xfer_reqs[nreq++]);
      }
      for (size_t s = 0; s < n_send_runs; s++)
      {
        MPI_Isend(ptr_cmplx_fourier_grid + send_runs[s].offset_planes * plane_cmplx_count,
                  send_runs[s].count_planes, plane_mpi_type,
                  send_runs[s].peer, send_runs[s].tag, MPI_COMM_WORLD, &xfer_reqs[nreq++]);
      }
      if (nreq > 0)
        MPI_Waitall(nreq, xfer_reqs, MPI_STATUSES_IGNORE);
    }

    if (USE_SMALL_FFT)
    {
      // Convolve and combine only kept modes directly into transposed target layout
#ifdef USE_OPENMP
#pragma omp parallel for collapse(3)
#endif
      for (size_t yi = 0; yi < n_kept_y; yi++)
        for (size_t xi = 0; xi < n_kept_x; xi++)
          for (size_t iz = 0; iz < half_keep; iz++)
          {
            size_t i_loc = kept_y_local[yi];
            size_t iy_f = kept_y_f[yi];
            size_t j_loc = kept_x_local[xi];
            size_t ix_f = kept_x_f[xi];
            size_t index1 = (i_loc * ntarget_dim + j_loc) * target_z_stride + iz;

            size_t ix_p = ix_f % npan_dim;
            size_t iy_p = iy_f % npan_dim;
            size_t iz_p = iz % npan_dim;
            int use_conj = 0;
            if (iz_p > npan_dim / 2)
            {
              ix_p = (ix_p == 0) ? 0 : npan_dim - ix_p;
              iy_p = (iy_p == 0) ? 0 : npan_dim - iy_p;
              iz_p = npan_dim - iz_p;
              use_conj = 1;
            }

            const FFTW_COMPLEX *src_cmplx;
            size_t base_index2;
            if (SHARED_FOUR_PAN_SPACE)
            {
              src_cmplx = ptr_cmplx_fourier_grid;
              base_index2 = (size_t)nchunk * ((ix_p * npan_dim + iy_p) * (npan_dim / 2 + 1) + iz_p);
            }
            else
            {
              src_cmplx = plane_ptr[iy_p];
              base_index2 = (size_t)nchunk * (ix_p * (npan_dim / 2 + 1) + iz_p);
            }

            FFTW_COMPLEX sum = return_field[index1];
            for (int m = 0; m < nchunk; m++)
            {
              double complex coeff = src_cmplx[base_index2 + (size_t)m];
              if (use_conj)
                coeff = conj(coeff);

              complex weight = sph_bessel_coeff[n4dimen * irank_p[0][m + moffset] + ix_f] *
                               sph_bessel_coeff[n4dimen * irank_p[1][m + moffset] + iy_f] *
                               sph_bessel_coeff[n4dimen * irank_p[2][m + moffset] + iz];
              sum += weight * coeff;
            }
            return_field[index1] = sum;
          }
    }
    else
    {
      // Fallback path: convolve and combine only kept modes on untransposed N_F grid
      size_t index1, index2;
      complex weight;

#ifdef USE_OPENMP
#pragma omp parallel for collapse(3) private(index1, index2, weight)
#endif
      for (int ix = 0; ix < local_0_target; ix++)
        for (size_t iy = 0; iy < nfft_dim; iy++)
          for (size_t iz = 0; iz < half_keep; iz++)
          {
            size_t ix_f = (size_t)(ix + local_0_start_target);
            if (!(ix_f < half_keep || ix_f > nfft_dim - half_keep) ||
                !(iy < half_keep || iy > nfft_dim - half_keep))
              continue;

            index1 = (size_t)ix * nfft_dim * (nfft_dim / 2 + 1) + iy * (nfft_dim / 2 + 1) + iz;
            FFTW_COMPLEX sum = return_field[index1];
            for (int m = 0; m < nchunk; m++)
            {
              index2 = (size_t)nchunk * index1 + (size_t)m;
              double complex coeff = ptr_cmplx_fourier_grid[index2];
              weight = sph_bessel_coeff[n4dimen * irank_p[0][m + moffset] + ix_f] *
                       sph_bessel_coeff[n4dimen * irank_p[1][m + moffset] + iy] *
                       sph_bessel_coeff[n4dimen * irank_p[2][m + moffset] + iz];
              sum += weight * coeff;
            }
            return_field[index1] = sum;
          }
    }

  }; // End loop over iter

  // Free coefficient and FFT/exchange buffers before normalisation
  FFTW_DESTROY_PLAN(output_coeff_forward_plan);
  if (!SHARED_FOUR_PAN_SPACE)
    FFTW_FREE(fourier_grids);
  if (recv_buf != NULL)
    FFTW_FREE(recv_buf);
  FFTW_FREE(panphasia_coefficients);

  free(plane_ptr);
  free(recv_runs);
  free(send_runs);
  free(xfer_reqs);
  if (have_mpi_types)
  {
    MPI_Type_free(&plane_mpi_type);
    MPI_Type_free(&row_mpi_type);
  }

  // Add phase shift and normalise kept modes (skipped modes are already zero)
  {
    const double pi = 4.0 * atan(1.0);

    if (USE_SMALL_FFT)
    {
#ifdef USE_OPENMP
#pragma omp parallel for collapse(3)
#endif
      for (size_t yi = 0; yi < n_kept_y; yi++)
        for (size_t xi = 0; xi < n_kept_x; xi++)
          for (size_t iz = 0; iz < half_keep; iz++)
          {
            size_t i_loc = kept_y_local[yi];
            size_t iy_f = kept_y_f[yi];
            size_t j_loc = kept_x_local[xi];
            size_t ix_f = kept_x_f[xi];
            size_t index1 = (i_loc * ntarget_dim + j_loc) * target_z_stride + iz;

            int kx = (ix_f > nfft_dim / 2) ? (int)ix_f - (int)nfft_dim : (int)ix_f;
            int ky = (iy_f > nfft_dim / 2) ? (int)iy_f - (int)nfft_dim : (int)iy_f;
            int kz = (int)iz;

            double complex phase_shift_and_scale;
            if ((ix_f == nfft_dim / 2) || (iy_f == nfft_dim / 2) || (iz == nfft_dim / 2))
            {
              phase_shift_and_scale = 0.0;
            }
            else
            {
              FFTW_REAL mode_weighting = compute_mode_weighting(sph_bessel_coeff, n4dimen, ncopy,
                                                                ix_f, iy_f, iz);
              phase_shift_and_scale = sqrt((double)(fdim * fdim * fdim)) *
                                      cexp((double)fdim * (-I) * pi * (double)(kx + ky + kz) / sqrt(mode_weighting) /
                                           (double)nfft_dim) /
                                      pow((double)nfft_dim, 1.5);
            }
            return_field[index1] *= phase_shift_and_scale;
          }
    }
    else
    {
#ifdef USE_OPENMP
#pragma omp parallel for collapse(3)
#endif
      for (int ix = 0; ix < local_0_target; ix++)
        for (size_t iy = 0; iy < nfft_dim; iy++)
          for (size_t iz = 0; iz < half_keep; iz++)
          {
            size_t ix_f = (size_t)(ix + local_0_start_target);
            if (!(ix_f < half_keep || ix_f > nfft_dim - half_keep) ||
                !(iy < half_keep || iy > nfft_dim - half_keep))
              continue;

            size_t index1 = (size_t)ix * nfft_dim * (nfft_dim / 2 + 1) + iy * (nfft_dim / 2 + 1) + iz;
            int kx = (ix_f > nfft_dim / 2) ? (int)ix_f - (int)nfft_dim : (int)ix_f;
            int ky = (iy > nfft_dim / 2) ? (int)iy - (int)nfft_dim : (int)iy;
            int kz = (int)iz;

            double complex phase_shift_and_scale;
            if ((ix_f == nfft_dim / 2) || (iy == nfft_dim / 2) || (iz == nfft_dim / 2))
            {
              phase_shift_and_scale = 0.0;
            }
            else
            {
              FFTW_REAL mode_weighting = compute_mode_weighting(sph_bessel_coeff, n4dimen, ncopy,
                                                                ix_f, iy, iz);
              phase_shift_and_scale = sqrt((double)(fdim * fdim * fdim)) *
                                      cexp((double)fdim * (-I) * pi * (double)(kx + ky + kz) / sqrt(mode_weighting) /
                                           (double)nfft_dim) /
                                      pow((double)nfft_dim, 1.5);
            }
            return_field[index1] *= phase_shift_and_scale;
          }
    }
  }

  // Rescale selected Fourier modes to unit amplitude when KK limit is set
  if (descriptor_kk_limit > 0)
  {
    if (USE_SMALL_FFT)
    {
#ifdef USE_OPENMP
#pragma omp parallel for collapse(3)
#endif
      for (size_t yi = 0; yi < n_kept_y; yi++)
        for (size_t xi = 0; xi < n_kept_x; xi++)
          for (size_t iz = 0; iz < half_keep; iz++)
          {
            size_t iy_f = kept_y_f[yi];
            size_t ix_f = kept_x_f[xi];
            if ((ix_f == nfft_dim / 2) || (iy_f == nfft_dim / 2) || (iz == nfft_dim / 2))
              continue;
            int kx = (ix_f > nfft_dim / 2) ? (int)ix_f - (int)nfft_dim : (int)ix_f;
            int ky = (iy_f > nfft_dim / 2) ? (int)iy_f - (int)nfft_dim : (int)iy_f;
            int kz = (int)iz;
            size_t ksquared = (size_t)(kx * kx + ky * ky + kz * kz);
            if ((ksquared <= descriptor_kk_limit) && (ksquared != 0))
            {
              size_t index1 = (kept_y_local[yi] * ntarget_dim + kept_x_local[xi]) * target_z_stride + iz;
              complex weight = cabs(return_field[index1]);
              return_field[index1] /= weight;
            }
          }
    }
    else
    {
#ifdef USE_OPENMP
#pragma omp parallel for collapse(3)
#endif
      for (int ix = 0; ix < local_0_target; ix++)
        for (size_t iy = 0; iy < nfft_dim; iy++)
          for (size_t iz = 0; iz < half_keep; iz++)
          {
            size_t ix_f = (size_t)(ix + local_0_start_target);
            if (!(ix_f < half_keep || ix_f > nfft_dim - half_keep) ||
                !(iy < half_keep || iy > nfft_dim - half_keep))
              continue;
            if ((ix_f == nfft_dim / 2) || (iy == nfft_dim / 2) || (iz == nfft_dim / 2))
              continue;
            int kx = (ix_f > nfft_dim / 2) ? (int)ix_f - (int)nfft_dim : (int)ix_f;
            int ky = (iy > nfft_dim / 2) ? (int)iy - (int)nfft_dim : (int)iy;
            int kz = (int)iz;
            size_t ksquared = (size_t)(kx * kx + ky * ky + kz * kz);
            if ((ksquared <= descriptor_kk_limit) && (ksquared != 0))
            {
              size_t index1 = (size_t)ix * nfft_dim * (nfft_dim / 2 + 1) + iy * (nfft_dim / 2 + 1) + iz;
              complex weight = cabs(return_field[index1]);
              return_field[index1] /= weight;
            }
          }
    }
  }

  free(kept_y_local);
  free(kept_y_f);
  free(kept_x_local);
  free(kept_x_f);

  // Transpose output field only for the fallback path (USE_SMALL_FFT writes directly in transposed order)
  if (!USE_SMALL_FFT)
  {
    FFTW_PLAN transpose_output_plan;
    unsigned flags = FFTW_ESTIMATE;
    FFTW_REAL *ptr_return_as_real_field = (FFTW_REAL *)return_field;
    ptrdiff_t howmany = N0_fourier_grid + 2;

    transpose_output_plan = FFTW_MPI_PLAN_MANY_TRANSPOSE(N0_fourier_grid, N0_fourier_grid,
                                                         howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                                         ptr_return_as_real_field, ptr_return_as_real_field,
                                                         MPI_COMM_WORLD, flags);

    FFTW_MPI_EXECUTE_R2R(transpose_output_plan, ptr_return_as_real_field, ptr_return_as_real_field);
    FFTW_DESTROY_PLAN(transpose_output_plan);
  }

  FFTW_FREE(sph_bessel_coeff);

  return (0);
}

int PANPHASIA_compute_kspace_field_(size_t relative_level, ptrdiff_t N0_fourier_grid,
                                    ptrdiff_t local_n0_fourier_return, ptrdiff_t local_0_start_fourier_return,
                                    FFTW_COMPLEX *return_field)
{
  return PANPHASIA_compute_kspace_field_cropped_(relative_level, N0_fourier_grid,
                                                 N0_fourier_grid, N0_fourier_grid,
                                                 local_n0_fourier_return, local_0_start_fourier_return,
                                                 return_field);
}

//==========================================================================================
//==========================================================================================
