#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>
#include <math.h>
#include <stdio.h>
#include <qmp.h>

#if defined(HAVE_LAPACK)
#  include <lapack.h>
#elif defined(HAVE_GSL)
#  include <gsl/gsl_linalg.h>
#  include <gsl/gsl_sort_double.h>
#  include <gsl/gsl_sort_vector_double.h>
#else
#  error "no linear algebra library"
#endif

/* macros to reuse workspace */
#define eps_reortho 2e-8
#define n_reortho   3
#define cur_v       (d->work_c_1)
#define cur_Av      (d->work_c_2)

/* XXX note that q(df_inject orthogonalizes vec wrt the vectors already in d->U
   perhaps you want to copy vec before invoking q(df_inject) */
int
q(df_inject)(struct Q(Deflator) *d,
             struct MxM_workspaceF *ws,
             latvec_c vec)
{
  double v_norm2_min = eps_reortho * eps_reortho;
  int i;

  if (d->usize == d->umax)
    return 0;

  if (0 < d->usize) {
    int i_reortho;
    /* reorthogonalize n_reortho times */
    for (i_reortho = n_reortho; i_reortho--; ) {
      latmat_c cur_U = q(latmat_c_submat_col)(d->U, 0, d->usize);
      q(lat_lmH_dot_lv)(d->usize, 
                        cur_U, 
                        vec, 
                        d->zwork);
      q(lat_lm_dot_zv)(d->usize, 
                       cur_U, 
                       d->zwork, 
                       cur_Av);
      q(lat_c_axpy_d)(-1., cur_Av, vec);
    }
  }

  double v_norm2 = q(lat_c_nrm2)(vec);
  if (v_norm2 < v_norm2_min)
    return 0;

  /* normalize the vector */
  q(lat_c_scal_d)(1. / sqrt(v_norm2), vec);
  q(latmat_c_insert_col)(d->U, d->usize, vec);


  /* compute (U^dag . A . vec) */
  latvec_c_linop(cur_Av, vec, ws);
  latmat_c cur_U = q(latmat_c_submat_col)(d->U, 0, d->usize + 1);
  q(lat_lmH_dot_lv)(d->usize + 1, cur_U, cur_Av, d->H + d->usize * d->umax);

  /* copy conjugated column into respective row */
  for (i = 0; i < d->usize; i++) {
    doublecomplex *p1 = d->H + i + d->usize * d->umax,
      *p2 = d->H + d->usize + i * d->umax;
    p2->r =  p1->r;
    p2->i = -p1->i;
  }
  d->H[d->usize * (d->umax + 1)].i = 0.;

  /* finished adding one new vector */
  d->usize ++;
  return 1;
}

void
q(df_rebuild)(struct Q(Deflator) *d)
{
  int i;

  /* compute Cholesky decomposition */
  memcpy(d->C, d->H, d->usize * d->umax * sizeof(d->C[0]));

#if HAVE_LAPACK
  long int usize  = d->usize,
           umax   = d->umax,
           info   = 0;
  char cU = 'U';
  zpotrf_(&cU, &usize, d->C, &umax, &info, 1);
  assert(0 == info);
#elif HAVE_GSL
  gsl_matrix_complex_view gsl_C = gsl_matrix_complex_view_array_with_tda(
            (double *)d->C, d->usize, d->usize, d->umax);
  /* transpose to convert matrix Fortran[row][col] -> C[col][row] */
  CHECK_GSL_STATUS(gsl_matrix_complex_transpose(&gsl_C.matrix));
  CHECK_GSL_STATUS(gsl_linalg_complex_cholesky_decomp(&gsl_C.matrix));
#else
#  error "no linear algebra library"
#endif
    
  /* XXX broadcast Cholesky matrix from the master node
     to keep deflation consistent
     FIXME check how much these matrices are different on other nodes */
  /* FIXME the matrix is broadcast in a single bucket umax*uzise of 
     doublecomplex; (umax-usize)*usize is garbage; is it more efficient to have 
     usize broadcasts? */
  for (i = 0 ; i < d->usize; i++)
    QMP_broadcast((void *)(d->C + i * d->umax), d->usize * sizeof(d->C[0]));
}

int
q(df_postamble)(
        struct Q(State)       *s,
        struct Q(Deflator)    *d,
        struct MxM_workspaceF *ws,
        unsigned int           options)
{
    int i;

    if (NULL == s ||
            NULL == d ||
            d->frozen || 
            d->umax <= d->usize ||
            d->vsize < d->nev)
        return 0;

    int unew = 0,
        i_v = 0;
    long int usize_old = d->usize;
    while ((d->usize < d->umax) && (i_v < d->nev)) {
/* macros to reuse workspace */
        q(latmat_c_get_col)(d->V, i_v, cur_v);
        if (q(df_inject)(d, ws, cur_v)) {
          unew ++;
        }
        i_v++;
    }
    assert(usize_old + unew == d->usize);


    if (options & QOP_MDWF_LOG_EIG_POSTAMBLE) {
        memcpy(d->C, d->H, d->usize * d->umax * sizeof(d->C[0]));
#if HAVE_LAPACK
        {
            long int usize  = d->usize,
                     umax   = d->umax,
                     info   = 0;
            char cU = 'U',
                cN = 'N';
            zheev_(&cN, &cU, &usize, d->C, &umax, d->debug_hevals, 
                   d->debug_zwork, &(d->debug_lwork), d->debug_rwork,
                   &info, 1, 1);
            assert(0 == info);
            for (i = 0; i < usize; i++)
                qf(zprint)(s, "postamble", "U %4d %17.9e",
                           i, d->debug_hevals[i]);
        }
#elif HAVE_GSL
        {
            long int usize  = d->usize;
            gsl_vector_view gsl_hevals = gsl_vector_subvector(
                d->debug_gsl_hevals, 0, d->usize);
            gsl_matrix_complex_view gsl_C =
                gsl_matrix_complex_view_array_with_tda(
                    (double *)d->C, d->usize, d->usize, d->umax);
            /* transpose to convert matrix Fortran[row][col] -> C[col][row] */
            CHECK_GSL_STATUS(gsl_matrix_complex_transpose(&gsl_C.matrix));
            CHECK_GSL_STATUS(gsl_eigen_herm(&gsl_C.matrix, &gsl_hevals.vector,
                                            d->debug_gsl_wkspace));
            gsl_sort_vector(&gsl_hevals.vector);
            for (i = 0; i < usize; i++)
                qf(zprint)(s, "postamble", "U %4d %17.9e", i,
                           gsl_vector_get(&gsl_hevals.vector, i));
        }
#else
#  error "no linear algebra library"
#endif
    }

    q(df_rebuild)(d);
    return unew;
}
