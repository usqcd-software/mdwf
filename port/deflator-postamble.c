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
#define cur_v       (df->work_c_1)
#define cur_Av      (df->work_c_2)

/* XXX note that q(df_inject orthogonalizes vec wrt the vectors already in df->U
   perhaps you want to copy vec before invoking q(df_inject) */
int
q(df_inject)(struct Q(Deflator) *df,
             struct MxM_workspaceF *ws,
             latvec_c vec)
{
  double v_norm2_min = eps_reortho * eps_reortho;
  int i;

  if (df->usize == df->umax)
    return 0;

  if (0 < df->usize) {
    int i_reortho;
    /* reorthogonalize n_reortho times */
    for (i_reortho = n_reortho; i_reortho--; ) {
      latmat_c cur_U = q(latmat_c_submat_col)(df->U, 0, df->usize);
      q(lat_lmH_dot_lv)(df->usize, 
                        cur_U, 
                        vec, 
                        df->zwork);
      q(lat_lm_dot_zv)(df->usize,
                       cur_U, 
                       df->zwork, 
                       cur_Av);
      q(lat_c_axpy_d)(-1., cur_Av, vec);
    }
  }

  double v_norm2 = q(lat_c_nrm2)(vec);
  if (v_norm2 < v_norm2_min)
    return 0;

  /* normalize the vector */
  q(lat_c_scal_d)(1. / sqrt(v_norm2), vec);
  q(latmat_c_insert_col)(df->U, df->usize, vec);


  /* compute (U^dag . A . vec) */
  latvec_c_linop(cur_Av, vec, ws);
  latmat_c cur_U = q(latmat_c_submat_col)(df->U, 0, df->usize + 1);
  q(lat_lmH_dot_lv)(df->usize + 1, cur_U, cur_Av, df->H + df->usize * df->umax);

  /* copy conjugated column into respective row */
  for (i = 0; i < df->usize; i++) {
    doublecomplex *p1 = df->H + i + df->usize * df->umax,
      *p2 = df->H + df->usize + i * df->umax;
    p2->r =  p1->r;
    p2->i = -p1->i;
  }
  df->H[df->usize * (df->umax + 1)].i = 0.;

  /* finished adding one new vector */
  df->usize ++;
  return 1;
}

int
q(df_rebuild)(struct Q(Deflator) *df)
{
    int status;
    int i;

    /* compute Cholesky decomposition */
    memcpy(df->C, df->H, df->usize * df->umax * sizeof(df->C[0]));

#if HAVE_LAPACK
    long int usize  = df->usize,
           umax   = df->umax,
           info   = 0;
    char cU = 'U';
    zpotrf_(&cU, &usize, df->C, &umax, &info, 1);
    assert(0 == info);
#elif HAVE_GSL
    gsl_matrix_complex_view gsl_C = gsl_matrix_complex_view_array_with_tda(
            (double *)df->C, df->usize, df->usize, df->umax);
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
    for (i = 0 ; i < df->usize; i++)
        QMP_broadcast((void *)(df->C + i * df->umax), df->usize * sizeof(df->C[0]));

    if (0 != (status = q(deflator_calc_evals)(df)))
        return status;

    return 0;
}

int
q(df_postamble)(
        struct Q(State)       *s,
        struct Q(Deflator)    *df,
        struct MxM_workspaceF *ws,
        unsigned int           options)
{
    int i;

    if (NULL == df ||
            NULL == df->state ||
            df->df_eigcg.frozen || 
            df->umax <= df->usize ||
            df->df_eigcg.vsize < df->df_eigcg.nev)
        return 0;
    struct q(DeflatorEigcg) *d_e = &(df->df_eigcg);

    int unew = 0,
        i_v = 0;
    long int usize_old = df->usize;
    while ((df->usize < df->umax) && (i_v < d_e->nev)) {
        q(latmat_c_get_col)(d_e->V, i_v, cur_v);
        if (q(df_inject)(df, ws, cur_v)) {
          unew ++;
        }
        i_v++;
    }
    assert(usize_old + unew == df->usize);

    if (options & QOP_MDWF_LOG_EIG_POSTAMBLE) {
        for (i = 0; i < df->usize; i++)
            qf(zprint)(s, "postamble", "U %4d %17.9e",
                       i, df->hevals[i]);
    }
    q(df_rebuild)(df);
    return unew;
}
