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

/* XXX uses WK: df->work_c_2, df->zwork */
int 
q(df_ortho_uspace)(latvec_c vec,
                   struct Q(Deflator) *df,
                   int u_lo, int u_hi)
{
    if (NULL == df || latvec_c_is_null(&vec))
        return q(set_error)(df->state, 0, "df_ortho_uspace(): null pointer(s)");
    
    int fl = 0,
        u_len = u_hi - u_lo;
    if (u_len <= 0)
        return 0; /* relax */

    latmat_c cur_U = q(latmat_c_submat_col)(df->U, u_lo, u_len);
    latvec_c ws_vec = df->work_c_2;

    /*fl += */q(lat_lmH_dot_lv)(u_len, cur_U, 
                            vec, df->zwork);
    /*fl += */q(lat_lm_dot_zv)(u_len, cur_U,
                           df->zwork, ws_vec);
    /*fl += */q(lat_c_axpy_d)(-1., ws_vec, vec);
    df->state->flops += fl;

    return 0;
}

/* TODO add flop counting */
/* TODO change latvec_c_copy(x, y) to (y, x) because it does y<-x */

/* orthogonalize wrt [u_lo: u_hi) vectors, normalize, inject into [u_pos]
   and compute H[u_lo:u_hi, vpos]
   XXX uses WK: (children: df->work_c_2, df->zwork)
       vector `vec' is orthogonal to U and |vec|=1 on return */
int
q(df_inject)(struct Q(Deflator) *df,
             struct MxM_workspaceF *ws,
             int u_lo, int u_hi, int u_pos,
             latvec_c vec)
{
    double v_norm2_min = eps_reortho * eps_reortho;
    int i, //fl = 0, 
        u_len = u_hi - u_lo;

    /* reorthogonalize n_reortho times */
    for (i = n_reortho; i--; )
        q(df_ortho_uspace)(vec, df, u_lo, u_hi);

    double v_norm2 = q(lat_c_nrm2)(vec);
    if (v_norm2 < v_norm2_min)
        return 0;

    /* normalize & insert the vector */
    q(lat_c_scal_d)(1. / sqrt(v_norm2), vec);
    q(latmat_c_insert_col)(df->U, u_pos, vec);

    latvec_c ws_vec = df->work_c_2;
    latvec_c_linop(ws_vec, vec, ws);

    if (0 < u_len) {
        /* compute and store (U^H . A . vec) */
        latmat_c cur_U = q(latmat_c_submat_col)(df->U, u_lo, u_len);
        q(lat_lmH_dot_lv)(u_len, cur_U, ws_vec, df->H + u_lo + u_pos * df->umax);
        for (i = u_lo ; i < u_hi ; i++) {
            doublecomplex *p1 = df->H + i + u_pos * df->umax,
                          *p2 = df->H + u_pos + i * df->umax;
            p2->r =  p1->r;
            p2->i = -p1->i;
        }
    }
    /* compute and store (vec^H . A . vec) */
    doublecomplex vAv;
    vAv = q(lat_c_dotu)(vec, ws_vec);
    df->H[u_pos * (df->umax + 1)].r = vAv.r;
    df->H[u_pos * (df->umax + 1)].i = 0.;

    /* finished adding one new vector */
    return 1;
}

/* orthogonalize, normalize, add 'vec' to U
   XXX uses WK: (children: df->work_c_2, df->zwork)
       vector `vec' is orthogonal to U and |vec|=1 on return */
int 
q(df_inject_back)(struct Q(Deflator) *df,
                  struct MxM_workspaceF *ws,
                  latvec_c vec)
{
    if (df->usize < df->umax) {
        q(df_inject)(df, ws, 0, df->usize, df->usize, vec);
        df->usize ++;
        return 1;
    } else
        return 0;
}
/* TODO replace all usages of df_inject -> df_inject_back */

/* use df_inject below : extract vector; in ject vector */
/* recompute H[:usize, :usize]; all vectors are reorthogonalized 
   with respect to preceding vectors in U
   XXX uses WK: df->work_c_1 (children: df->work_c_2, df->zwork) */
int 
q(df_recalc_mat) (struct Q(Deflator) *df,
                  struct MxM_workspaceF *ws)
{
    double v_norm2_min = eps_reortho * eps_reortho;
    int i, u_pos;

    latvec_c vec = df->work_c_1;

    for (i = 0 ; i < df->usize ; i++) {
        q(latmat_c_get_col)(df->U, i, vec);
        /*increment only if have an indep. vector */
        u_pos += q(df_inject)(df, ws, 0, u_pos, u_pos, vec);
    }
  
    return u_pos;
}

/* recompute Cholesky factorization for U-space solving */
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

    /* broadcast Cholesky matrix from the master node for consistency */
    for (i = 0 ; i < df->usize; i++)
        QMP_broadcast((void *)(df->C + i * df->umax), 
                      df->usize * sizeof(df->C[0]));

    if (0 != (status = q(deflator_calc_evals)(df)))
        return status;

    return 0;
}

/* XXX uses WK: df->work_c_1 (children: df->work_c_2, df->zwork) */
int
q(df_postamble)(
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
        if (q(df_inject_back)(df, ws, cur_v)) {
          unew ++;
        }
        i_v++;
    }
    assert(usize_old + unew == df->usize);

    if (options & QOP_MDWF_LOG_EIG_POSTAMBLE) {
        for (i = 0; i < df->usize; i++)
            qf(zprint)(df->state, "postamble", "U %4d %17.9e",
                       i, df->hevals[i]);
    }
    q(df_rebuild)(df);
    return unew;
}
