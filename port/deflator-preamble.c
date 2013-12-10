#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

#include <math.h>
#if defined(HAVE_LAPACK)
#  include <lapack.h>
#elif defined(HAVE_GSL)
#  include <gsl/gsl_linalg.h>
#else
#  error "no linear algebra library"
#endif

int
q(df_solve_in_eigenspace)(
        struct Q(State)         *s, /*remove; use d->state */
        struct Q(Deflator)      *df,
        struct FermionF         *x,
        struct FermionF         *b)
{
    assert(NULL != df &&
            NULL != df->state &&
            NULL != x &&
            NULL != b);
    latvec_c lv_x   = q(latvec_c_view)(df->state, x);
    if (df->usize <= 0) {
        q(latvec_c_zero)(lv_x);
        return 0;
    }    
    latvec_c lv_b   = q(latvec_c_view)(df->state, b);
    latmat_c cur_U  = q(latmat_c_submat_col)(df->U, 0, df->usize);
    q(lat_lmH_dot_lv)(df->usize, cur_U, lv_b, df->zwork);

#if defined(HAVE_LAPACK)
    long int usize  = df->usize,
             ONE    = 1,
             umax   = df->umax,
             info   = 0;
    char cU = 'U';
    zpotrs_(&cU, &usize, &ONE, df->C, &umax, df->zwork, &umax, &info, 1);
    assert(0 == info);
#elif defined(HAVE_GSL)

    gsl_matrix_complex_view gsl_C = gsl_matrix_complex_view_array_with_tda(
            (double *)df->C, df->usize, df->usize, df->umax);
    gsl_vector_complex_view gsl_zwork = gsl_vector_complex_view_array(
            (double *)df->zwork, df->usize);
    CHECK_GSL_STATUS(gsl_linalg_complex_cholesky_svx(
            &gsl_C.matrix, 
            &gsl_zwork.vector));
#else
#  error "no linear algebra library"
#endif

    q(lat_lm_dot_zv)(df->usize, cur_U, df->zwork, lv_x);

    return 0;
}

int
q(df_preamble)(
        struct Q(State)           *s,
        struct Q(Deflator)        *df,
        struct FermionF           *x,
        struct FermionF           *r,
        double                    *r_norm2,
        struct FermionF           *b,  /* const! */
        struct MxM_workspaceF     *ws,
        unsigned int               options)
{
    struct q(DeflatorEigcg) *d_e;
    assert(NULL != df &&
            NULL != df->state &&
            NULL != x &&
            NULL != b);

    latvec_c lv_x = q(latvec_c_view)(s, x);
    latvec_c lv_b = q(latvec_c_view)(s, b);
    latvec_c lv_r = q(latvec_c_view)(s, r);

#define cur_r       (df->work_c_1)
#define cur_Ax      (df->work_c_2)

    if (df->usize <= 0) {
        q(latvec_c_zero)(lv_x);
        q(latvec_c_copy)(lv_b, lv_r);
    } else {
        if (q(df_solve_in_eigenspace)(s, df, x, b))
            return 1;
        /* compute residual */
        latvec_c_linop(cur_Ax, lv_x, ws);
        q(latvec_c_copy)(lv_b, lv_r);
        q(lat_c_axpy_d)(-1., cur_Ax, lv_r);
    }
    
    *r_norm2 = q(lat_c_nrm2)(lv_r);

    if ( ! df->do_eigcg)
        return 0;
    
    d_e = &(df->df_eigcg);

    Q(deflator_reset)(df);
    if (d_e->vsize != 0) {
        q(set_error)(s, 0, "df_preamble: deflator in non-initial state");
        return -1;
    }

    if (d_e->frozen || *r_norm2 <= 0.)
        return 0;

    /* save normalized residual as the first vector */
    q(latvec_c_copy)(lv_r, cur_r);
    q(lat_c_scal_d)(1. / sqrt(*r_norm2), cur_r);
    q(latmat_c_insert_col)(d_e->V, 0, cur_r);
    
    /* init eigenvalue search: vsize, T, V */
    memset(d_e->T, 0, d_e->vmax * d_e->vmax * sizeof(d_e->T[0]));
    d_e->vsize = 1;

    /* init stopping threshold */
    d_e->resid_norm_sq_min = d_e->eps * d_e->eps * (*r_norm2);
    
    return 0;
}
