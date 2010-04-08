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
        struct Q(State)         *s,
        struct Q(Deflator)      *d,
        struct FermionF         *x,
        struct FermionF         *b)
{
    assert(NULL != s &&
            NULL != d &&
            NULL != x &&
            NULL != b);
    latvec_c lv_x   = q(latvec_c_view)(d->dim, d->Ls, x);
    if (d->usize <= 0) {
        q(latvec_c_zero)(lv_x);
        return 0;
    }    
    latvec_c lv_b   = q(latvec_c_view)(d->dim, d->Ls, b);
    latmat_c cur_U  = q(latmat_c_submat_col)(d->U, 0, d->usize);
    q(lat_lmH_dot_lv)(d->usize, cur_U, lv_b, d->zwork);

#if defined(HAVE_LAPACK)
    long int usize  = d->usize,
             ONE    = 1,
             umax   = d->umax,
             info   = 0;
    char cU = 'U';
    zpotrs_(&cU, &usize, &ONE, d->C, &umax, d->zwork, &umax, &info, 1);
    assert(0 == info);
#elif defined(HAVE_GSL)

    gsl_matrix_complex_view gsl_C = gsl_matrix_complex_view_array_with_tda(
            (double *)d->C, d->usize, d->usize, d->umax);
    gsl_vector_complex_view gsl_zwork = gsl_vector_complex_view_array(
            (double *)d->zwork, d->usize);
    CHECK_GSL_STATUS(gsl_linalg_complex_cholesky_svx(
            &gsl_C.matrix, 
            &gsl_zwork.vector));
#else
#  error "no linear algebra library"
#endif

    q(lat_lm_dot_zv)(d->usize, cur_U, d->zwork, lv_x);

    return 0;
}

int
q(df_preamble)(
        struct Q(State)           *s,
        struct Q(Deflator)        *d,
        struct FermionF           *x,
        struct FermionF           *r,
        double                    *r_norm2,
        struct FermionF           *b,  /* const! */
        struct MxM_workspaceF     *ws,
        unsigned int               options)
{
    assert(NULL != s &&
            NULL != x &&
            NULL != b);

    const int dim = DEFLATOR_VEC_SIZE(s);
    int Ls = s->Ls;
    latvec_c lv_x = q(latvec_c_view)(dim, Ls, x);
    latvec_c lv_b = q(latvec_c_view)(dim, Ls, b);
    latvec_c lv_r = q(latvec_c_view)(dim, Ls, r);


    if (d == NULL) {
        q(latvec_c_zero)(lv_x);
        q(latvec_c_copy)(lv_b, lv_r);
        *r_norm2 = q(lat_c_nrm2)(lv_r);
        return 0;
    }

    assert(dim == d->dim);
    Q(deflator_reset)(d);

    if (d->vsize != 0) {
        q(set_error)(s, 0, "df_preamble: deflator in non-initial state");
        return -1;
    }

#define cur_r       (d->work_c_1)
#define cur_Ax      (d->work_c_2)
    if (d->usize <= 0) {
        q(latvec_c_zero)(lv_x);
        q(latvec_c_copy)(lv_b, lv_r);
    } else {
        if (q(df_solve_in_eigenspace)(s, d, x, b))
            return 1;
        /* compute residual */
        latvec_c_linop(cur_Ax, lv_x, ws);
        /* FIXME optimize the code below with a special primitive;
           lv_r <- lv_b - cur_Ax */
        q(latvec_c_copy)(lv_b, lv_r);
        q(lat_c_axpy_d)(-1., cur_Ax, lv_r);
    }

    *r_norm2 = q(lat_c_nrm2)(lv_r);
    if (d->frozen || *r_norm2 <= 0.)
        return 0;

    /* save normalized residual as the first vector */
    q(latvec_c_copy)(lv_r, cur_r);
    q(lat_c_scal_d)(1. / sqrt(*r_norm2), cur_r);
    q(latmat_c_insert_col)(d->V, 0, cur_r);
    
    /* init eigenvalue search: vsize, T, V */
    memset(d->T, 0, d->vmax * d->vmax * sizeof(d->T[0]));
    d->vsize = 1;

    /* init stopping threshold */
    d->resid_norm_sq_min = d->eps * d->eps * (*r_norm2);
    
    return 0;
}
