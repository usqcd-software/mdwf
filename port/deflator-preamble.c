#include <math.h>
#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>



/* XXX uses WK: df->work_c_1, df->work_c_2 (children: df->zwork) */
int
qx(defl_eigcg_preamble)(
        struct QX(Deflator)       *df,
        struct Fermion            *x,
        struct Fermion            *r,
        double                    *r_norm2,
        struct Fermion            *b,  /* const! */
        struct qx(MxM_workspace)  *ws,
        unsigned int               options)
{
    long long fl = 0;
    struct qx(DeflatorEigcg) *d_e = NULL;
    if (NULL == df
            || NULL == df->state
            || NULL == x
            || NULL == b)
        return q(set_error)(df->state, 0, "defl_eigcg_preamble(): null pointer");

    qx(defl_vec) lv_x = qx(defl_vec_view)(df->state, x);
    qx(defl_vec) lv_b = qx(defl_vec_view)(df->state, b);
    qx(defl_vec) lv_r = qx(defl_vec_view)(df->state, r);

    qx(defl_vec) ws_vec = df->work_c_1;

    if (df->usize <= 0) {
        qx(defl_vec_zero)(lv_x);
        qx(defl_vec_copy)(lv_b, lv_r);
    } else {
        if (qx(defl_solve_in_uspace)(df, x, b))
            return 1;
        /* compute residual lv_r <- b - A.x */
        qx(defl_vec_linop)(ws_vec, lv_x, ws);
        qx(defl_vec_copy)(lv_b, lv_r);
        fl += qx(defl_vec_axpy)(-1., ws_vec, lv_r);
    }
    
    fl += qx(defl_vec_nrm2)(r_norm2, lv_r);

    if ( ! df->do_eigcg)
        return 0;
    
    d_e = &(df->df_eigcg);

    QX(deflator_eigcg_reset)(df);
    if (d_e->vsize != 0) {
        q(set_error)(df->state, 0, "defl_eigcg_preamble: deflator in non-initial state");
        return -1;
    }

    if (d_e->frozen || *r_norm2 <= 0.)
        return 0;

    /* save normalized residual as the first vector */
    qx(defl_vec_copy)(lv_r, ws_vec);
    fl += qx(defl_vec_scal)(1. / sqrt(*r_norm2), ws_vec);
    fl += qx(defl_mat_insert_col)(d_e->V, 0, ws_vec);
    
    /* init eigenvalue search: vsize, T, V */
    memset(d_e->T, 0, d_e->vmax * d_e->vmax * sizeof(d_e->T[0]));
    d_e->vsize = 1;

    /* init stopping threshold */
    d_e->resid_norm_sq_min = d_e->eps * d_e->eps * (*r_norm2);
    
    df->state->flops += fl;
    return 0;
}
