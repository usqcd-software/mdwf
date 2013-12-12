#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <math.h>
#include <assert.h>
#include <mdwf.h>

/* XXX uses WK: df->df_eigcg.work_c_1 */
int
qx(defl_eigcg_update0)(
        struct QX(Deflator)     *df,
        double                   alpha, 
        double                   beta, 
        double                   alpha_prev, 
        double                   beta_prev, 
        double                   resid_norm_sq, 
        struct Fermion          *resid,
        unsigned int             options)
{
    long long fl = 0;
    double resid_norm;
    struct qx(DeflatorEigcg) *d_e = NULL;
    if (NULL == df || 
            NULL == df->state ||
            ! df->do_eigcg ||
            df->df_eigcg.frozen)
        return 0;       /* skip call */

    d_e = &(df->df_eigcg);   /* shorthand */

    if (d_e->vmax <= d_e->vsize) 
        return -1;      /* should call update1 instead */

    assert(0 < d_e->vsize);
    if (0 == d_e->vsize) 
        return -1;      /* should never get here; 
                           otherwise, the preamble has not been called */

    /* check eig convergence */
    if (resid_norm_sq < d_e->resid_norm_sq_min)
        return 3;   /* eig converged */
    resid_norm = sqrt(resid_norm_sq);
    
    /* [vsize-1, vsize-1] elem */
    if (1 == d_e->vsize) {
        d_e->T[0].r = 1. / alpha;
        d_e->T[0].i = 0.0;
    } else {
        doublecomplex *pT = d_e->T + (d_e->vsize - 1 ) * (1 + d_e->vmax);
        pT->r = 1. / alpha + beta_prev / alpha_prev;
        pT->i = 0.0;
    }   
    assert(d_e->vsize < d_e->vmax);

    /* [vsize-1, vsize], [vsize, vsize-1] elems */
    doublecomplex *pT1 = d_e->T + d_e->vsize - 1 + d_e->vsize * d_e->vmax,
                  *pT2 = d_e->T + d_e->vsize + (d_e->vsize - 1) * d_e->vmax;
    pT1->r = pT2->r = -sqrt(beta) / alpha;
    pT1->i = pT2->i = 0.0;

    /* remember the vector ||resid|| */
    qx(defl_vec) cur_r = d_e->work_c_1;
    qx(defl_vec_copy)(qx(defl_vec_view)(df->state, resid),
                      cur_r);
    fl += qx(defl_vec_scal)(1. / resid_norm, cur_r);
    fl += qx(defl_mat_insert_col)(d_e->V, d_e->vsize, cur_r);

    d_e->vsize += 1;

    df->state->flops += fl;    
    return 0; /* normal */
}
