#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>
#include <math.h>

int
q(df_update0)(
        struct Q(State)         *s,
        struct Q(Deflator)      *d,
        double                   alpha, 
        double                   beta, 
        double                   alpha_prev, 
        double                   beta_prev, 
        double                   resid_norm_sq, 
        struct FermionF         *resid,
        unsigned int             options)
{
    if (NULL == d || 
            NULL == s ||
            d->frozen)
        return 0;       /* skip call */

    if (d->vmax <= d->vsize) 
        return -1;      /* should call update1 instead */

    assert(0 < d->vsize);
    if (0 == d->vsize) 
        return -1;      /* should never get here; 
                           otherwise, the preamble has not been called */

    /* check eig convergence */
    if (resid_norm_sq < d->resid_norm_sq_min)
        return 3;   /* eig converged */
    double resid_norm = sqrt(resid_norm_sq);
    
    /* [vsize-1, vsize-1] elem */
    if (1 == d->vsize) {
        d->T[0].r = 1. / alpha;
        d->T[0].i = 0.0;
    } else {
        doublecomplex *pT = d->T + (d->vsize - 1 ) * (1 + d->vmax);
        pT->r = 1. / alpha + beta_prev / alpha_prev;
        pT->i = 0.0;
    }   
    assert(d->vsize < d->vmax);

    /* [vsize-1, vsize], [vsize, vsize-1] elems */
    doublecomplex *pT1 = d->T + d->vsize - 1 + d->vsize * d->vmax,
                  *pT2 = d->T + d->vsize + (d->vsize - 1) * d->vmax;
    pT1->r = pT2->r = -sqrt(beta) / alpha;
    pT1->i = pT2->i = 0.0;

    /* remember the vector ||resid|| */
#define cur_r   (d->work_c_1)
    q(latvec_c_copy)(q(latvec_c_view)(d->dim, d->Ls, resid), 
                  cur_r);
    q(lat_c_scal_d)(1. / resid_norm, cur_r);
    q(latmat_c_insert_col)(d->V, d->vsize, cur_r);
#undef cur_r

    d->vsize += 1;
    
    return 0; /* normal */
}
