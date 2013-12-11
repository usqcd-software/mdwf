#include <assert.h>
#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

/* FIXME total space occupied by umax will be NCV>NEV from Lanczos;
    1) realloc space 
        - what about alignment?
    2) reuse space for vmax vectors for EigCG 
        - will we realistically use EigCG after Lanczos?
*/
/* move vectors in latmat from hfm_ptr to deflator; 
   hfm_ptr can be safely free'd after that */

int
QX(create_deflator_inplace)(
        struct Q(Deflator) **deflator_ptr,
        const struct Q(Parameters)  *params,
        const struct QX(Gauge)      *gauge,
        struct QX(HalfFermionMat) *hfm_ptr, 
        int hfm_nev,    /* take `hfm_nev' first vectors from the matrix */
        int eigcg_vmax, 
        int eigcg_nev,
        double eigcg_eps,
        int eigcg_umax)
{
    DECLARE_STATE;
    struct Q(Deflator) *df;
    int status;

    CHECK_ARG0(hfm_ptr);
    CHECK_ARGn(params, "create_deflator_inplace");
    CHECK_ARGn(gauge, "create_deflator_inplace");


    if (NULL == hfm_ptr || NULL == hfm_ptr)
        return q(set_error)(state, 0, 
                "create_deflator_inplace(): NULL matrix");
    if (deflator_ptr == NULL)
        return q(set_error)(state, 0, 
                "create_deflator_inplace(): NULL pointer");

    df = q(malloc)(state, sizeof (struct Q(Deflator)));
    if (NULL == df)
        return q(set_error)(state, 0, 
                "create_deflator_inplace(): not enough memory");
    
    BEGIN_TIMING(state);

    /* init operator workspace */
    struct MxM_workspace  ws;
    void *temps = NULL, 
         *temps_ptr = NULL;
    size_t temps_size = 0;
    long long flops = 0, sent = 0, received = 0;
    temps_ptr = qx(allocate_eo)(state, &temps_size, &temps, 0, 2, 1);
    if (NULL == temps_ptr) {
        status = q(set_error)(state, 0, 
                "create_deflator_inplace(): not enough memory");
        goto clearerr_1;
    }

    ws.state     = state;
    ws.params    = params;
    ws.gauge     = gauge->data;
    ws.tmp_e     = temps;
    ws.tmp2_e    = temps = qx(step_even)(state, temps);
    ws.tmp_o     = temps = qx(step_even)(state, temps);
    ws.flops     = &flops;
    ws.sent      = &sent;
    ws.received  = &received;


    /* create deflator with existing vectors, leave door open to EigCG */
    latmat_c *m = &(hfm_ptr->m);
    int umax    = m->len; /* `m' already allocated ; cannot extend */
    int do_eigcg = (0 < eigcg_vmax) && (0 < eigcg_nev);
    if (do_eigcg && eigcg_umax < umax)
        umax = eigcg_umax;  /* limit umax to fill by EigCG */
    if (hfm_nev < umax)
        hfm_nev = umax;

    if (0 != (status = q(init_deflator)(df, state, umax, m, hfm_nev, do_eigcg,
                                        eigcg_vmax, eigcg_nev, eigcg_eps)))
        goto clearerr_2;

    /* initialize deflator internal matrices */
    if (0 != (status = q(df_recalc_mat)(df, &ws)))
        goto clearerr_2;

    if (0 != (status = q(df_rebuild)(df)))
        goto clearerr_2;
    
    df->df_eigcg.frozen = 1;    /* set EigCG to "frozen" condition */
    df->loading = 0;            /* set deflator to "ready" */
    *deflator_ptr = df;

    END_TIMING(state, flops, sent, received);
    
    q(free)(state, temps_ptr, temps_size);
    return 0;

clearerr_2:
    q(free)(state, temps_ptr, temps_size);
clearerr_1:
    q(free)(state, df, sizeof (struct Q(Deflator)));
clearerr_0:
    return status;
}
