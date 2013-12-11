#include <assert.h>
#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

/* FIXME total space occupied by umax will be NCV>NEV from Lanczos;
    1) realloc space 
        - what about alignment?
    2) reuse space for vmax vectors for EigCG 
        - will we realistically use EigCG after Lanczos?
*/

int
Q(create_deflator_inplace)(
        struct Q(Deflator) **deflator_ptr,
        struct QX(HalfFermionMat) **hfm_ptr, 
        int hfm_nev,
        int eigcg_vmax, int eigcg_nev,
        double eigcg_eps, int eigcg_umax)
{
    DECLARE_STATE;
    struct Q(Deflator) *d;
    int dim, Ls;
    int status;

    CHECK_ARG0(*hfm_ptr);

    if (NULL == hfm_ptr || NULL == *hfm_ptr)
        return q(set_error)(state, 0, 
                "create_deflator_inplace_half_fermion_matrix(): NULL matrix");
    if (deflator_ptr == NULL)
        return q(set_error)(state, 0, 
                "create_deflator_inplace_half_fermion_matrix(): NULL pointer");

    d = q(malloc)(state, sizeof (struct Q(Deflator)));
    if (NULL == d)
        return q(set_error)(state, 0, 
                "create_deflator_inplace_half_fermion_matrix(): not enough memory");

    int do_eigcg = (0 < eigcg_vmax) && (0 < eigcg_nev);
    int umax = hfm_nev;
    if (do_eigcg && umax < eigcg_umax)
        umax = eigcg_umax;

    latmat_c *m = &((*hfm_ptr)->m);
    int m_size = m->len;
    if (m_size < umax)
        umax = m_size;
    if (hfm_nev < umax)
        hfm_nev = umax;

    if (0 != (status = q(init_deflator)(d, state, umax, m, hfm_nev, do_eigcg,
                                        eigcg_vmax, eigcg_nev, eigcg_eps)))
        goto clearerr_1;
#if 0
    /* TODO set EigCG to "full" condition */
    if (0 != (status = q(df_recalc_mat)(/*TODO implement and call calc U^H.A.U*/)))
        goto clearerr_1;

    if (0 != (status = q(df_rebuild)(/*TODO pass params*/)))
        goto clearerr_1;
#endif
    return 0;

clearerr_1:
    q(free)(state, d, sizeof (struct Q(Deflator)));
clearerr_0:
    return status;
}
