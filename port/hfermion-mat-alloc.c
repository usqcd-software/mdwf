#include <assert.h>
#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

int
QX(allocate_half_fermion_matrix)(
        struct QX(HalfFermionMat) **hfm_ptr,
        REAL **blas_ptr,
        size_t *blas_ld,
        struct Q(State) *state,
        int ncol)
{
    latmat_c *m;
    int Ls, dim;
    struct QX(HalfFermionMat) *hfm;
    struct vFermion *fv;
    
    if (NULL == hfm_ptr)
        return q(set_error)(state, 0, "allocate_half_fermion_matrix(): NULL pointer");
    *hfm_ptr = NULL;

    if (NULL == state || state->error_latched)
        return 1;
    Ls  = state->Ls;
    dim = state->even.full_size;

    hfm = q(malloc)(state, sizeof(struct QX(HalfFermionMat)));
    if (NULL == hfm)
        return q(set_error)(state, 0, "allocate_half_fermion_matrix(): not enough memory");

    hfm->mem_size = sizeof(struct QX(HalfFermionMat));
    hfm->m = q(latmat_c_alloc)(state, ncol);
    if (latmat_c_is_null(&(hfm->m))) {
        q(free)(state, hfm, hfm->mem_size);
        return q(set_error)(state, 0, "allocate_half_fermion_matrix(): not enough memory");
    }

    *hfm_ptr = hfm;
    if (NULL != blas_ptr)
        *blas_ptr = (REAL *)m->fv;
    if (NULL != blas_ld) {
        assert(0 == m->stride % (2 * sizeof(REAL)));
        *blas_ld = m->stride / (2 * sizeof(REAL));
    }

    return 0;
}
