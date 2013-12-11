#include <assert.h>
#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

int
QX(alloc_half_fermion_matrix)(
        struct QX(HalfFermionMat) **hfm_ptr,
        struct Q(State) *state,
        int ncol)
{
    struct QX(HalfFermionMat) *hfm;
    struct vFermion *fv;
    
    if (NULL == hfm_ptr)
        return q(set_error)(state, 0, "allocate_half_fermion_matrix(): NULL pointer");
    *hfm_ptr = NULL;

    if (NULL == state || state->error_latched)
        return 1;

    hfm = q(malloc)(state, sizeof(struct QX(HalfFermionMat)));
    if (NULL == hfm)
        return q(set_error)(state, 0, "allocate_half_fermion_matrix(): not enough memory");

    hfm->m = qx(defl_mat_alloc)(state, ncol);
    if (defl_mat_is_null(&(hfm->m))) {
        q(free)(state, hfm, hfm->mem_size);
        return q(set_error)(state, 0, "allocate_half_fermion_matrix(): not enough memory");
    }

    hfm->state      = state;
    hfm->mem_size   = sizeof(struct QX(HalfFermionMat));

    *hfm_ptr = hfm;

    return 0;
}
int QX(blas_view_half_fermion_matrix)(
        struct QX(HalfFermionMat) *hfm_ptr,
        int *nrow_loc,
        int *ncol,
        REAL **blas_ptr,
        size_t *blas_ld)
{
    if (NULL == hfm_ptr)
        return q(set_error)(hfm_ptr->state, 0, "blas_view_half_fermion_matrix(): null pointer");

    if (NULL != nrow_loc) 
        *nrow_loc = hfm_ptr->state->Ls * hfm_ptr->state->even.full_size;
    if (NULL != ncol) 
        *ncol = hfm_ptr->m.len;
    if (NULL != blas_ptr)
        *blas_ptr = (REAL *)(hfm_ptr->m.fv);
    if (NULL != blas_ld) {
        assert(0 == hfm_ptr->m.stride % (2 * sizeof(REAL)));
        *blas_ld = hfm_ptr->m.stride / (2 * sizeof(REAL));
    }
    return 0;
}
