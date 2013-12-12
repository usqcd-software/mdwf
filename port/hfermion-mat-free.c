#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

void
QX(free_half_fermion_matrix)(struct QX(HalfFermionMat) **hfm_ptr)
{
    struct Q(State) *state;

    if (hfm_ptr == 0 || *hfm_ptr == 0)
        return;
    state = (*hfm_ptr)->state;

    qx(defl_mat_free)(state, &((*hfm_ptr)->m));
    q(free)(state, *hfm_ptr, (*hfm_ptr)->mem_size);
    *hfm_ptr = 0;
}

