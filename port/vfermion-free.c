#include <mdwf.h>
void
QX(free_vector_fermion)(struct QX(VectorFermion) **vfermion_ptr)
{
    struct Q(State) *state;

    if (vfermion_ptr == 0 || *vfermion_ptr == 0)
	return;
    state = (*vfermion_ptr)->state;
    BEGIN_TIMING(state);
    q(free)(state, *vfermion_ptr, (*vfermion_ptr)->size);
    END_TIMING(state, 0, 0, 0);
    *vfermion_ptr = 0;
}
