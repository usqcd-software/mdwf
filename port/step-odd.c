#include <mdwf.h>

void *
qx(step_odd)(struct Q(State) *state, void *aligned_ptr)
{
    int size = 2 * Q(COLORS) * Q(FERMION_DIM) * sizeof (REAL);

    if (state == 0 || aligned_ptr == 0)
        return NULL;
    return ALIGN(aligned_ptr, state->Ls * state->odd.full_size * size);
}
