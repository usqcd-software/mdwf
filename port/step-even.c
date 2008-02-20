#include <mdwf.h>

void *
q(step_even)(struct Q(State) *state, void *aligned_ptr, size_t fsize)
{
  int size = fsize * 2 * Q(COLORS) * Q(FERMION_DIM);

  if (state == 0 || aligned_ptr == 0)
    return NULL;
  return ALIGN(aligned_ptr + state->Ls * state->even.full_size * size);
}
