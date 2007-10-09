#include <mdwf.h>

void *
q(step_odd)(struct Q(State) *state, void *aligned_ptr, size_t fsize)
{
  if (state == 0 || aligned_ptr == 0)
    return NULL;
  return ALIGN(aligned_ptr + state->Ls * state->odd.full_size * fsize);
}