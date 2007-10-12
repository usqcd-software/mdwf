#include <mdwf.h>
void
QX(free_half_fermion)(struct QX(HalfFermion) **hf_ptr)
{
  struct Q(State) *state;

  if (hf_ptr == 0 || *hf_ptr == 0)
    return;
  state = (*hf_ptr)->state;
  BEGIN_TIMING(state);
  q(free)(state, *hf_ptr, (*hf_ptr)->size);
  END_TIMING(state, 0, 0, 0);
  Q(fini)(&state);
}
