#include <mdwf.h>
void
QX(free_half_fermion)(struct QX(HalfFermion) **hf_ptr)
{
  struct Q(State) *state;

  if (hf_ptr == 0 || *hf_ptr == 0)
    return;
  state = (*hf_ptr)->state;
  q(free)(state, *hf_ptr, (*hf_ptr)->size);
  Q(fini)(&state);
}
