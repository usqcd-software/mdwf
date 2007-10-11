#include <mdwf.h>

int
QX(allocate_half_fermion)(struct QX(HalfFermion) **hf_ptr,
   		          struct Q(State) *state)
{
  struct QX(HalfFermion) *hf;
  void *even;
  size_t size;

  if (state == 0 || state->error_latched)
    return 1;

  *hf_ptr = NULL;
  hf = q(allocate_eo)(state, &size, &even,
		      sizeof (struct QX(Fermion)), 1, 0,
		      sizeof (REAL));
  if (hf == 0)
    return q(set_error)(state, 0, "allocate_half_fermion(): not enough memory");

  state->used++;
  hf->state = state;
  hf->size = size;
  hf->even = even;
  *hf_ptr = hf;

  return 0;
}
