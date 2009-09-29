#include <mdwf.h>

int
QX(allocate_half_fermion)(struct QX(HalfFermion) **hf_ptr,
                          struct Q(State) *state)
{
  struct QX(HalfFermion) *hf;
  void *even;
  size_t size;

  if (state == NULL || state->error_latched)
    return 1;

  if (hf_ptr == NULL)
    return q(set_error)(state, 0, "allocate_half_fermion(): NULL pointer");

  *hf_ptr = NULL;
  hf = qx(allocate_eo)(state, &size, &even,
                       sizeof (struct QX(Fermion)), 1, 0);
  if (hf == 0)
    return q(set_error)(state, 0, "allocate_half_fermion(): not enough memory");

  BEGIN_TIMING(state);
  hf->state = state;
  hf->size = size;
  hf->even = even;
  *hf_ptr = hf;
  END_TIMING(state, 0, 0, 0);

  return 0;
}
