#include <mdwf.h>

int
QX(allocate_fermion)(struct QX(Fermion) **fermion_ptr,
		     struct Q(State) *state)
{
  struct QX(Fermion) *fermion;
  void *even;
  size_t size;

  if (state == 0 || state->error_latched)
    return 1;

  *fermion_ptr = NULL;
  fermion = q(allocate_eo)(state, &size, &even,
			   sizeof (struct QX(Fermion)), 1, 1,
			   sizeof (REAL));
  if (fermion == 0)
    return q(set_error)(state, 0, "allocate_fermion(): not enough memory");

  state->used++;
  fermion->state = state;
  fermion->size = size;
  fermion->even = even;
  fermion->odd = q(step_even)(state, even, sizeof (REAL));
  *fermion_ptr = fermion;

  return 0;
}
