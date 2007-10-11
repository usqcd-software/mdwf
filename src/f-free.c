#include <mdwf.h>
void
QX(free_fermion)(struct QX(Fermion) **fermion_ptr)
{
  struct Q(State) *state;

  if (fermion_ptr == 0 || *fermion_ptr == 0)
    return;
  state = (*fermion_ptr)->state;
  q(free)(state, *fermion_ptr, (*fermion_ptr)->size);
  Q(fini)(&state);
}
