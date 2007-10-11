#include <mdwf.h>
void
QX(free_gauge)(struct QX(Gauge) **g_ptr)
{
  struct Q(State) *state;

  if (g_ptr == 0 || *g_ptr == 0)
    return;
  state = (*g_ptr)->state;
  q(free)(state, *g_ptr, (*g_ptr)->size);
  Q(fini)(&state);
}
