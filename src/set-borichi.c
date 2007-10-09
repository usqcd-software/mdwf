#include <mdwf.h>

int
Q(set_Borichi)(struct Q(Parameters) **param_ptr,
               struct Q(State) *state,
               double a_5,
               double M_5,
               double m)
{
  int i;
  double *b;
  int status;

  if (state == 0 || state->error_latched)
    return 1;

  b = q(malloc)(state, state->Ls * sizeof (double));
  if (b == 0)
    return q(set_error)(state, 0, "set_Borichi(): Not enough space");
  for (i = 0; i < state->Ls; i++) {
    b[i] = a_5;
  }
  status = Q(set_generic)(param_ptr, state, b, b, M_5, m);
  q(free)(state, b, state->Ls * sizeof (double));
  return status;
}
