#include <mdwf.h>

int
Q(set_Moebius)(struct Q(Parameters) **param_ptr,
               struct Q(State) *state,
               const double b_5[],
               double kappa,
               double M_5,
               double m)
{
  int i;
  double *c;
  int status;

  if (state == NULL || state->error_latched)
    return 1;

  if (param_ptr == NULL)
    return q(set_error)(state, 0, "set_Moebius(): NULL pointer");

  c = q(malloc)(state, state->Ls * sizeof (double));
  if (c == 0)
    return q(set_error)(state, 0, "set_Moebius(): Not enough space");
  for (i = 0; i < state->Ls; i++)
    c[i] = kappa - b_5[i];
  status = Q(set_generic)(param_ptr, state, b_5, c, M_5, m);
  q(free)(state, c, state->Ls * sizeof (double));
  return status;
}
