#include <mdwf.h>

int
Q(set_Chiu)(struct Q(Parameters) **param_ptr,
            struct Q(State) *state,
            const double *a_5[],
            double M_5,
            double m)
{
  int i;
  int status;

  if (state == 0 || state->error_latched)
    return 1;

  status = Q(set_generic)(param_ptr, state, a_5, a_5, M_5, m);
  return status;
}
