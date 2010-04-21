#include <mdwf.h>

int
QX(midpoint_pseudo)(void (*writer)(const int pos[Q(DIM)],
                                   double value,
                                   void *env),
                    void *env,
                    const struct QX(Fermion) *fermion)
{
  if (fermion == 0)
    return 1;

  struct Q(State) *state = fermion->state;
  int size = state->Ls * Q(FERMION_DIM) * Q(COLORS) * 2;
  double m[size];

  BEGIN_TIMING(state);
  qx(x_midpoint)(&state->even, m, fermion->even, writer, env);
  qx(x_midpoint)(&state->odd, m, fermion->odd, writer, env);
  END_TIMING(state, 0, 0, 0);

  return 0;
}
