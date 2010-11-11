#include <mdwf.h>

int
QX(export_4d_fermion)(void (*writer)(const int pos[Q(DIM)],
                                     int color,
                                     int dirac,
                                     double val_re,
                                     double val_im,
                                     void *env),
                      void *env,
                      const struct QX(Fermion) *fermion)
{
  if (fermion == 0)
    return 1;

  struct Q(State) *state = fermion->state;
  int size = state->Ls * Q(FERMION_DIM) * Q(COLORS) * 2 * sizeof (double);
  double *m = q(malloc)(state, size);

  if (m == 0)
	  return 1;

  BEGIN_TIMING(state);
  qx(x4_export)(&state->even, m, fermion->even, writer, env);
  qx(x4_export)(&state->odd, m, fermion->odd, writer, env);
  END_TIMING(state, 0, 0, 0);
  q(free)(state, m, size);

  return 0;
}
