#include <mdwf.h>

int
QX(export_fermion)(void (*writer)(const int pos[5],
				  int color,
				  int dirac,
				  int re_im,
				  double value,
				  void *env),
		   void *env,
		   const struct QX(Fermion) *fermion)
{
  struct Q(State) *state;
  double *m;
  int size;

  if (fermion == 0)
    return 1;

  state = fermion->state;
  size = state->Ls * Q(FERMION_DIM) * Q(COLORS) * 2 * sizeof (double);
  m = q(malloc)(state, size);
  if (m == 0) {
    return q(set_error)(state, 0, "export_fermion(): not enough space");
  }
  BEGIN_TIMING(state);
  q(x_export)(&state->even, m, fermion->even, writer, env);
  q(x_export)(&state->odd, m, fermion->odd, writer, env);
  END_TIMING(state, 0, 0, 0);

  q(free)(state, m, size);
  return 0;
}
