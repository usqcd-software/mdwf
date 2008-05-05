#include <mdwf.h>

int
QX(import_fermion)(struct QX(Fermion) **fermion_ptr,
		   struct Q(State) *state,
		   double (*reader)(const int pos[5],
				    int color,
				    int dirac,
				    int re_im,
				    void *env),
		   void *env)
{
  double *m;
  int size;

  if (QX(allocate_fermion)(fermion_ptr, state) != 0)
    return 1;

  size = state->Ls * Q(FERMION_DIM) * Q(COLORS) * 2 * sizeof (double);
  m = q(malloc)(state, size);
  if (m == 0) {
    QX(free_fermion)(fermion_ptr);
    return q(set_error)(state, 0, "import_fermion(): not enough space");
  }
  BEGIN_TIMING(state);
  qx(x_import)(&state->even, m, (*fermion_ptr)->even, reader, env);
  qx(x_import)(&state->odd, m, (*fermion_ptr)->odd, reader, env);
  END_TIMING(state, 0, 0, 0);

  q(free)(state, m, size);
  return 0;
}
