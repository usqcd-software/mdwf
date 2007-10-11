#include <mdwf.h>

int
QX(import_gauge)(struct QX(Gauge) **gauge_ptr,
		 struct Q(State) *state,
		 double (*reader)(int dir,
				  const int pos[Q(DIM)],
				  int a,
				  int b,
				  int re_im,
				  void *env),
		 void *env)
{
  struct QX(Gauge) *gauge;
  void *ptr;
  size_t size;
  double vr, vi;
  int p, q, i, d, a, b;
  int x[Q(DIM)], dx[Q(DIM)];

  if (state == 0 || state->error_latched)
    return 1;

  *gauge_ptr = NULL;
  gauge = q(allocate_aligned)(state, &size, &ptr,
			      sizeof (struct QX(Gauge)),
			      Q(DIM) * state->volume * sizeof (struct SUn));
  if (gauge == 0)
    return q(set_error)(state, 0, "import_gauge(): not enough memory");
  state->used++;
  gauge->state = state;
  gauge->size = size;
  gauge->data = ptr;
  *gauge_ptr = gauge;

  for (i = 0; i < Q(DIM); i++)
    dx[i] = state->sublattice.hi[i] - state->sublattice.lo[i];

  for (p = 0; p < state->volume; p++) {
    for (q = state->layout2linear[p], i = 0; i < Q(DIM); i++) {
      x[i] = q % dx[i] + state->sublattice.lo[i];
      q = q / dx[i];
    }
    for (d = 0; d < Q(DIM); d++) {
      for (a = 0; a < Q(COLORS); a++) {
	for (b = 0; b < Q(COLORS); b++) {
	  vr = - 0.5 * reader(d, x, a, b, 0, env);
	  vi = - 0.5 * reader(d, x, a, b, 1, env);
	  qx(set_gauge)(gauge->data, Q(DIM) * p + d, a, b, vr, vi);
	}
      }
    }
  }
  return 0;
}
