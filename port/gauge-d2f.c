#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>
#undef QOP_MDWF_DEFAULT_PRECISION
#define QOP_MDWF_DEFAULT_PRECISION 'D'
#include <mdwf.h>

int
Q(gauge_float_from_double)(struct QF(Gauge) **result_ptr,
                           struct QD(Gauge) *gauge_d)
{
  DECLARE_STATE;
  struct QF(Gauge) *g;
  int us;
  void *ptr;
  size_t size;

  CHECK_ARG0(gauge_d);
  *result_ptr = NULL;
  us = qx(sizeof_gauge)(state->volume);
  g = q(allocate_aligned)(state, &size, &ptr,
                          sizeof (struct QF(Gauge)), us);
  if (g == 0)
    return q(set_error)(state, 0, "gauge_float_from_double(): not enough memory");
  g->state = state;
  g->size = size;
  g->data = ptr;
  *result_ptr = g;
  
  q(g_f_eq_d)(g->data, Q(DIM) * state->volume, gauge_d->data);

  return 0;
}
