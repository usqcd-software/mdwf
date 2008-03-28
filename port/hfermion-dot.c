#include <mdwf.h>

int
QX(dot_half_fermion)(double *v_r, double *v_i,
		     const struct QX(HalfFermion) *a,
		     const struct QX(HalfFermion) *b)
{
  long long flops;
  DECLARE_STATE;

  CHECK_ARG0(a);
  CHECK_ARGn(b, "dot_half_fermion");
  CHECK_POINTER(v_r, "dot_half_fermion");
  CHECK_POINTER(v_i, "dot_half_fermion");

  BEGIN_TIMING(state);
  flops = qx(f_dot)(v_r, v_i,
		    state->even.full_size, state->Ls,
		    a->even, b->even);
  QMP_sum_double(v_r);
  QMP_sum_double(v_i);
  END_TIMING(state, flops, 2 * sizeof (double), 2 * sizeof (double));
  return 0;
  
}
