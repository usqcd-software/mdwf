#include <mdwf.h>

int
QX(dot_fermion)(double *v_r, double *v_i,
		const struct QX(Fermion) *a,
		const struct QX(Fermion) *b)
{
  long long flops;
  DECLARE_STATE;

  CHECK_ARG0(a);
  CHECK_ARGn(b, "dot_fermion");
  CHECK_POINTER(v_r, "dot_fermion");
  CHECK_POINTER(v_i, "dot_fermion");

  BEGIN_TIMING(state);
  flops = qx(dot_fermion)(v_r, v_i,
			  state->even.full_size, state->Ls,
			  a->even, b->even);
  flops += qx(dot_fermion)(v_r, v_i,
			   state->odd.full_size, state->Ls,
			   a->odd, b->odd);
  QMP_sum_double(v_r);
  QMP_sum_double(v_i);
  END_TIMING(state, flops, 2 * sizeof (double), 2 * sizeof (double));
  return 0;
  
}
