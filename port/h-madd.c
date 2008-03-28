#include <mdwf.h>

int
QX(madd_half_fermion)(struct QX(HalfFermion) *r,
		      const struct QX(HalfFermion) *a,
		      double alpha,
		      const struct QX(HalfFermion) *b)
{
  long long flops;
  DECLARE_STATE;

  CHECK_ARG0(r);
  CHECK_ARGn(a, "madd_half_fermion");
  CHECK_ARGn(b, "madd_half_fermion");

  BEGIN_TIMING(state);
  flops = qx(f_add3)(r->even,
		     state->even.full_size, state->Ls,
		     a->even, alpha, b->even);
  END_TIMING(state, flops, 0, 0);
  return 0;
  
}
