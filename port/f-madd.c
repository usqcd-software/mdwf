#include <mdwf.h>

int
QX(madd_fermion)(struct QX(Fermion) *r,
		 const struct QX(Fermion) *a,
		 double alpha,
		 const struct QX(Fermion) *b)
{
  long long flops;
  DECLARE_STATE;

  CHECK_ARG0(r);
  CHECK_ARGn(a, "madd_fermion");
  CHECK_ARGn(b, "madd_fermion");

  BEGIN_TIMING(state);
  flops = qx(f_add3)(r->even,
		     state->even.full_size, state->Ls,
		     a->even, alpha, b->even);
  flops += qx(f_add3)(r->odd,
		      state->odd.full_size, state->Ls,
		      a->odd, alpha, b->odd);
  END_TIMING(state, flops, 0, 0);
  return 0;
  
}
