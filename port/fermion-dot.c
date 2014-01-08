#include <mdwf.h>

int
QX(dot_fermion)(double *v_r, double *v_i,
                const struct QX(Fermion) *a,
                const struct QX(Fermion) *b)
{
  long long flops;
  double e_r = 0;
  double e_i = 0;
  double o_r = 0;
  double o_i = 0;
  DECLARE_STATE;

  CHECK_ARG0(a);
  CHECK_ARGn(b, "dot_fermion");
  CHECK_POINTER(v_r, "dot_fermion");
  CHECK_POINTER(v_i, "dot_fermion");

  BEGIN_TIMING(state);
  flops = qx(f_dot)(&e_r, &e_i,
                    state->even.full_size, state->Ls,
                    a->even, b->even);
  flops += qx(f_dot)(&o_r, &o_i,
                     state->odd.full_size, state->Ls,
                     a->odd, b->odd);
  *v_r = e_r + o_r;
  *v_i = e_i + o_i;
  QMP_sum_double(v_r);
  QMP_sum_double(v_i);
  END_TIMING(state, flops, 2 * sizeof (double), 2 * sizeof (double));
  return 0;
  
}
