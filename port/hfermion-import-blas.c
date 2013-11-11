#include <mdwf.h>

int
QX(half_fermion_from_blas)(struct QX(HalfFermion) *half_fermion,
                           const REAL *data,
                           int data_size)
{
  struct Q(State) *state;
  int size;

  if (half_fermion == 0)
    return 1;

  if (data == 0)
    return 1;

  state = half_fermion->state;
  size = state->Ls * state->even.full_size * Q(FERMION_DIM) * Q(COLORS) * 2;
  if (data_size < size)
    return q(set_error(state, 0, "half_fermion_from_blas(): not enough space"));

  BEGIN_TIMING(state);
  qx(blas2fermion)(half_fermion->even, state->even.full_size, state->Ls, data);
  END_TIMING(state, 0, 0, 0);
  return 0;
}
