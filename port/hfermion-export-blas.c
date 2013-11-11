#include <mdwf.h>

int
QX(blas_from_half_fermion)(REAL *data,
                           int data_size,
                           const struct QX(HalfFermion) *half_fermion)
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
    return q(set_error(state, 0, "blass_from_half_fermion(): not enough space"));

  BEGIN_TIMING(state);
  qx(fermion2blas)(data, half_fermion->even, state->even.full_size, state->Ls);
  END_TIMING(state, 0, 0, 0);
  return 0;
}
