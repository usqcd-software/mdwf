#include <mdwf.h>

/* NOTE: not the beset overlap of comm/comm */

int
QX(DDW_operator)(struct QX(Fermion) *result,
		 const struct Q(Parameters) *params,
		 const struct QX(Gauge) *gauge,
		 const struct QX(Fermion) *source)
{
  DECLARE_STATE;
  long long flops = 0;
  long long sent = 0;
  long long received = 0;
  size_t alloc_size;
  void *alloc_ptr;
  void *aligned_ptr;
  struct Fermion *even_b;
  struct Fermion *odd_b;

  CHECK_ARG0(result);
  CHECK_ARGn(params, "DDW_operator");
  CHECK_ARGn(gauge, "DDW_operator");
  CHECK_ARGn(source, "DDW_operator");

  if (q(setup_comm)(state, sizeof (REAL)))
      return q(set_error)(state, 0, "DDW_operator(): communication setup failed");

  alloc_ptr = q(allocate_eo)(state, &alloc_size, &aligned_ptr,
			     0, 1, 1, sizeof (REAL));
  if (alloc_ptr == 0)
      return q(set_error)(state, 0, "DDW_operator(): not enough memory");

  even_b = aligned_ptr;
  odd_b = q(step_even)(state, even_b, sizeof (REAL));

  BEGIN_TIMING(state);
  qx(op_B)(even_b, &state->even, params, source->even, &flops);
  qx(op_B)(odd_b, &state->odd, params, source->odd, &flops);
  qx(op_ApF)(result->even, &state->even, params,
	     gauge->data, source->even, odd_b,
	     &flops, &sent, &received);
  qx(op_ApF)(result->odd, &state->odd, params,
	     gauge->data, source->odd, even_b,
	     &flops, &sent, &received);
  END_TIMING(state, flops, sent, received);

  q(free)(state, alloc_ptr, alloc_size);

  return 0;
}	
