#include <mdwf.h>

int
QX(M_operator)(struct QX(HalfFermion) *result,
	       const struct Q(Parameters) *params,
	       const struct QX(Gauge) *gauge,
	       const struct QX(HalfFermion) *fermion)
{
    DECLARE_STATE;
    long long flops = 0;
    long long sent = 0;
    long long received = 0;
    size_t alloc_size;
    void *alloc_ptr;
    void *aligned_ptr;
    struct Fermion *odd_b;
    
    CHECK_ARG0(result);
    CHECK_ARGn(params, "M_operator");
    CHECK_ARGn(gauge, "M_operator");
    CHECK_ARGn(fermion, "M_operator");
    
    if (q(setup_comm)(state, sizeof (REAL)))
	return q(set_error)(state, 0,
			    "M_operator(): communication setup failed");

    alloc_ptr = q(allocate_eo)(state, &alloc_size, &aligned_ptr,
			       0, 0, 1, sizeof (REAL));
    if (alloc_ptr == 0)
	return q(set_error)(state, 0, "DDW_operator(): not enough memory");
    odd_b = aligned_ptr;
    
    BEGIN_TIMING(state);
    qx(op_M)(result->even, &state->even, params,
	     gauge->data, fermion->even,
	     &flops, &sent, &received,
	     odd_b);
    END_TIMING(state, flops, sent, received);

    q(free)(state, alloc_ptr, alloc_size);

    return 0;
}	
