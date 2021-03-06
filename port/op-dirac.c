#include <mdwf.h>

/* NOTE: not the best overlap of comm/comm */

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
      return q(set_error)(state, 0,
                          "DDW_operator(): communication setup failed");
    
    alloc_ptr = qx(allocate_eo)(state, &alloc_size, &aligned_ptr, 0, 1, 1);
    if (alloc_ptr == 0)
        return q(set_error)(state, 0, "DDW_operator(): not enough memory");
    
    even_b = aligned_ptr;
    odd_b = qx(step_even)(state, even_b);
    
    BEGIN_TIMING(state);
    qx(op_D)(result->even, &state->even, &state->odd, params, gauge->data,
             source->even, source->odd,
             &flops, &sent, &received, odd_b);
    qx(op_D)(result->odd, &state->odd, &state->even, params, gauge->data,
             source->odd, source->even,
             &flops, &sent, &received, even_b);
    END_TIMING(state, flops, sent, received);
    
    q(free)(state, alloc_ptr, alloc_size);

    return 0;
}       
