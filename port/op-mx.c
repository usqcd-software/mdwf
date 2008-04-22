#include <mdwf.h>

void
qx(op_even_Mx)(struct Fermion *r_e,
	       struct Q(State) *state,
	       const struct Q(Parameters) *params,
	       const struct SUn *gauge,
	       const struct Fermion *s_e,
	       long long *flops,
	       long long *sent,
	       long long *received,
	       struct Fermion *t_e,
	       struct Fermion *t_o)
{
    qx(op_A1xBx)(t_e, &state->even, params, s_e, flops);
    qx(op_A1xBxFx)(t_o, &state->odd, params, gauge, t_e, flops, sent, received);
    qx(op_1mFx)(r_e, &state->even, gauge, s_e, t_o, flops, sent, received);
}
