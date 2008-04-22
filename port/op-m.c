#include <mdwf.h>

void
qx(op_even_M)(struct Fermion *r_e,
	      struct Q(State) *state,
	      const struct Q(Parameters) *params,
	      const struct SUn *gauge,
	      const struct Fermion *s_e,
	      long long *flops,
	      long long *sent,
	      long long *received,
	      struct Fermion *t_o)
{
    qx(op_BA1F)(t_o, &state->odd, params,
		gauge, s_e, flops, sent, received);
    qx(op_1mBA1F)(r_e, &state->even, params,
		  gauge, s_e, t_o, flops, sent, received);
}
