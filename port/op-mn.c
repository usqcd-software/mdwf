#include <mdwf.h>

void
qx(op_even_Mn)(struct Fermion *r_e,
	       double *global_norm,
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
    qx(op_1mBA1F_norm)(r_e, global_norm, &state->even, params,
		       gauge, s_e, t_o, flops, sent, received);
    QMP_sum_double(global_norm);
}
