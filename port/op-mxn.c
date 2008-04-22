#include <mdwf.h>

void
qx(op_even_Mxn)(struct Fermion *r_e,
		double *global_norm,
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
    *global_norm = 0.0;
    qx(op_A1xBx)(t_e, &state->even, params, s_e, flops);
    qx(op_A1xBxFx)(t_o, &state->odd, params, gauge, t_e, flops, sent, received);
    qx(op_1mFx_norm)(r_e, global_norm, &state->even, gauge, s_e, t_o,
		     flops, sent, received);
    QMP_sum_double(global_norm);
}
