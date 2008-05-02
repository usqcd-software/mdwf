#include <mdwf.h>

void
qx(cg_precondition)(struct Fermion *xi0_e,
		    struct Fermion *chi_e,
		    struct Q(State) *state,
		    const struct Q(Parameters) *params,
		    const struct SUn *U,
		    const struct Fermion *psi0_e,
		    const struct Fermion *eta_e,
		    const struct Fermion *eta_o,
		    long long *flops,
		    long long *sent,
		    long long *received,
		    struct Fermion *t0_e,
		    struct Fermion *t1_e,
		    struct Fermion *t0_o)
{
    qx(op_BA1)(t0_o, &state->odd, params, eta_o, flops);
    qx(op_1mF)(t0_e, &state->even, U, eta_e, t0_o, flops, sent, received);
    qx(op_BA1)(t1_e, &state->even, params, t0_e, flops);
    qx(op_even_Mx)(chi_e, state, params, U, t1_e, flops, sent, received,
		   t0_e, t0_o);
    qx(op_B)(xi0_e, &state->even, params, psi0_e, flops);
}
