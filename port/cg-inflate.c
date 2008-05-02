#include <mdwf.h>

void
qx(cg_inflate)(struct Fermion *psi_e,
	       struct Fermion *psi_o,
	       struct Q(State) *state,
	       const struct Q(Parameters) *params,
	       const struct SUn *U,
	       const struct Fermion *eta_o,
	       const struct Fermion *xi_e,
	       long long *flops,
	       long long *sent,
	       long long *received,
	       struct Fermion *t_o)
{
    qx(op_B1)(psi_e, &state->even, params, xi_e, flops);
    qx(op_1mF)(t_o, &state->odd, U, eta_o, xi_e, flops, sent, received);
    qx(op_A1)(psi_o, &state->odd, params, t_o, flops);
}

