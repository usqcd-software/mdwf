#include <mdwf.h>

double
qx(cg_dirac_error)(const struct Fermion *psi_e,
		   const struct Fermion *psi_o,
		   struct Q(State) *state,
		   const struct Q(Parameters) *params,
		   const struct SUn *U,
		   const struct Fermion *eta_e,
		   const struct Fermion *eta_o,
		   long long *flops,
		   long long *sent,
		   long long *received,
		   struct Fermion *t0_e,
		   struct Fermion *t1_e,
		   struct Fermion *t0_o)
{
    double e_norm, o_norm, norm;
    
    qx(op_D)(t0_e, &state->even, &state->odd, params, U, psi_e, psi_o,
	     flops, sent, received,
	     t0_o);
    qx(op_D)(t0_o, &state->odd, &state->even, params, U, psi_o, psi_e,
	     flops, sent, received,
	     t1_e);
    *flops += qx(f_diff_norm)(&e_norm, state->even.full_size, state->Ls,
			      t0_e, eta_e);
    *flops += qx(f_diff_norm)(&o_norm, state->odd.full_size, state->Ls,
			      t0_o, eta_o);
    
    norm = e_norm + o_norm;
    *flops += 1; /* every flop counts ... */
    QMP_sum_double(&norm);

    return norm;
}

