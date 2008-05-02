#include <mdwf.h>

int
qx(cg_solver)(struct Fermion *xi_e, const char *source,
	      int *out_iter,
	      double *out_epsilon,
	      struct Q(State) *state,
	      const struct Q(Parameters) *params,
	      const struct SUn *U,
	      const struct Fermion *chi_e,
	      const struct Fermion *eta_e,
	      const struct Fermion *eta_o,
	      int max_iter,
	      double epsilon,
	      unsigned options,
	      long long *flops,
	      long long *sent,
	      long long *received,
	      struct Fermion *rho_e,
	      struct Fermion *pi_e,
	      struct Fermion *zeta_e,
	      struct Fermion *t0_e,
	      struct Fermion *t1_e,
	      struct Fermion *t2_e,
	      struct Fermion *t0_o,
	      struct Fermion *t1_o)
{
    int e_size = state->even.full_size;
    int Ls = state->Ls;
    double a, b, g, r, norm_omega;
    int i;

    qx(op_even_M)(t1_e, state, params, U, xi_e, flops, sent, received,
		  t0_o);
    qx(op_even_Mx)(t0_e, state, params, U, t1_e, flops, sent, received,
		   t2_e, t0_o);
    *flops += qx(f_add3)(rho_e, e_size, Ls, chi_e, -1.0, t0_e);
    *flops += qx(f_norm)(&r, e_size, Ls, rho_e);
    QMP_sum_double(&r);
    qx(f_copy)(pi_e, e_size, Ls, rho_e);
    if (r < epsilon) {
	i = 0;
	goto end;
    }
    for (i = 0; i < max_iter; i++) {
	qx(op_even_Mn)(t0_e, &norm_omega, state, params, U, pi_e,
		       flops, sent, received,
		       t0_o);
	qx(op_even_Mx)(zeta_e, state, params, U, t0_e,
		       flops, sent, received,
		       t1_e, t0_o);
	/* could z_n == 0? -- if so, we are stuck on a zero mode */
	a = r / norm_omega;
	*flops += qx(f_add2_norm)(rho_e, &g, e_size, Ls, -a, zeta_e);
	QMP_sum_double(&g);
	if (g < epsilon) {
	    *flops += qx(f_add2)(xi_e, e_size, Ls, a, pi_e);
	    r = g;
	    break;
	}
	b = g / r;
	r = g;
	qx(cg_xp)(xi_e, pi_e, e_size, Ls, a, b, rho_e);
	if (options)
	    qx(cg_log)(r, source,
		       i, xi_e, state, params, U, chi_e, eta_e, eta_o,
		       flops, sent, received,
		       options,
		       t0_e, t1_e, t2_e, t0_o, t1_o);
    }
end:
    *out_iter = i;
    *out_epsilon = r;
    if (i == max_iter)
	return 1;
    return 0;
}
