#include <stdarg.h>
static void
qx(cg_show)(struct Q(State) *state,
	    const char *fmt,
	    ...)
{
    va_list va;
    char buffer[4096];

    if (state->master_p) {
	va_start(va, fmt);
	vsnprintf(buffer, sizeof (buffer) - 1, fmt, va);
	va_end(va);
	QMP_printf("MDWF CG: %s\n", buffer);
    }
	
}


static void
qx(cg_project)(struct Fermion *psi_e,
	       struct Fermion *rhs_e,
	       struct Q(State) *state,
	       const struct Q(Parameters) *params,
	       const struct SUn *U,
	       const struct Fermion *guess_e,
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
    qx(op_even_Mx)(rhs_e, state, params, U, t1_e, flops, sent, received,
		   t0_e, t0_o);
    qx(op_B)(psi_e, &state->even, params, guess_e, flops);
}

static void
qx(cg_unproject)(struct Fermion *psi_e,
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

static double
qx(cg_true_residual)(const struct Fermion *psi_e,
		     struct Q(State) *state,
		     const struct Q(Parameters) *params,
		     const struct SUn *U,
		     const struct Fermion *rhs_e,
		     long long *flops,
		     long long *sent,
		     long long *received,
		     struct Fermion *t0_e,
		     struct Fermion *t1_e,
		     struct Fermion *t2_e,
		     struct Fermion *t0_o)
{
    double norm;

    qx(op_even_M)(t0_e, state, params, U, psi_e, flops, sent, received,
		  t0_o);
    qx(op_even_Mx)(t1_e, state, params, U, t0_e, flops, sent, received,
		   t2_e, t0_o);
    *flops += qx(f_diff_norm)(&norm, state->even.full_size, state->Ls,
			      t1_e, rhs_e);
    QMP_sum_double(&norm);

    return norm;
}

static double
qx(cg_dirac_residual)(const struct Fermion *psi_e,
		      struct Q(State) *state,
		      const struct Q(Parameters) *params,
		      const struct SUn *U,
		      const struct Fermion *rhs_e,
		      const struct Fermion *rhs_o,
		      long long *flops,
		      long long *sent,
		      long long *received,
		      struct Fermion *t0_e,
		      struct Fermion *t1_e,
		      struct Fermion *t2_e,
		      struct Fermion *t0_o,
		      struct Fermion *t1_o)
{
    double e_norm, o_norm, norm;

    qx(cg_unproject)(t0_e, t0_o, state, params, U, rhs_o, psi_e,
		     flops, sent, received,
		     t1_o);
    qx(op_D)(t1_e, &state->even, &state->odd, params, U, t0_e, t0_o,
	     flops, sent, received,
	     t1_o);
    qx(op_D)(t1_o, &state->odd, &state->odd, params, U, t0_e, t0_o,
	     flops, sent, received,
	     t2_e);
    *flops += qx(f_diff_norm)(&e_norm, state->even.full_size, state->Ls,
			      t1_e, rhs_e);
    *flops += qx(f_diff_norm)(&o_norm, state->odd.full_size, state->Ls,
			      t1_o, rhs_o);
    
    norm = e_norm + o_norm;
    QMP_sum_double(&norm);

    return norm;
}

static void
qx(cg_log)(double cg_res, int iter,
	   const struct Fermion *psi_e,
	   struct Q(State) *state,
	   const struct Q(Parameters) *params,
	   const struct SUn *U,
	   const struct Fermion *rhs_e,
	   const struct Fermion *xi_e,
	   const struct Fermion *xi_o,
	   long long *flops,
	   long long *sent,
	   long long *received,
	   unsigned int options,
	   struct Fermion *t0_e,
	   struct Fermion *t1_e,
	   struct Fermion *t2_e,
	   struct Fermion *t0_o,
	   struct Fermion *t1_o)
{
    double true_res = 0.0;
    double dirac_res = 0.0;

    if (options & Q(LOG_TRUE_RESIDUAL)) {
	true_res = qx(cg_true_residual)(psi_e, state, params, U, rhs_e,
					flops, sent, received,
					t0_e, t1_e, t2_e, t0_o);
    }
    if (options & Q(LOG_DIRAC_RESIDUAL)) {
	dirac_res = qx(cg_dirac_residual)(psi_e, state, params, U, xi_e, xi_o,
					  flops, sent, received,
					  t0_e, t1_e, t2_e, t0_o, t1_o);
    }
    switch (options & Q(LOG_EVERYTHING)) {
    default:
	break;
    case Q(LOG_CG_RESIDUAL):
	qx(cg_show)(state,
		    "cg step %5d"
		    "  CG residual %20.10e",
		    iter, cg_res);
	break;
    case Q(LOG_TRUE_RESIDUAL):
	qx(cg_show)(state,
		    "cg step %5d"
		    "  true residual %20.10e",
		    iter, true_res);
	break;
    case Q(LOG_CG_RESIDUAL) | Q(LOG_TRUE_RESIDUAL):
	qx(cg_show)(state, 
		    "cg step %5d"
		    "  CG residual %20.10e"
		    "  true residual %20.10e",
		    iter, cg_res, true_res);
	break;
    case Q(LOG_DIRAC_RESIDUAL):
	qx(cg_show)(state,
		    "cg step %5d"
		    "  Dirac residual %20.10e",
		    iter, dirac_res);
	break;
    case Q(LOG_DIRAC_RESIDUAL) | Q(LOG_CG_RESIDUAL):
	qx(cg_show)(state,
		    "cg step %5d"
		    "  Dirac residual %20.10e"
		    "  CG residual %20.10e",
		    iter, dirac_res, cg_res);
	break;
    case Q(LOG_DIRAC_RESIDUAL) | Q(LOG_TRUE_RESIDUAL):
	qx(cg_show)(state,
		    "cg step %5d"
		    "  Dirac residual %20.10e"
		    "  true residual %20.10e",
		    iter, dirac_res, true_res);
	break;
    case Q(LOG_DIRAC_RESIDUAL) | Q(LOG_CG_RESIDUAL) | Q(LOG_TRUE_RESIDUAL):
	qx(cg_show)(state,
		    "cg step %5d"
		    "  Dirac residual %20.10e"
		    "  CG residual %20.10e"
		    "  true residual %20.10e",
		    iter, dirac_res, cg_res, true_res);
	break;
    }
}

static int
qx(cg_solver)(struct Fermion *psi_e,
	      int *out_iter,
	      double *out_epsilon,
	      struct Q(State) *state,
	      const struct Q(Parameters) *params,
	      const struct SUn *U,
	      const struct Fermion *rhs_e,
	      const struct Fermion *xi_e,
	      const struct Fermion *xi_o,
	      int max_iter,
	      double epsilon,
	      unsigned options,
	      long long *flops,
	      long long *sent,
	      long long *received,
	      struct Fermion *r_e,
	      struct Fermion *p_e,
	      struct Fermion *q_e,
	      struct Fermion *t0_e,
	      struct Fermion *t1_e,
	      struct Fermion *t2_e,
	      struct Fermion *t0_o,
	      struct Fermion *t1_o)
{
    int e_size = state->even.full_size;
    int Ls = state->Ls;
    double alpha, beta, gamma, rho, z_n;
    int i;

    qx(op_even_M)(p_e, state, params, U, psi_e, flops, sent, received,
		  t0_o);
    qx(op_even_Mx)(t0_e, state, params, U, p_e, flops, sent, received,
		   t1_e, t0_o);
    *flops += qx(f_add3)(r_e, e_size, Ls, rhs_e, -1.0, t0_e);
    *flops += qx(f_norm)(&rho, e_size, Ls, r_e);
    QMP_sum_double(&rho);
    qx(f_copy)(p_e, e_size, Ls, r_e);
    if (rho < epsilon) {
	i = 0;
	goto end;
    }
    for (i = 0; i < max_iter; i++) {
	qx(op_even_Mn)(t0_e, &z_n, state, params, U, p_e, flops, sent, received,
		       t0_o);
	qx(op_even_Mx)(q_e, state, params, U, t0_e, flops, sent, received,
		       t1_e, t0_o);
	/* could z_n == 0? -- if so, we are stuck */
	alpha = rho / z_n;
	*flops += qx(f_add2_norm)(r_e, &gamma, e_size, Ls, -alpha, q_e);
	QMP_sum_double(&gamma);
	if (gamma < epsilon) {
	    *flops += qx(f_add2)(psi_e, e_size, Ls, alpha, p_e);
	    rho = gamma;
	    break;
	}
	beta = gamma / rho;
	rho = gamma;
	qx(cg_xp)(psi_e, p_e, e_size, Ls, alpha, beta, r_e);
	if (options)
	    qx(cg_log)(gamma, i, psi_e, state, params, U, rhs_e, xi_e, xi_o,
		       flops, sent, received,
		       options,
		       t0_e, t1_e, t2_e, t0_o, t1_o);
    }
end:
    *out_iter = i;
    *out_epsilon = rho;
    if (i == max_iter)
	return 1;
    return 0;
}

	      

int
cg(struct QX(Fermion) *result,
   int *out_iterations,
   double *out_epsilon,
   const struct QOP_MDWF_Parameters *params,
   const struct QX(Fermion) *guess,
   const struct QX(Gauge) *gauge,
   const struct QX(Fermion) *rhs,
   int max_iteration,
   double epsilon,
   unsigned int options)
{
    /* XXX */
    *out_iterations = 0;
    *out_epsilon = -1;
    return 1;
}
