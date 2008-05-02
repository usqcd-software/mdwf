/* Follow conventions of the mdwf.pdf */
#include <stdarg.h>
static void
zprint(struct Q(State) *state,
       const char *fmt,
       ...)
{
    va_list va;
    char buffer[4096];

    if (state->master_p) {
	va_start(va, fmt);
	vsnprintf(buffer, sizeof (buffer) - 1, fmt, va);
	va_end(va);
	QMP_printf("MDWF: %s\n", buffer);
    }
	
}


static void
cg_precondition(struct Fermion *xi0_e,
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

static void
cg_inflate(struct Fermion *psi_e,
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
cg_true_residual(const struct Fermion *xi_e,
		 struct Q(State) *state,
		 const struct Q(Parameters) *params,
		 const struct SUn *U,
		 const struct Fermion *chi_e,
		 long long *flops,
		 long long *sent,
		 long long *received,
		 struct Fermion *t0_e,
		 struct Fermion *t1_e,
		 struct Fermion *t2_e,
		 struct Fermion *t0_o)
{
    double norm;

    qx(op_even_M)(t0_e, state, params, U, xi_e, flops, sent, received,
		  t0_o);
    qx(op_even_Mx)(t1_e, state, params, U, t0_e, flops, sent, received,
		   t2_e, t0_o);
    *flops += qx(f_diff_norm)(&norm, state->even.full_size, state->Ls,
			      t1_e, chi_e);
    QMP_sum_double(&norm);

    return norm;
}

static double
cg_dirac_error(const struct Fermion *psi_e,
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

static double
cg_dirac_residual(const struct Fermion *xi_e,
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
		  struct Fermion *t2_e,
		  struct Fermion *t0_o,
		  struct Fermion *t1_o)
{
    cg_inflate(t0_e, t0_o, state, params, U, eta_o, xi_e,
	       flops, sent, received,
	       t1_o);
    return cg_dirac_error(t0_e, t0_o, state, params, U, eta_e, eta_o,
			  flops, sent, received,
			  t1_e, t2_e, t1_o);
}

static void
cg_log(double cg_res, int iter,
       const struct Fermion *xi_e,
       struct Q(State) *state,
       const struct Q(Parameters) *params,
       const struct SUn *U,
       const struct Fermion *chi_e,
       const struct Fermion *eta_e,
       const struct Fermion *eta_o,
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
	true_res = cg_true_residual(xi_e, state, params, U, chi_e,
				    flops, sent, received,
				    t0_e, t1_e, t2_e, t0_o);
    }
    if (options & Q(LOG_DIRAC_RESIDUAL)) {
	dirac_res = cg_dirac_residual(xi_e, state, params, U, eta_e, eta_o,
				      flops, sent, received,
				      t0_e, t1_e, t2_e, t0_o, t1_o);
    }
#define ITER_LOG (Q(LOG_CG_RESIDUAL) |	  \
                  Q(LOG_TRUE_RESIDUAL) |  \
                  Q(LOG_DIRAC_RESIDUAL))
    switch (options & ITER_LOG) {
#undef ITER_LOG
    default:
	break;
    case Q(LOG_CG_RESIDUAL):
	zprint(state,
	       "CG step %5d"
	       "  CG residual %11.4e",
	       iter, cg_res);
	break;
    case Q(LOG_TRUE_RESIDUAL):
	zprint(state,
	       "CG step %5d"
	       "  true residual %11.4e",
	       iter, true_res);
	break;
    case Q(LOG_CG_RESIDUAL) | Q(LOG_TRUE_RESIDUAL):
	zprint(state, 
	       "CG step %5d"
	       "  CG residual %11.4e"
	       "  true residual %11.4e",
	       iter, cg_res, true_res);
	break;
    case Q(LOG_DIRAC_RESIDUAL):
	zprint(state,
	       "CG step %5d"
	       "  Dirac residual %11.4e",
	       iter, dirac_res);
	break;
    case Q(LOG_DIRAC_RESIDUAL) | Q(LOG_CG_RESIDUAL):
	zprint(state,
	       "CG step %5d"
	       "  Dirac residual %11.4e"
	       "  CG residual %11.4e",
	       iter, dirac_res, cg_res);
	break;
    case Q(LOG_DIRAC_RESIDUAL) | Q(LOG_TRUE_RESIDUAL):
	zprint(state,
	       "CG step %5d"
	       "  Dirac residual %11.4e"
	       "  true residual %11.4e",
	       iter, dirac_res, true_res);
	break;
    case Q(LOG_DIRAC_RESIDUAL) | Q(LOG_CG_RESIDUAL) | Q(LOG_TRUE_RESIDUAL):
	zprint(state,
	       "CG step %5d"
	       "  Dirac residual %11.4e"
	       "  CG residual %11.4e"
	       "  true residual %11.4e",
	       iter, dirac_res, cg_res, true_res);
	break;
    }
}

static int
cg_solver(struct Fermion *xi_e,
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
	    cg_log(r, i, xi_e, state, params, U, chi_e, eta_e, eta_o,
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

	      
int
cg(struct QX(Fermion)                *psi,       /* result                  */
   int                               *out_iters, /* number of cg steps used */
   double                            *out_eps,   /* final precision         */
   const struct QOP_MDWF_Parameters  *params,    /* MDWF parameters         */
   const struct QX(Fermion)          *psi0,      /* initial guess for psi   */
   const struct QX(Gauge)            *gauge,     /* gauge field             */
   const struct QX(Fermion)          *eta,       /* Dirac eq. rhs           */
   int                                max_iters, /* maximum iterations      */
   double                             epsilon,   /* desired precision       */
   unsigned int                       options)   /* log options             */
{
    DECLARE_STATE;
    long long flops = 0;
    long long sent = 0;
    long long received = 0;
    void *ptr;
    size_t ptr_size;
    void *temps;
    struct SUn *U;
    int status;
    struct Fermion *xi_e;
    struct Fermion *chi_e;
    struct Fermion *rho_e;
    struct Fermion *pi_e;
    struct Fermion *zeta_e;
    /* XXX other pointers to temps */
    struct Fermion *t0_e;
    struct Fermion *t1_e;
    struct Fermion *t2_e;
    struct Fermion *t0_o;
    struct Fermion *t1_o;
    double dirac_residual = 0;
    double rhs_norm = 0;

    /* check arguments */
    CHECK_ARG0(psi);
    CHECK_ARGn(params, "DDW_CG");
    CHECK_ARGn(psi0, "DDW_CG");
    CHECK_ARGn(gauge, "DDW_CG");
    CHECK_ARGn(eta, "DDW_CG");

    /* setup communication */
    if (q(setup_comm)(state, sizeof (REAL))) {
	return q(set_error)(state, 0, "MDWF_CG(): communication setup failed");
    }

    /* allocate locals */
    ptr = q(allocate_eo)(state, &ptr_size, &temps,
			 0, /* header */
			 8, /* evens */
			 2, /* odds */
			 sizeof (REAL));
    if (ptr == 0) {
	return q(set_error)(state, 0, "MDWF_CG(): not enough memory");
    }
    U = gauge->data;
    t0_e  = temps;
    t1_e  = temps = q(step_even)(state, temps, sizeof (REAL));
    t2_e  = temps = q(step_even)(state, temps, sizeof (REAL));
    xi_e  = temps = q(step_even)(state, temps, sizeof (REAL));
    chi_e = temps = q(step_even)(state, temps, sizeof (REAL));
    rho_e = temps = q(step_even)(state, temps, sizeof (REAL));
    pi_e = temps = q(step_even)(state, temps, sizeof (REAL));
    zeta_e = temps = q(step_even)(state, temps, sizeof (REAL));
    t0_o  = temps = q(step_even)(state, temps, sizeof (REAL));
    t1_o  = temps = q(step_odd)(state, temps, sizeof (REAL));

    BEGIN_TIMING(state);
    /* precondition */
    cg_precondition(xi_e, chi_e, state, params,
		    U, psi0->even, eta->even, eta->odd,
		    &flops, &sent, &received,
		    t0_e, t1_e, t0_o);
    /* solve */
    status = cg_solver(xi_e, out_iters, out_eps,
		       state, params, U,
		       chi_e, eta->even, eta->odd,
		       max_iters, epsilon, options,
		       &flops, &sent, &received,
		       rho_e, pi_e, zeta_e,
		       t0_e, t1_e, t2_e, t0_o, t1_o);
    
    /* unprecondition */
    cg_inflate(psi->even, psi->odd,
	       state, params, U, eta->odd, xi_e,
	       &flops, &sent, &received,
	       t0_o);
    if (options & (Q(FINAL_DIRAC_RESIDUAL) | Q(LOG_DIRAC_RESIDUAL))) {
	dirac_residual = cg_dirac_error(psi->even, psi->odd,
					state, params, U,
					eta->even, eta->odd,
					&flops, &sent, &received,
					t0_e, t1_e, t0_o);
	flops += qx(op_norm2)(&rhs_norm, eta, state);
    }
    END_TIMING(state, flops, sent, received);

    zprint(state, "CG results: status %d, iterations %d, residual %e",
	   status, *out_iters, *out_eps);
    
    if (options & (Q(FINAL_DIRAC_RESIDUAL) | Q(LOG_DIRAC_RESIDUAL))) {
	zprint(state, "CG Dirac residual: %e rhs norm %e",
	       dirac_residual, rhs_norm);
    }

    q(free)(state, ptr, ptr_size);
    return status;
}
