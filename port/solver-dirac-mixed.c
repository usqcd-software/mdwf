#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>
#undef QOP_MDWF_DEFAULT_PRECISION
#define QOP_MDWF_DEFAULT_PRECISION 'D'
#include <mdwf.h>

/* Solve
 *   D_dw psi = eta
 *
 * with psi_0 as an initial guess
 *
 */ 

#define MAX_OPTIONS (Q(LOG_CG_RESIDUAL)     | \
                     Q(FINAL_DIRAC_RESIDUAL))

static void *
allocate_fgauge(struct Q(State) *state,
                struct QF(Gauge) *fg,
                const struct QD(Gauge) *dg)
{
    void *pgf;
    void *ptr;
    int u_s;
    size_t size;

    fg->state = state;
    fg->size = 0;
    u_s = qf(sizeof_gauge)(state->volume);
    pgf = q(allocate_aligned)(state, &size, &ptr, 0, u_s);
    if (pgf == NULL)
        return NULL;

    fg->size = size;
    fg->data = ptr;
    q(g_f_eq_d)(fg->data, Q(DIM) * state->volume, dg->data);

    return pgf;
}

int
Q(mixed_D_CG)(struct QD(Fermion)          *psi,
              int                         *out_iterations,
              double                      *out_epsilon,
              const struct Q(Parameters)  *params,
              const struct QD(Fermion)    *psi_0,
              const struct QD(Gauge)      *gauge,
              const struct QD(Fermion)    *eta,
              int                          f_iter,
              int                          max_iterations,
              double                       min_epsilon,
              unsigned int                 options)
{
    static const char *comm_failed = "mDDW_CG(): communication setup failed";
    static const char *no_mem = "mDDW_CG(): not enough memory";
    DECLARE_STATE;
    long long flops = 0;
    long long sent = 0;
    long long received = 0;
    void *ptr = 0;
    void *ptrF = 0;
    void *ptr_g = 0;
    size_t ptr_size = 0;
    size_t ptrF_size = 0;
    void *temps = 0;
    int status = 1;
    struct FermionD *t0_e = 0;
    struct FermionD *t1_e = 0;
    struct FermionD *t2_e = 0;
    struct FermionD *t0_o = 0;
    struct FermionD *t1_o = 0;
    struct FermionD *chi_e = 0;
    struct FermionD *xi_e = 0;
    struct FermionF *delta_Fe = 0;
    struct FermionF *xi_Fe = 0;
    struct FermionF *rho_Fe = 0;
    struct FermionF *pi_Fe = 0;
    struct FermionF *zeta_Fe = 0;
    struct FermionF *t0_Fe = 0;
    struct FermionF *t1_Fe = 0;
    struct FermionF *t2_Fe = 0;
    struct FermionF *t0_Fo = 0;
    struct FermionF *t1_Fo = 0;
    struct QF(Gauge) gauge_F;
    double dirac_residual = 0;
    double rhs_norm = 0;
    double scaled_eps;
    int iter_left;
    const char *cg_error = 0;
#define CG_ERROR(msg) do { cg_error = msg; goto end; } while (0)
#define CG_ERROR_T(msg) do { \
        END_TIMING(state, flops, sent, received); \
        cg_error = msg; goto end; } while (0)

    gauge_F.data = 0;
    gauge_F.size = 0;

    /* check arguments */
    CHECK_ARG0(psi);
    CHECK_ARGn(params, "mDDW_CG");
    CHECK_ARGn(psi_0, "mDDW_CG");
    CHECK_ARGn(gauge, "mDDW_CG");
    CHECK_ARGn(eta, "mDDW_CG");

    /* setup communication */
    if (q(setup_comm)(state, sizeof (double)))
        CG_ERROR(comm_failed);

    /* allocate locals */
    ptr = qd(allocate_eo)(state, &ptr_size, &temps,
                          0, /* header */
                          5, /* evens */
                          2); /* odds */
    if (ptr == 0)
        CG_ERROR(no_mem);
    t0_e  = temps;
    t1_e  = temps = qd(step_even)(state, temps);
    t2_e  = temps = qd(step_even)(state, temps);
    chi_e = temps = qd(step_even)(state, temps);
    xi_e  = temps = qd(step_even)(state, temps);
    t0_o  = temps = qd(step_even)(state, temps);
    t1_o  = temps = qd(step_odd)(state, temps);

    ptrF = qf(allocate_eo)(state, &ptrF_size, &temps,
                           0, /* header */
                           8, /* evens */
                           2); /* odds */
    if (ptrF == 0)
        CG_ERROR(no_mem);
    delta_Fe = temps;
    xi_Fe    = temps = qf(step_even)(state, temps);
    rho_Fe   = temps = qf(step_even)(state, temps);
    pi_Fe    = temps = qf(step_even)(state, temps);
    zeta_Fe  = temps = qf(step_even)(state, temps);
    t0_Fe    = temps = qf(step_even)(state, temps);
    t1_Fe    = temps = qf(step_even)(state, temps);
    t2_Fe    = temps = qf(step_even)(state, temps);
    t0_Fo    = temps = qf(step_even)(state, temps);
    t1_Fo    = temps = qf(step_odd)(state, temps);

    ptr_g = allocate_fgauge(state, &gauge_F, gauge);
    if (ptr_g == 0)
        CG_ERROR(no_mem);

    /* clear bits we do not understand */
    options = options & MAX_OPTIONS;

    BEGIN_TIMING(state);
    /* compute t1{eo} = eta - D psi_0 */
    qd(op_D)(t0_e, &state->even, &state->odd, params, gauge->data,
             psi_0->even, psi_0->odd,
             &flops, &sent, &received, t1_o);
    qd(op_D)(t0_o, &state->odd, &state->even, params, gauge->data,
             psi_0->odd, psi_0->even,
             &flops, &sent, &received, t1_e);
    qd(f_add3)(t1_e, state->even.full_size, state->Ls, eta->even, -1, t0_e);
    qd(f_add3)(t1_o, state->odd.full_size, state->Ls, eta->odd, -1, t0_o);

    /* compute chi_e from t1{eo} */
    qd(op_BA1)(t0_o, &state->odd, params, t1_o,
               &flops);
    qd(op_1mF)(t0_e, &state->even, gauge->data, t1_e, t0_o,
               &flops, &sent, &received);
    /* t1_o is eta_o in the notes */
    qd(op_BA1)(t1_e, &state->even, params, t0_e, &flops);
    qd(op_even_Mx)(chi_e, state, params, gauge->data, t1_e,
                   &flops, &sent, &received,
                   t0_e, t0_o);
    qd(f_zero)(xi_e, state->even.full_size, state->Ls);
    /* norm of the normal equation right hand side */
    flops += qd(f_norm)(&rhs_norm, state->even.full_size, state->Ls, chi_e);
    scaled_eps = min_epsilon * rhs_norm;
    if (options)
        qd(zprint)(state, "mDDW CG", "rhs norm %e normalized epsilon %e",
                   rhs_norm, scaled_eps);
    /* solve in MxM xi_e = chi_e in mixed precision with zeta_e(in) = 0 */
    for (iter_left = max_iterations; iter_left > 0 && status; ) {
        int here_iters;
        
        qd(op_even_M)(t0_e, state, params, gauge->data, xi_e,
                      &flops, &sent, &received, t0_o);
        qd(op_even_Mx)(t1_e, state, params, gauge->data, t0_e,
                       &flops, &sent, &received, t2_e, t0_o);
        q(f_f_eq_dmd)(delta_Fe, state->even.full_size, state->Ls, chi_e, t1_e);
        qf(f_zero)(xi_Fe, state->even.full_size, state->Ls);
        
        if (q(setup_comm)(state, sizeof (float)))
            CG_ERROR_T(comm_failed);

        status = qf(cg_solver)(xi_Fe, "mDDW CG", &here_iters, out_epsilon,
                               state, params, gauge_F.data,
                               delta_Fe, NULL, NULL,
                               iter_left > f_iter? f_iter: iter_left,
                               scaled_eps, options,
                               &flops, &sent, &received,
                               rho_Fe, pi_Fe, zeta_Fe,
                               t0_Fe, t1_Fe, t2_Fe, t0_Fo, t1_Fo);
        if (q(setup_comm)(state, sizeof (double)))
            CG_ERROR_T(comm_failed);

        flops += q(f_d_peq_f)(xi_e, state->even.full_size, state->Ls, xi_Fe);

        if (status > 1)
            goto end;
    }

    /* compute psi = psi0 + .... */
    qd(op_B1)(psi->even, &state->even, params, xi_e, &flops);
    qd(op_1mF)(t0_o, &state->odd, gauge->data, t1_o, xi_e,
               &flops, &sent, &received);
    qd(op_A1)(psi->odd, &state->odd, params, t0_o, &flops);
    qd(f_add2)(psi->odd, state->odd.full_size, state->Ls, +1, psi_0->odd);
    qd(f_add2)(psi->even, state->even.full_size, state->Ls, +1, psi_0->even);
    
    if (options & (Q(FINAL_DIRAC_RESIDUAL) | Q(LOG_DIRAC_RESIDUAL))) {
        dirac_residual = qx(cg_dirac_error)(psi->even, psi->odd,
                                            state, params, gauge->data,
                                            eta->even, eta->odd,
                                            &flops, &sent, &received,
                                            t0_e, t1_e, t0_o);
    }

    END_TIMING(state, flops, sent, received);

    /* output final residuals if desired */
    if (options) {
        qx(zprint)(state, "mDDW CG", "status %d, total iterations %d",
                  status, *out_iterations);
    }
    if (options & (Q(FINAL_CG_RESIDUAL) | Q(LOG_CG_RESIDUAL))) {
        double norm = rhs_norm == 0? 1: rhs_norm;

        qx(zprint)(state, "mDDW CG", "solver residual %e normalized %e",
                   *out_epsilon, *out_epsilon / norm);
    }
    if (options & (Q(FINAL_DIRAC_RESIDUAL) | Q(LOG_DIRAC_RESIDUAL))) {
        double norm = rhs_norm == 0? 1: rhs_norm;

        qx(zprint)(state, "mDDW CG", "Dirac residual %e normalized %e",
                   dirac_residual, dirac_residual / norm);
    }
    if (rhs_norm != 0.0)
        *out_epsilon = *out_epsilon / rhs_norm;

end:
    /* free memory */
    if (ptr)
        q(free)(state, ptr, ptr_size);
    if (ptrF)
        q(free)(state, ptrF, ptrF_size);
    if (ptr_g)
        q(free)(state, ptr_g, gauge_F.size);
    if (cg_error)
        q(set_error)(state, 0, cg_error);
    else if (status != 0)
        q(set_error)(state, 0, "mDDW_CG() solver failed to converge");

    return status;
}

