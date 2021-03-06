#include <mdwf.h>

/* Solve
 *   (M^\dagget M + s[n]) v_psi[n] = eta
 *     for n = 0 ... count - 1
 *     and
 *   M^\dagger M psi = eta
 *
 *
 */ 

#define MAX_OPTIONS (Q(LOG_CG_RESIDUAL)     | \
                     Q(LOG_TRUE_RESIDUAL)   | \
                     Q(FINAL_CG_RESIDUAL))

int
QX(MxM_SCG)(struct QX(VectorFermion)       *v_psi,
            struct QX(HalfFermion)         *psi,
            int                             *out_iterations,
            double                          *out_epsilon,
            const struct Q(Parameters)      *params,
            const double                     shift[],
            const struct QX(Gauge)          *gauge,
            const struct QX(HalfFermion)    *eta,
            int                              max_iterations,
            double                           min_epsilon,
            unsigned int                     options)
{
    DECLARE_STATE;
    int count = 0;
    long long flops = 0;
    long long sent = 0;
    long long received = 0;
    void *sptr = 0;
    size_t sptr_size = 0;
    double *v = 0;
    double *w = 0;
    double *adn = 0;
    double *bdd = 0;
    void *pptr = 0;
    size_t pptr_size = 0;
    void *temps = 0;
    struct SUn *U = 0;
    int status = 1;
    struct VectorFermion *vpi_e = 0;
    struct Fermion *rho_e = 0;
    struct Fermion *pi_e = 0;
    struct Fermion *zeta_e = 0;
    struct Fermion *t0_e = 0;
    struct Fermion *t1_e = 0;
    struct Fermion *t2_e = 0;
    struct Fermion *t0_o = 0;
    double rhs_norm = 0;

    /* check arguments */
    CHECK_ARG0(v_psi);
    CHECK_ARGn(psi, "MxM_SCG");
    CHECK_ARGn(params, "MxM_SCG");
    CHECK_POINTER(shift, "MxM_SCG");
    CHECK_ARGn(gauge, "MxM_SCG");
    CHECK_ARGn(eta, "MxM_SCG");
    count = v_psi->count;

    /* setup communication */
    if (q(setup_comm)(state, sizeof (REAL))) {
        return q(set_error)(state, 0, "MxM_SCG(): communication setup failed");
    }

    sptr_size = count * 4 * sizeof (double);
    sptr = q(malloc)(state, sptr_size);
    if (sptr == 0) {
        return q(set_error)(state, 0, "MxM_SCG(): not enough memory");
    }
    adn = sptr;
    bdd = adn + count;
    v = bdd + count;
    w = v + count;

    /* allocate locals */
    pptr = qx(allocate_eo)(state, &pptr_size, &temps,
                           0, /* header */
                           6 + count, /* evens */
                           1); /* odds */
    if (pptr == 0) {
        q(free)(state, sptr, sptr_size);
        return q(set_error)(state, 0, "MxM_CG(): not enough memory");
    }
    U = gauge->data;
    t0_e   = temps;
    t1_e   = temps = qx(step_even)(state, temps);
    t2_e   = temps = qx(step_even)(state, temps);
    rho_e  = temps = qx(step_even)(state, temps);
    pi_e   = temps = qx(step_even)(state, temps);
    zeta_e = temps = qx(step_even)(state, temps);
    t0_o   = temps = qx(step_even)(state, temps);
    vpi_e  = temps = qx(step_odd)(state, temps);

    /* clear bits we do not understand */
    options = options & MAX_OPTIONS;

    BEGIN_TIMING(state);
    /* compute the norm of the RHS */
    flops += qx(f_norm)(&rhs_norm, state->even.full_size, state->Ls, eta->even);
    QMP_sum_double(&rhs_norm);
    if (options) {
        qx(zprint)(state, "MxM SCG", "rhs norm %e normalized epsilon %e",
                   rhs_norm, min_epsilon * rhs_norm);
    }
    
    /* solve */
    status = qx(scg_solver)(v_psi->even, psi->even, count, "MxM SCG",
                            out_iterations, out_epsilon,
                            state, params, shift, U, eta->even, 
                            max_iterations, min_epsilon * rhs_norm, options,
                            &flops, &sent, &received,
                            v, w, adn, bdd,
                            rho_e, vpi_e, pi_e, zeta_e,
                            t0_e, t1_e, t2_e, t0_o);

    END_TIMING(state, flops, sent, received);

    /* handle zero mode properly */
    if (status > 1)
        goto end;

    /* output final residuals if desired */
    if (options) {
        qx(zprint)(state, "MxM SCG", "status %d, total iterations %d",
                  status, *out_iterations);
    }
    if (options & (Q(FINAL_CG_RESIDUAL) | Q(LOG_CG_RESIDUAL))) {
        double norm = rhs_norm == 0? 1: rhs_norm;

        qx(zprint)(state, "MxM SCG",
                   "zero shift residual %e normalized %e",
                   *out_epsilon, *out_epsilon / norm);
    }
    if (rhs_norm != 0.0)
        *out_epsilon = *out_epsilon / rhs_norm;

end:
    /* free memory */
    q(free)(state, pptr, pptr_size);
    q(free)(state, sptr, sptr_size);
    if (status != 0) {
        q(set_error)(state, 0, "MxM_SCG() solver failed to converge");
    }
    return status;
}
