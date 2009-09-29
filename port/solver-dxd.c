#include <mdwf.h>

/* Solve
 *   M^\dagget M psi = eta
 *
 * with psi as an initial guess
 *
 */ 

#define MAX_OPTIONS (Q(LOG_CG_RESIDUAL)     | \
                     Q(LOG_TRUE_RESIDUAL)   | \
                     Q(FINAL_CG_RESIDUAL))

static void
qx(op_d)(struct QX(Fermion)          *result,
         struct Q(State)             *state,
         const struct Q(Parameters)  *params,
         const struct QX(Gauge)      *gauge,
         const struct QX(Fermion)    *source,
         long long                   *flops,
         long long                   *sent,
         long long                   *received,
         struct Fermion              *tmp_e,
         struct Fermion              *tmp_o)
{
    qx(op_D)(result->even, &state->even, &state->odd, params, gauge->data,
             source->even, source->odd,
             flops, sent, received, tmp_o);
    qx(op_D)(result->odd, &state->odd, &state->even, params, gauge->data,
             source->odd, source->even,
             flops, sent, received, tmp_e);
}

static void
qx(op_dx)(struct QX(Fermion)          *result,
          struct Q(State)             *state,
          const struct Q(Parameters)  *params,
          const struct QX(Gauge)      *gauge,
          const struct QX(Fermion)    *source,
          long long                   *flops,
          long long                   *sent,
          long long                   *received)
{
    qx(op_AxpBxFx)(result->even, &state->even, params,
                   gauge->data, source->even, source->odd,
                   flops, sent, received);
    qx(op_AxpBxFx)(result->odd, &state->odd, params,
                   gauge->data, source->odd, source->even,
                   flops, sent, received);
}

static void
qx(op_dn)(struct QX(Fermion)          *result,
          double                      *norm,
          struct Q(State)             *state,
          const struct Q(Parameters)  *params,
          const struct QX(Gauge)      *gauge,
          const struct QX(Fermion)    *source,
          long long                   *flops,
          long long                   *sent,
          long long                   *received,
          struct Fermion              *tmp_e,
          struct Fermion              *tmp_o)      
{
    double n_e, n_o;

    qx(op_D_norm)(result->even, &n_e,
                  &state->even, &state->odd, params, gauge->data,
                  source->even, source->odd,
                  flops, sent, received, tmp_o);
    qx(op_D_norm)(result->odd, &n_o,
                  &state->odd, &state->even, params, gauge->data,
                  source->odd, source->even,
                  flops, sent, received, tmp_e);
    *norm = n_e + n_o;
    QMP_sum_double(norm);
}

static void
qx(F_copy)(struct QX(Fermion)      *result,
           struct Q(State)         *state,
           const struct Q(Fermion) *source)
{
    int Ls = state->Ls;

    qx(f_copy)(result->even, state->even.full_size, Ls, source->even);
    qx(f_copy)(result->odd, state->odd.full_size, Ls, source->odd);
}

static int
qx(F_norm2)(const struct QX(Fermion) *fermion,
            double                   *n_f,
            struct Q(State)          *state)
{
    int Ls = state->Ls;
    int flops = 0;
    double n_e, n_o;

    flops += qx(f_norm)(&n_e, state->even.full_size, Ls, fermion->even);
    flops += qx(f_norm)(&n_o, state->odd.full_size, Ls, fermion->odd);
    *n_f = n_e + n_o;
    QMP_sum_double(n_f);

    return flops;
}

static int
qx(F_madd3)(struct QX(Fermion)        *r,
            struct Q(State)           *state,
            const struct QX(Fermion)  *a,
            double                     s,
            const struct QX(Fermion)  *b)
{
    int flops = 0;
    int Ls = state->Ls;

    flops += qx(f_add3)(r->even, state->even.full_size, Ls, a->even, s, b->even);
    flops += qx(f_add3)(r->odd, state->odd.full_size, Ls, a->odd, s, b->odd);

    return flops;
}

static int
qx(F_madd2)(struct QX(Fermion)        *r,
            struct Q(State)           *state,
            double                     s,
            const struct QX(Fermion)  *b)
{
    int flops = 0;
    int Ls = state->Ls;

    flops += qx(f_add2)(r->even, state->even.full_size, Ls, s, b->even);
    flops += qx(f_add2)(r->odd, state->odd.full_size, Ls, s, b->odd);
    
    return flops;
}

static int
qx(F_madd2_norm)(struct QX(Fermion)         *r,
                 double                     *norm,
                 struct Q(State)            *state,
                 double                     s,
                 const struct QX(Fermion)   *b)
{
    int flops = 0;
    int Ls = state->Ls;
    double n_e, n_o;

    flops += qx(f_add2_norm)(r->even, &n_e, state->even.full_size, Ls, s, b->even);
    flops += qx(f_add2_norm)(r->odd, &n_o, state->odd.full_size, Ls, s, b->odd);
    *norm = n_e + n_o;
    QMP_sum_double(norm);

    return flops;
}

static int
qx(F_xmadd)(struct QX(Fermion) *xi,
            struct QX(Fermion) *pi,
            struct Q(State)    *state,
            double a,
            double b,
            const struct QX(Fermion) *rho)
{
    int flops = 0;
    int Ls = state->Ls;

    flops += qx(cg_xp)(xi->even, pi->even, state->even.full_size, Ls, a, b, rho->even);
    flops += qx(cg_xp)(xi->odd, pi->odd, state->odd.full_size, Ls, a, b, rho->odd);

    return flops;
}

int
QX(DxD_CG)(struct QX(Fermion)            *psi,
           int                           *out_iterations,
           double                        *out_epsilon,
           const struct Q(Parameters)    *params,
           const struct QX(Fermion)      *psi_0,
           const struct QX(Gauge)        *gauge,
           const struct QX(Fermion)      *eta,
           int                            max_iterations,
           double                         min_epsilon,
           unsigned int                   options)
{
    DECLARE_STATE;
    long long flops = 0;
    long long sent = 0;
    long long received = 0;
    void *ptr = 0;
    size_t ptr_size = 0;
    void *temps = 0;
    struct SUn *U = 0;
    int status = 1;
    struct QX(Fermion) *rho = 0;
    struct QX(Fermion) *pi = 0;
    struct QX(Fermion) *omega = 0;
    struct QX(Fermion) *zeta = 0;
    struct Fermion *tmp_e = 0;
    struct Fermion *tmp_o = 0;
    double rhs_norm = 0;
    static char *no_mem = "DxD_CG(): not enough memory";
    char *e_msg = 0;
    int k;
    double r, g, a, b, omega_norm;
#define CG_ERROR(msg) do { e_msg = msg; goto end; } while (0)

    /* check arguments */
    CHECK_ARG0(psi);
    CHECK_ARGn(psi_0, "DxD_CG");
    CHECK_ARGn(params, "DxD_CG");
    CHECK_ARGn(gauge, "DxD_CG");
    CHECK_ARGn(eta, "DxD_CG");

    /* setup communication */
    if (q(setup_comm)(state, sizeof (REAL))) {
        return q(set_error)(state, 0, "DxD_CG(): communication setup failed");
    }

    /* allocate locals */
    ptr = qx(allocate_eo)(state, &ptr_size, &temps,
                          0, /* header */
                          1, /* evens */
                          1); /* odds */
    if (ptr == 0)
        CG_ERROR(no_mem);

    tmp_e  = temps;
    tmp_o  = temps = qx(step_even)(state, temps);

    if (QX(allocate_fermion)(&rho, state) ||
        QX(allocate_fermion)(&pi, state) ||
        QX(allocate_fermion)(&rho, state) ||
        QX(allocate_fermion)(&zeta, state))
        CG_ERROR(no_mem);

    U = gauge->data;
    /* clear bits we do not understand */
    options = options & MAX_OPTIONS;

    BEGIN_TIMING(state);
    /* compute the norm of the RHS */
    flops += qx(f_norm)(&rhs_norm, state->even.full_size, state->Ls, eta->even);
    if (options) {
        qx(zprint)(state, "DxD CG", "rhs norm %e normalized epsilon %e",
                   rhs_norm, min_epsilon * rhs_norm);
    }

    qx(F_copy)(psi, state, psi_0);
    qx(op_d)(pi, state, params, gauge, psi,
             &flops, &sent, &received, tmp_e, tmp_o);
    qx(op_dx)(omega, state, params, gauge, pi,
              &flops, &sent, &received);
    flops += qx(F_madd3)(rho, state, eta, -1.0, omega);
    qx(F_copy)(pi, state, rho);
    flops += qx(F_norm2)(rho, &r, state);
    for (k = 0; k < max_iterations; k++) {
        qx(op_dn)(omega, &omega_norm, state, params, gauge, pi,
                  &flops, &sent, &received, tmp_e, tmp_o);
        qx(op_dx)(zeta, state, params, gauge, omega,
                  &flops, &sent, &received);
        if (omega_norm == 0) {
            *out_iterations = k;
            *out_epsilon = r;
            status = 2;
            CG_ERROR("DxD_CG(): hit zero mode");
        }
        a = r / omega_norm;
        flops += qx(F_madd2_norm)(rho, &g, state, -a, zeta);
        if (g < rhs_norm * min_epsilon) {
            flops += qx(F_madd2)(psi, state, a, pi);
            r = g;
            status = 0;
            break;
        }
        b = g / r;
        r = g;
        flops += qx(F_xmadd)(psi, pi, state, a, b, rho);
    }
    *out_iterations = k;
    *out_epsilon = r;
    END_TIMING(state, flops, sent, received);

    /* output final residuals if desired */
    if (options) {
        qx(zprint)(state, "DxD CG", "status %d, total iterations %d",
                  status, *out_iterations);
    }
    if (options & (Q(FINAL_CG_RESIDUAL) | Q(LOG_CG_RESIDUAL))) {
        double norm = rhs_norm == 0? 1: rhs_norm;

        qx(zprint)(state, "DxD CG", "solver residual %e normalized %e",
                   *out_epsilon, *out_epsilon / norm);
    }
    if (rhs_norm != 0.0)
        *out_epsilon = *out_epsilon / rhs_norm;

end:
    /* free memory */
    q(free)(state, ptr, ptr_size);
    QX(free_fermion)(&rho);
    QX(free_fermion)(&pi);
    QX(free_fermion)(&omega);
    QX(free_fermion)(&zeta);
    if (status != 0) {
        q(set_error)(state, 0, e_msg);
    }
    return status;
}
