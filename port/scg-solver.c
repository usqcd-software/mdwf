#include <mdwf.h>

int
qx(scg_solver)(struct VectorFermion *v_xi_e,
               struct Fermion *xi_e,
               int count,
               const char *source,
               int *out_iterations,
               double *out_epsilon,
               struct Q(State) *state,
               const struct Q(Parameters) *params,
               const double shift[],
               const struct SUn *U,
               const struct Fermion *chi_e,
               int max_iterations,
               double min_epsilon,
               unsigned options,
               long long *flops,
               long long *sent,
               long long *received,
               double dp[],
               double d[],
               double dn[],
               double ad[],
               double bdd[],
               struct Fermion *rho_e,
               struct VectorFermion *v_pi_e,
               struct Fermion *pi_e,
               struct Fermion *zeta_e,
               struct Fermion *t0_e,
               struct Fermion *t1_e,
               struct Fermion *t2_e,
               struct Fermion *t0_o)
{
    int e_size = state->even.full_size;
    int Ls = state->Ls;
    int k, j;
    double r, z, ap, bp, a, g, b;


    /* setup */
    qx(f_zero)(xi_e, e_size, Ls);
    qx(f_copy)(rho_e, e_size, Ls, chi_e);
    qx(f_copy)(pi_e, e_size, Ls, chi_e);
    qx(fv_zero)(v_xi_e, e_size, Ls, count);
    qx(fv_copy)(v_pi_e, e_size, Ls, count, rho_e);
    for (j = 0; j < count; j++) {
        dp[j] = d[j] = 1;
    }
    *flops += qx(f_norm)(&r, e_size, Ls, rho_e);
    ap = bp = 1;
    if (r < min_epsilon) {
        k = 0;
        goto end;
    }
    /* loop for convergence */
    for (k = 0; k < max_iterations; k++) {
        qx(op_even_Mn)(t0_e, &z, state, params, U, pi_e,
                       flops, sent, received,
                       t0_o);
        qx(op_even_Mx)(zeta_e, state, params, U, t0_e,
                       flops, sent, received,
                       t1_e, t0_o);
        if (z == 0.0) {
            *out_iterations = k;
            *out_epsilon = r;
            q(set_error)(state, 0, "scg_solver() hit zero mode");
            return 2;
        }
        a = r / z;
        *flops += qx(f_add2_norm)(rho_e, &g, e_size, Ls, -a, zeta_e);
        QMP_sum_double(&g);
        b = g / r;
        r = g;
        for (j = 0; j < count; j++) {
            dn[j] = d[j] * (1.0 + a * shift[j]) + a * bp * (d[j] - dp[j]) / ap;
            ad[j] = a / dn[j];
            bdd[j] = b * d[j] / dn[j];
        }
        if (g < min_epsilon) {
            qx(scg_madd)(xi_e, v_xi_e, e_size, Ls, count, a, ad, pi_e);
            break;
        }
        qx(scg_xp)(xi_e, pi_e,
                   v_xi_e, v_pi_e,
                   e_size, Ls, count,
                   a, b,
                   ad, bdd,
                   rho_e);
        for (j = 0; j < count; j++) {
            dp[j] = d[j];
            d[j] = dn[j];
        }
        bp = b;
        ap = a;
        if (options)
            qx(cg_log)(r, "MxM SCG",
                       k, xi_e, state, params, U, chi_e,
                       flops, sent, received,
                       options,
                       t0_e, t1_e, t2_e, t0_o);
    }
end:
    *out_iterations = k;
    *out_epsilon = r;
    if (k == max_iterations)
        return 1;
    return 0;
}
