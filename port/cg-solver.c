#include <mdwf.h>

#if QOP_MDWF_DEFAULT_PRECISION == 'F'
#define DF_PREAMBLE(psi_e, rho_e, r, chi_e) do {                        \
        if (q(df_preamble)(state, deflator, psi_e, rho_e, r, chi_e,     \
                           &ws, options)) {                             \
            q(set_error)(state, 0, "cg_solver() not enough memory");    \
            return CG_NOEMEM;                                           \
        } } while (0)
#define DF_UPDATE0(a1,b1,a0,b0,r,rho)                           \
    q(df_update0)(state, deflator, a1, b1, a0, b0, r, rho, options)
#define DF_UPDATE1(a1,b1,a0,b0,r,rho,A_rho)                             \
    q(df_update1)(state, deflator, a1, b1, a0, b0, r, rho, A_rho, options)
#define DF_POSTAMBLE() \
    do { q(df_postamble)(state, deflator, &ws, options); } while (0)
#else
#define DF_PREAMBLE(psi_e, rho_e, r, chi_e) do {        \
        qx(f_zero)(psi_e, Ls, e_size);                  \
        qx(f_copy)(rho_e, Ls, e_size, chi_e);           \
        qx(f_norm)(r, Ls, e_size, rho_e);               \
        QMP_sum_double(r);                              \
    } while (0)
#define DF_UPDATE0(a1,b1,a0,b0,r,rho)  0
#define DF_UPDATE1(a1,b1,a0,b0,r,rho,A_rho)  0
#define DF_POSTAMBLE()  do {} while (0)
#endif

void
qx(cg_operator)(struct Fermion           *res_e,
                const struct Fermion     *psi_e,
                struct MxM_workspace     *ws)
{
    qx(op_even_M)(ws->tmp_e, ws->state, ws->params, ws->gauge, psi_e,
                  ws->flops, ws->sent, ws->received,
                  ws->tmp_o);
    qx(op_even_Mx)(res_e, ws->state, ws->params, ws->gauge, ws->tmp_e,
                   ws->flops, ws->sent, ws->received,
                   ws->tmp2_e, ws->tmp_o);
}

CG_STATUS
qx(cg_solver)(struct Fermion              *xi_e,
              const char                  *source,
              int                         *out_iter,
              double                      *out_epsilon,
              struct Q(State)             *state,
              const struct Q(Parameters)  *params,
              const struct SUn            *U,
              const struct Fermion        *chi_e,
              struct Q(Deflator)          *deflator,
              int                          max_iter,
              double                       epsilon,
              unsigned                     options,
              long long                   *flops,
              long long                   *sent,
              long long                   *received,
              struct Fermion              *rho_e,
              struct Fermion              *pi_e,
              struct Fermion              *zeta_e,
              struct Fermion              *t0_e,
              struct Fermion              *t1_e,
              struct Fermion              *t2_e,
              struct Fermion              *t0_o)
{
#if QOP_MDWF_DEFAULT_PRECISION == 'F'
    double a0 = 1, b0 = 0;
#endif /*  QOP_MDWF_DEFAULT_PRECISION == 'F' */
    int df_status;
    int e_size = state->even.full_size;
    int Ls = state->Ls;
    double a, b, g, r, norm_omega;
    int i;
    struct MxM_workspace  ws;

    ws.state     = state;
    ws.params    = params;
    ws.gauge     = U;
    ws.tmp_e     = t0_e;
    ws.tmp2_e    = t2_e;
    ws.tmp_o     = t0_o;
    ws.flops     = flops;
    ws.sent      = sent;
    ws.received  = received;

    DF_PREAMBLE(xi_e, rho_e, &r, (struct Fermion *) chi_e);
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
        if (norm_omega == 0.0) {
            *out_iter = i;
            *out_epsilon = r;
            q(set_error)(state, 0, "cg_solver() hit zero mode");
            DF_POSTAMBLE();
            return CG_ZEROMODE;
        }
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
        df_status = DF_UPDATE0(a, b, a0, b0, g, rho_e);
        if (-1 == df_status) {
            qx(cg_operator)(zeta_e, rho_e, &ws);
            df_status = DF_UPDATE1(a, b, a0, b0, g, rho_e, zeta_e);
        }
        if (3 == df_status) {
            *out_iter = i;
            *out_epsilon = r;
            DF_POSTAMBLE();
            return CG_EIGCONV;
        }
#if QOP_MDWF_DEFAULT_PRECISION == 'F'
        a0 = a;
        b0 = b;
#endif /*  QOP_MDWF_DEFAULT_PRECISION == 'F' */
        if (options)
            qx(cg_log)(r, source,
                       i, xi_e, state, params, U, chi_e,
                       flops, sent, received,
                       options,
                       t0_e, t1_e, t2_e, t0_o);
    }
end:
    *out_iter = i;
    *out_epsilon = r;
    DF_POSTAMBLE();
    if (i == max_iter)
        return CG_MAXITER;
    return CG_SUCCESS;
}
