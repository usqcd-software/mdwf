#include <mdwf.h>

/* Solve
 *   D_dw psi = eta
 *
 * with psi_0 as an initial guess
 *
 */ 

int
Q(deflated_mixed_D_CG)(struct QD(Fermion)          *psi,
                       int                         *out_iterations,
                       double                      *out_epsilon,
                       const struct Q(Parameters)  *params,
                       const struct QD(Fermion)    *psi_0,
                       const struct QD(Gauge)      *gauge,
                       const struct QD(Fermion)    *eta,
                       struct QF(Deflator)         *deflator,
                       int                          f_iter,
                       double                       f_epsilon,
                       int                          max_iterations,
                       double                       min_epsilon,
                       unsigned int                 options)
{
    DECLARE_STATE;

    /* check arguments */
    CHECK_ARG0(psi);
    CHECK_ARGn(params, "eigDDW_CG");
    CHECK_ARGn(psi_0, "eigDDW_CG");
    CHECK_ARGn(gauge, "eigDDW_CG");
    CHECK_ARGn(eta, "eigDDW_CG");
    CHECK_ARGn(deflator, "eigDDW_CG");

    if (deflator->loading)
      return q(set_error)(state, 0, "eigDDW_CG is loading vectors");

    return q(mixed_cg)(state, "eigDDW_CG", params,
                       psi, out_iterations, out_epsilon,
                       psi_0, gauge, eta, deflator,
                       f_iter, f_epsilon, max_iterations, min_epsilon,
                       options);
}

