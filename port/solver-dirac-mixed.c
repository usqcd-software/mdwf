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

int
Q(mixed_DDW_CG)(struct QD(Fermion)          *psi,
                int                         *out_iterations,
                double                      *out_epsilon,
                const struct Q(Parameters)  *params,
                const struct QD(Fermion)    *psi_0,
                const struct QD(Gauge)      *gauge,
                const struct QD(Fermion)    *eta,
                int                          f_iter,
                double                       f_epsilon,
                int                          max_iterations,
                double                       min_epsilon,
                unsigned int                 options)
{
    DECLARE_STATE;

    /* check arguments */
    CHECK_ARG0(psi);
    CHECK_ARGn(params, "mDDW_CG");
    CHECK_ARGn(psi_0, "mDDW_CG");
    CHECK_ARGn(gauge, "mDDW_CG");
    CHECK_ARGn(eta, "mDDW_CG");

    return q(mixed_cg)(state, "mDDW_CG", params,
                       psi, out_iterations, out_epsilon,
                       psi_0, gauge, eta, NULL,
                       f_iter, f_epsilon, max_iterations, min_epsilon,
                       options);
}

