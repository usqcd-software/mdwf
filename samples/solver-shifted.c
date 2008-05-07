/* Shamir fermions shifted solver for the preconditioned equation.
 * arguments:
 *   Nx Ny Nz Nt
 *   Lx Ly Lz Lt Ls
 *   M
 *   m_5
 *   kappa
 *   guess-seed
 *   src-seed
 *   U-seed
 *   U-scale
 *   max-iterations
 *   min-epsilon
 *   verbose?
 *   shift ...
 */
#include <qop-mdwf3.h>
#include <stdlib.h>  /* for NULL */
#include "common.h"
#define SELF "MxM"

void
check_solution(const struct QOP_MDWF_HalfFermion *rhs,
	       const struct QOP_MDWF_HalfFermion *sol,
	       const struct QOP_MDWF_Parameters *params,
	       const struct QOP_MDWF_Gauge *gauge,
	       struct QOP_MDWF_HalfFermion *t1,
	       struct QOP_MDWF_HalfFermion *t2,
	       double shift,
	       int index)
{
    double norm;

    QOP_MDWF_M_operator(t1, params, gauge, sol);
    QOP_MDWF_M_operator_conjugated(t2, params, gauge, t1);
    QOP_MDWF_madd_half_fermion(t1, t2, shift, sol);
    QOP_MDWF_madd_half_fermion(t2, t1, -1.0, rhs);
    QOP_MDWF_norm2_half_fermion(&norm, t2);
    if (index < 0) {
	zprint(SELF, "base solution epsilon %e", norm);
    } else {
	zprint(SELF, "solution at shift[%d]=%8.6f epsilon %e",
	       index, shift, norm);
    }
}

int
main(int argc, char *argv[])
{
    struct QOP_MDWF_State *state = NULL;
    struct QOP_MDWF_Parameters *params = NULL;
    struct QOP_MDWF_HalfFermion *rhs = NULL;
    struct QOP_MDWF_HalfFermion *sol = NULL;
    struct QOP_MDWF_HalfFermion *t1 = NULL;
    struct QOP_MDWF_HalfFermion *t2 = NULL;
    struct QOP_MDWF_HalfFermion *ssol = NULL;
    struct QOP_MDWF_VectorFermion *vsol = NULL;
    struct QOP_MDWF_Gauge *gauge = NULL;
    int status;
    int out_iterations;
    double out_epsilon;
    double *shift;
    int count;
    int i;

    /* begin substrate */
    if (init_qmp(argc, argv, SELF " shifted solver example",
		 QOP_MDWF_DEFAULT_PRECISION, &count, &shift))
	goto end;

    /* start MDWF */
    if (QOP_MDWF_init(&state, lattice, network, this_node, primary_p,
		      get_sublattice, NULL)) {
	zprint("MDWF failed to init: %s", QOP_MDWF_error(state));
	goto end;
    }

    /* create parameters */
    if (QOP_MDWF_set_Shamir(&params, state, kappa, M, m_5)) {
	zprint("set_Shamir() failed: %s", QOP_MDWF_error(state));
	goto end_no_params;
    }

    /* Initialize fields */
    if (QOP_MDWF_import_gauge(&gauge, state, read_gauge, &U_seed) ||
	QOP_MDWF_import_half_fermion(&rhs, state, read_fermion, &rhs_seed) ||
	QOP_MDWF_allocate_half_fermion(&sol, state) ||
	QOP_MDWF_allocate_half_fermion(&t1, state) ||
	QOP_MDWF_allocate_half_fermion(&t2, state) ||
	QOP_MDWF_allocate_half_fermion(&ssol, state) ||
	QOP_MDWF_allocate_vector_fermion(&vsol, state, count)) {
	goto no_memory;
    }

    /* call the solver */
    status = QOP_MDWF_MxM_SCG(vsol, ssol, &out_iterations, &out_epsilon,
			      params, shift, gauge, rhs,
			      max_iterations, min_epsilon,
			      options);

    zprint(SELF, "status %d, iterations %d, epsilon %e",
	   status, out_iterations, out_epsilon);

    /* In real life one would like to export sol to this side now ... */
    check_solution(rhs, ssol, params, gauge, t1, t2, 0.0, -1);
    for (i = 0; i < count; i++) {
	QOP_MDWF_get_vector_fermion(sol, vsol, i);
	check_solution(rhs, sol, params, gauge, t1, t2, shift[i], i);
    }

    /* Get statistics from the MDWF layer */
    report_performance(state, SELF);

    /* free fields -- it's ok to free a NULL */
no_memory:
    QOP_MDWF_free_vector_fermion(&vsol);
    QOP_MDWF_free_half_fermion(&t1);
    QOP_MDWF_free_half_fermion(&t2);
    QOP_MDWF_free_half_fermion(&ssol);
    QOP_MDWF_free_half_fermion(&sol);
    QOP_MDWF_free_half_fermion(&rhs);
    QOP_MDWF_free_gauge(&gauge);

    /* free parameters */
    QOP_MDWF_free_parameters(&params);

end_no_params:
    /* end MDWF */
    QOP_MDWF_fini(&state);
    
    /* end substrate */
end:
    fini_qmp();

    return 0;
}
