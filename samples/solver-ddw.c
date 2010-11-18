/* Shamir fermions solver for the preconditioned equation.
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
 */
#include <qop-mdwf3.h>
#include <stdlib.h>  /* for NULL */
#include <qmp.h>
#include "common.h"
#define SELF "DDW"

int
main(int argc, char *argv[])
{
    struct QOP_MDWF_State *state = NULL;
    struct QOP_MDWF_Parameters *params = NULL;
    struct QOP_MDWF_Fermion *rhs = NULL;
    struct QOP_MDWF_Fermion *sol = NULL;
    struct QOP_MDWF_Fermion *guess = NULL;
    struct QOP_MDWF_Gauge *gauge = NULL;
	struct QOP_MDWF_Config config;
    int status;
    int out_iterations;
    double out_epsilon;

    /* begin substrate */
    if (init_qmp(argc, argv, SELF " solver example",
		 QOP_MDWF_DEFAULT_PRECISION, NULL, NULL))
	goto end;

    /* start MDWF */
	config.self = QMP_get_node_number();
	config.master_p = QMP_is_primary_node();
    config.rank = 4;
    config.lat = lattice;
    config.ls = lattice[4];
    config.net = network;
    config.neighbor_up = neighbor_up;
    config.neighbor_down = neighbor_down;
    config.sublattice = get_sublattice;
    config.env = NULL;
	
    if (QOP_MDWF_init(&state, &config)) {
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
	QOP_MDWF_import_fermion(&rhs, state, read_fermion, &rhs_seed) ||
	QOP_MDWF_import_fermion(&guess, state, read_fermion, &sol_seed) ||
	QOP_MDWF_allocate_fermion(&sol, state)) {
	goto no_memory;
    }

    /* call the solver */
    status = QOP_MDWF_DDW_CG(sol, &out_iterations, &out_epsilon,
			     params, guess, gauge, rhs,
			     max_iterations, min_epsilon,
			     options);

    zprint(SELF, "status %d, iterations %d, epsilon %e",
	   status, out_iterations, out_epsilon);

    /* In real life one would like to export sol to this side now ... */
    if (status) {
	zprint(SELF, "CG erorr %s", QOP_MDWF_error(state));
    }

    /* Get statistics from the MDWF layer */
    report_performance(state, SELF);

    /* free fields -- it's ok to free a NULL */
no_memory:
    QOP_MDWF_free_fermion(&sol);
    QOP_MDWF_free_fermion(&guess);
    QOP_MDWF_free_fermion(&rhs);
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
