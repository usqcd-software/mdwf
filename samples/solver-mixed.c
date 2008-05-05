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
#include "common.h"
#define SELF "DDW mixed"

int
main(int argc, char *argv[])
{
    struct QOP_MDWF_State *state = NULL;
    struct QOP_MDWF_Parameters *params = NULL;
    struct QOP_F3_MDWF_Fermion *f_rhs = NULL;
    struct QOP_F3_MDWF_Fermion *f_sol = NULL;
    struct QOP_F3_MDWF_Fermion *f_guess = NULL;
    struct QOP_F3_MDWF_Gauge *f_gauge = NULL;
    struct QOP_D3_MDWF_Fermion *d_rhs = NULL;
    struct QOP_D3_MDWF_Fermion *d_sol = NULL;
    struct QOP_D3_MDWF_Fermion *d_guess = NULL;
    struct QOP_D3_MDWF_Gauge *d_gauge = NULL;
    int status;
    int out_iterations;
    double out_epsilon;

    /* begin substrate */
    if (init_qmp(argc, argv, SELF "  solver example", 'm'))
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

    /* single precision start */
    /* Initialize fields */
    if (QOP_F3_MDWF_import_gauge(&f_gauge, state, read_gauge, &U_seed) ||
	QOP_F3_MDWF_import_fermion(&f_rhs, state, read_fermion, &rhs_seed) ||
	QOP_F3_MDWF_import_fermion(&f_guess, state, read_fermion, &sol_seed) ||
	QOP_F3_MDWF_allocate_fermion(&f_sol, state)) {
	goto no_memory;
    }

    /* call the solver */
    status = QOP_F3_MDWF_DDW_CG(f_sol, &out_iterations, &out_epsilon,
				params, f_guess, f_gauge, f_rhs,
				max_iterations, min_epsilon,
				options);

    zprint(SELF, "(signle) status %d, iterations %d, epsilon %e",
	   status, out_iterations, out_epsilon);

    /* In real life one would like to export sol to this side now ... */

    /* Get statistics from the MDWF layer */
    report_performance(state, SELF " (single)");

    /* export the solution */
    if (QOP_F3_MDWF_export_fermion(write_fermion, NULL, f_sol)) {
	zprint(SELF, "error exporting single solution");
	goto no_memory;
    }

    /* how long it takes to export a fermion? */
    report_time(state, SELF " (single) export_fermion()");

    /* release all parallel data */
    QOP_F3_MDWF_free_fermion(&f_sol);
    QOP_F3_MDWF_free_fermion(&f_guess);
    QOP_F3_MDWF_free_fermion(&f_rhs);
    QOP_F3_MDWF_free_gauge(&f_gauge);

    /* double precision polish */
    /* Initialize fields */
    if (QOP_D3_MDWF_import_gauge(&d_gauge, state, read_gauge, &U_seed) ||
	QOP_D3_MDWF_import_fermion(&d_rhs, state, read_fermion, &rhs_seed) ||
	QOP_D3_MDWF_import_fermion(&d_guess, state, read_fermion, &sol_seed) ||
	QOP_D3_MDWF_allocate_fermion(&d_sol, state)) {
	goto no_memory;
    }

    /* call the solver */
    status = QOP_D3_MDWF_DDW_CG(d_sol, &out_iterations, &out_epsilon,
				params, d_guess, d_gauge, d_rhs,
				max_iterations, min_epsilon,
				options);

    zprint(SELF, "(double) status %d, iterations %d, epsilon %e",
	   status, out_iterations, out_epsilon);

    /* In real life one would like to export sol to this side now ... */

    /* Get statistics from the MDWF layer */
    report_performance(state, SELF " (double)");

    /* export the solution */
    if (QOP_D3_MDWF_export_fermion(write_fermion, NULL, d_sol)) {
	zprint(SELF, "error exporting double solution");
	goto no_memory;
    }

    /* how long it takes to export a fermion? */
    report_time(state, SELF " (double) export_fermion()");

    /* release all parallel data */
    QOP_D3_MDWF_free_fermion(&d_sol);
    QOP_D3_MDWF_free_fermion(&d_guess);
    QOP_D3_MDWF_free_fermion(&d_rhs);
    QOP_D3_MDWF_free_gauge(&d_gauge);

no_memory:

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
