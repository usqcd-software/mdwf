#include <qop-mdwf3.h>

static int
do_xxx(struct QOP_MDWF_State *state, struct QOP_MDWF_Parameters *params,
       char *name)
{
    struct QOP_MDWF_Gauge *U;
    struct QOP_MDWF_Fermion *rhs;
    struct QOP_MDWF_Fermion *guess;
    struct QOP_MDWF_Fermion *R;
    int max_iter = 100;
    double eps = 1e-4;
    int out_iter;
    double out_eps;

    zprint("%s: starting conjugate gradient test in precision %c",
	   name, QOP_MDWF_DEFAULT_PRECISION);

    if (QOP_MDWF_import_gauge(&U, state, read_gauge, NULL)) {
	zprint("%s: import gauge failed", name);
	goto no_U;
    }
    if (QOP_MDWF_import_fermion(&rhs, state, read_fermion, NULL)) {
	zprint("%s: import rhs failed", name);
	goto no_rhs;
    }
    if (QOP_MDWF_import_fermion(&guess, state, read_fermion, NULL)) {
	zprint("%s: import initial guess failed", name);
	goto no_guess;
    }
    if (QOP_MDWF_allocate_fermion(&R, state)) {
	zprint("%s: allocating solution failed", name);
	goto no_R;
    }

    max_iter = 100;
    eps = 1e-4;
    do {
	zprint("%s: loop max_iter=%d, eps=%g", name, max_iter, eps);
	start_perf();
	QOP_MDWF_DDW_CG(R, &out_iter, &out_eps,
			params, guess, U, rhs, max_iter, eps);
	max_iter *= 2;
	eps *= 0.5;
    } while (report_perf(name, state, out_iter, out_eps) && eps != 0.0);

    QOP_MDWF_free_fermion(&R);
    QOP_MDWF_free_fermion(&guess);
    QOP_MDWF_free_fermion(&rhs);
    QOP_MDWF_free_gauge(&U);
    zprint("%s: end max_iter=%d, esp=%g", name, max_iter/2, eps * 2.0);
    return 0;

no_R:
    QOP_MDWF_free_fermion(&guess);
no_guess:
    QOP_MDWF_free_fermion(&rhs);
no_rhs:
    QOP_MDWF_free_gauge(&U);
no_U:
    return 1;
}

#undef do_xxx
#include "clear-mdwf.h"
