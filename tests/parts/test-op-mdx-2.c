#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"
#include "op-routines.h"

char *op_a_name = "lib:conj(DDW)";
char *op_b_name = "conj(DDW)";

double
read_gauge(int dir,
	   const int pos[4],
	   int a, int b,
	   int re_im,
	   void *env)
{
    int i;
    unsigned int v = sum_init(seed_u);

    v = sum_add(v, re_im);
    v = sum_add(v, a);
    v = sum_add(v, b);
    for (i = 0; i < 4; i++) {
	v = sum_add(v, pos[i]);
	v = sum_add(v, dir);
    }
    v = sum_add(v, seed_u);
    return sum_fini(v);
}

static double
read_fermion(unsigned int seed,
	     const int pos[5],
	     int c, int d,
	     int re_im)
{
    int i;
    unsigned int v = sum_init(seed);

    v = sum_add(v, c);
    v = sum_add(v, d);
    v = sum_add(v, re_im);
    for (i = 0; i < 5; i++)
	v = sum_add(v, pos[i]);
    v = sum_add(v, seed);
    return sum_fini(v);
}

double
read_fermion_a(const int pos[5],
	       int c, int d,
	       int re_im,
	       void *env)
{
    return read_fermion(seed_a, pos, c, d, re_im);
}

double
read_fermion_b(const int pos[5],
	       int c, int d,
	       int re_im,
	       void *env)
{
    return read_fermion(seed_b, pos, c, d, re_im);
}

int
operator_a(void)
{
    struct QX(Fermion) *fermion_x;
    double x, y;

    if (QOP_MDWF_allocate_fermion(&fermion_x, state)) {
	zprint("operator_b(): alloc failed on x");
	return 1;
    }

    if (QX(DDW_operator_conjugated)(fermion_x, params, gauge, fermion_a)) {
	zprint("operator_A(): failed");
	return 1;
    };

    dot_fermion(&x, &y, fermion_b, fermion_x);

    QOP_MDWF_free_fermion(&fermion_x);
    zprint("compose: %20.10e %20.10e", x, y);
    return 0;
}

/* same in parts */
static void
op_DDF_even(struct Fermion *res_e,
	    struct Q(State) *state,
	    struct Q(Parameters) *params,
	    const struct SUn *U,
	    struct Fermion *t0_e, struct Fermion *t1_e, struct Fermion *t0_o,
	    const struct Fermion *src_e, const struct Fermion *src_o)
{
    op_Fx_even(t0_e, state, U, src_o);
    op_Bx_even(t1_e, params, t0_e);
    op_Ax_even(t0_e, params, src_e);
    qx(f_add3)(res_e, state->even.full_size, state->even.Ls, t1_e, +1.0, t0_e);
}

static void
op_DDF_odd(struct Fermion *res_o,
	    struct Q(State) *state,
	    struct Q(Parameters) *params,
	    const struct SUn *U,
	    struct Fermion *t0_o, struct Fermion *t1_o, struct Fermion *t0_e,
	    const struct Fermion *src_o, const struct Fermion *src_e)
{
    op_Fx_odd(t0_o, state, U, src_e);
    op_Bx_odd(t1_o, params, t0_o);
    op_Ax_odd(t0_o, params, src_o);
    qx(f_add3)(res_o, state->odd.full_size, state->odd.Ls, t1_o, +1.0, t0_o);
}

int
operator_b(void)
{
    double x, y;
    struct Q(State) *state = gauge->state;
    struct QX(Fermion) *fermion_x;
    struct QX(Fermion) *fermion_y;
    struct QX(Fermion) *fermion_z;

    if (QOP_MDWF_allocate_fermion(&fermion_x, state)) {
	zprint("operator_b(): alloc failed on x");
	return 1;
    }

    if (QOP_MDWF_allocate_fermion(&fermion_y, state)) {
	zprint("operator_b(): alloc failed on y");
	return 1;
    }

    if (QOP_MDWF_allocate_fermion(&fermion_z, state)) {
	zprint("operator_b(): alloc failed on z");
	return 1;
    }

    op_DDF_even(fermion_x->even,
		state, params, gauge->data,
		fermion_y->even, fermion_z->even, fermion_y->odd,
		fermion_a->even, fermion_a->odd);
    op_DDF_odd(fermion_x->odd,
	       state, params, gauge->data,
	       fermion_y->odd, fermion_z->odd, fermion_y->even,
	       fermion_a->odd, fermion_a->even);

    dot_fermion(&x, &y, fermion_b, fermion_x);
    QOP_MDWF_free_fermion(&fermion_z);
    QOP_MDWF_free_fermion(&fermion_y);
    QOP_MDWF_free_fermion(&fermion_x);
    zprint("parts  : %20.10e %20.10e", x, y);
    return 0;
}
