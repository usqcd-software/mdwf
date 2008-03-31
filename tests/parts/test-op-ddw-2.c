#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "optest.h"
#include "op-routines.h"

char *op_name = "DDW in parts";

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

double
read_fermion(const int pos[5],
	     int c, int d,
	     int re_im,
	     void *env)
{
    int i;
    for (i = 0; i < 5; i++) {
	if (pos[i] != fermion_pos[i])
	    return 0.0;
    }
    if (c != fermion_color ||
	d != fermion_dirac ||
	re_im != fermion_reim)
	return 0.0;
    return 1.0;
}

void
write_fermion(const int pos[5],
	      int c, int d,
	      int re_im,
	      double value,
	      void *env)
{
    if (value != 0.0)
	xprint("   fermion[%d,%d,%d,%d;%d][%d,%d].%d = %12.8f",
	       pos[0], pos[1], pos[2], pos[3], pos[4],
	       c, d, re_im, value);
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
    op_B_odd(t0_o, params, src_o);
    op_F_even(t0_e, state, U, t0_o);
    op_A_even(t1_e, params, src_e);
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
    op_B_even(t0_e, params, src_e);
    op_F_odd(t0_o, state, U, t0_e);
    op_A_odd(t1_o, params, src_o);
    qx(f_add3)(res_o, state->odd.full_size, state->odd.Ls, t1_o, +1.0, t0_o);
}

int
operator(void)
{
    struct Q(State) *state = gauge->state;
    struct QX(Fermion) *fermion_y;
    struct QX(Fermion) *fermion_z;

    if (QOP_MDWF_allocate_fermion(&fermion_y, state)) {
	zprint("operator_b(): alloc failed on y");
	return 1;
    }

    if (QOP_MDWF_allocate_fermion(&fermion_z, state)) {
	zprint("operator_b(): alloc failed on z");
	return 1;
    }

    op_DDF_even(result->even,
		state, params, gauge->data,
		fermion_y->even, fermion_z->even, fermion_y->odd,
		fermion->even, fermion->odd);
    op_DDF_odd(result->odd,
	       state, params, gauge->data,
	       fermion_y->odd, fermion_z->odd, fermion_y->even,
	       fermion->odd, fermion->even);

    QOP_MDWF_free_fermion(&fermion_z);
    QOP_MDWF_free_fermion(&fermion_y);

    return 0;
}
